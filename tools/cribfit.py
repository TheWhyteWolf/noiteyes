#!/usr/bin/env python3
"""Step H: crib tester for the lossless short-step gak41.

Model (round 5 leader): c_t = D[a_t], a_t = a_(t-1) + b_(t-1) * k(p_t) mod 83,
b_t = b_(t-1) * m(p_t); k: plaintext letter -> step, INJECTIVE (decryption
needs it), in 1..K (short steps); m = 1 except on twist letters. Every message
starts from the same state (b, a) = (1, a0); position 0 is an off-chain
header. With pi = D^-1 a crib (plaintext p_1..p_L on a shared opening) gives
EXACT linear equations over GF(83):

    pi(c_t) - pi(c_(t-1)) - b_(t-1) * k(p_t) = 0,   pi(c_0) := X (start state)

with b_(t-1) known once the crib's twist letter tau and multiplier m are
fixed (tau = none, or one crib letter; m in the order-41 subgroup). These are
solved TOGETHER with Step F's isomorph equations pi(b) = beta_j pi(a) + alpha_j
(tools/affinefit.py blocks, excursions cut; beta_j searched per block). A
branch dies when the solution space forces two symbols to one pi value or two
crib letters to one step. At each surviving leaf the remaining freedom (after
the additive gauge) is enumerated when it is <= MAX_ENUM_DIM dimensional, and
the crib is ACCEPTED only if some point has every step in 1..K, all steps
distinct, and pi injective on the touched symbols.

Gates (synthetic, tournament5.make_lossy, known plaintext and pi; the block
set is thinned to the real corpus's density):
  true   -- the true opening as crib must be accepted, with pi recovered
  wrong  -- shuffled / shifted / random cribs: the false-accept rate is
            measured (this is the tester's power)

Usage (run from repo root; heavy searches belong on enter):
    python3 tools/cribfit.py gate
    python3 tools/cribfit.py real A "CRIBTEXT" [B "CRIBTEXT"]
    python3 tools/cribfit.py file cribs.txt      (lines: A|B<TAB>text)
"""
import itertools
import json
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from affinefit import Echelon, QR, blocks_from  # noqa: E402
from pi_solver import P, load_msgs  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
K = 27
MAX_ENUM_DIM = 3
PRUNE_DIM = 2                                 # range-check internal nodes up to this image dim
NODE_CAP = 200_000
FAMILIES = {'A': ('east-1', 25), 'B': ('east-4', 21)}     # shared [1, end)


def norm_crib(text, spaces):
    t = ''.join(ch for ch in text.upper() if ch.isalpha() or ch == ' ')
    return t.replace(' ', '_') if spaces else t.replace(' ', '')


class CribProblem:
    """Columns: pi[0..82], X (start state), alpha per block, k per crib letter."""

    def __init__(self, blocks, cribs):
        self.blocks = [sorted(set(b)) for b in blocks]
        self.cribs = cribs                      # list of (symbols c_1..c_L, text)
        self.letters = sorted({ch for _, t in cribs for ch in t})
        self.X = P
        self.a0 = P + 1
        self.k0 = self.a0 + len(self.blocks)
        self.ncols = self.k0 + len(self.letters)
        self.kcol = {ch: self.k0 + i for i, ch in enumerate(self.letters)}
        touched = set()
        for syms, _ in cribs:
            touched |= set(syms)
        self.crib_syms = touched
        order, rest = [], set(range(len(self.blocks)))
        while rest:                             # most overlap with what is known first
            j = max(rest, key=lambda j: (len(self.bl_syms(j) & touched), len(self.blocks[j]), -j))
            order.append(j); touched |= self.bl_syms(j); rest.discard(j)
        self.order = order

    def bl_syms(self, j):
        return {x for pr in self.blocks[j] for x in pr}

    def crib_rows(self, tau, m):
        rows = []
        for syms, text in self.cribs:
            b = 1
            for i, (c, ch) in enumerate(zip(syms, text)):
                r = np.zeros(self.ncols, dtype=np.int64)
                r[c] += 1
                if i == 0:
                    r[self.X] -= 1
                else:
                    r[syms[i - 1]] -= 1
                r[self.kcol[ch]] -= b
                rows.append(r % P)
                if ch == tau:
                    b = b * m % P
        return rows

    def block_rows(self, j, beta):
        out = []
        for a, b in self.blocks[j]:
            r = np.zeros(self.ncols, dtype=np.int64)
            r[b] += 1; r[a] -= beta; r[self.a0 + j] -= 1
            out.append(r % P)
        return out

    def distinct_ok(self, ech, syms):
        B = ech.null_basis()
        s = B[sorted(syms)]
        if len({row.tobytes() for row in s}) < len(syms):
            return False
        kk = B[[self.kcol[ch] for ch in self.letters]]
        return len({row.tobytes() for row in kk}) == len(self.letters)

    def leaf(self, ech, syms, max_dim=MAX_ENUM_DIM):
        """Enumerate the remaining freedom on the coordinates that matter (crib
        steps + touched symbols, additive gauge fixed): (n_points, example) or
        (None, dim) when the image is more than max_dim-dimensional. A point
        needs every step in 1..K, all steps distinct, pi injective on syms."""
        B = ech.null_basis()
        kidx = [self.kcol[ch] for ch in self.letters]
        sidx = sorted(syms)
        ref = min(self.crib_syms)
        M = (B[kidx + sidx] - B[ref]) % P       # gauge: pi[ref] = 0
        img = Echelon(M.shape[0])
        for i in range(M.shape[1]):
            img.add(M[:, i])
        Gi = img.R.T if img.R.shape[0] else np.zeros((M.shape[0], 0), dtype=np.int64)
        d = Gi.shape[1]
        if d > max_dim:
            return None, d
        nk = len(kidx)
        grids = np.array(list(itertools.product(range(P), repeat=d)), dtype=np.int64).T \
            if d else np.zeros((0, 1), dtype=np.int64)
        found, ex = 0, None
        for lo in range(0, grids.shape[1], 50_000):
            C = grids[:, lo:lo + 50_000]
            V = (Gi @ C) % P                    # rows: coordinates, cols: points
            kv = V[:nk]
            ok = (kv >= 1).all(0) & (kv <= K).all(0)
            if not ok.any():
                continue
            ks = np.sort(kv[:, ok], axis=0)
            dk = (np.diff(ks, axis=0) != 0).all(0) if nk > 1 else np.ones(ks.shape[1], bool)
            idx = np.nonzero(ok)[0][dk]
            for j in idx:
                sv = V[nk:, j]
                if len(set(sv.tolist())) == len(sv):
                    found += 1
                    if ex is None:
                        ex = {s_: int(v) for s_, v in zip(sidx, sv)}
                        ex['_k'] = {ch: int(v) for ch, v in zip(self.letters, V[:nk, j])}
        return found, ex

    def solve(self, tau, m, node_cap=NODE_CAP, max_leaves=20):
        ech = Echelon(self.ncols)
        for r in self.crib_rows(tau, m):
            ech.add(r)
        if not self.distinct_ok(ech, self.crib_syms):
            return {'status': 'reject-crib-only', 'nodes': 0, 'leaves': []}
        nodes = [0]; leaves = []

        def rec(depth, ech, syms, assign):
            if nodes[0] >= node_cap or len(leaves) >= max_leaves:
                return
            if depth == len(self.order):
                n, info = self.leaf(ech, syms)
                leaves.append({'assign': dict(assign), 'points': n,
                               'info': info if n else None})
                return
            j = self.order[depth]
            s2 = syms | self.bl_syms(j)
            for beta in QR:
                nodes[0] += 1
                e2 = ech.copy()
                for r in self.block_rows(j, beta):
                    e2.add(r)
                if not self.distinct_ok(e2, s2):
                    continue
                n, _ = self.leaf(e2, s2, max_dim=PRUNE_DIM)
                if n == 0:                      # range-infeasible: prune
                    continue
                assign[j] = beta
                rec(depth + 1, e2, s2, assign)
                del assign[j]

        rec(0, ech, set(self.crib_syms), {})
        acc = [lf for lf in leaves if lf['points'] is None or lf['points'] > 0]
        status = ('accept' if any(lf['points'] for lf in acc) else
                  'undetermined' if acc else
                  'capped' if nodes[0] >= node_cap else 'reject')
        return {'status': status, 'nodes': nodes[0], 'n_leaves': len(leaves),
                'n_accept_leaves': len(acc), 'leaves': acc[:3]}


def test_crib(blocks, cribs, twist=True, node_cap=NODE_CAP):
    """cribs: list of (symbols, text). Tries tau = none, then each crib letter
    (most frequent first) x each multiplier; stops at the first ACCEPT.
    Returns (status, accepted/undetermined/capped hypotheses, n_tried).
    status: accept > undetermined > capped > reject."""
    prob = CribProblem(blocks, cribs)
    hyps = [(None, 1)]
    if twist:
        freq = {ch: sum(t.count(ch) for _, t in cribs) for ch in prob.letters}
        hyps += [(t, m) for t in sorted(prob.letters, key=lambda c: -freq[c])
                 for m in QR if m != 1]
    out, n = [], 0
    for tau, m in hyps:
        n += 1
        r = prob.solve(tau, m, node_cap=node_cap)
        if r['status'] in ('accept', 'undetermined', 'capped'):
            out.append({'tau': tau, 'm': m, 'status': r['status'], 'nodes': r['nodes'],
                        'leaves': r.get('leaves', [])})
            if r['status'] == 'accept':
                break
    best = ('accept' if any(o['status'] == 'accept' for o in out) else
            'undetermined' if any(o['status'] == 'undetermined' for o in out) else
            'capped' if out else 'reject')
    return best, out, n


# ------------------------------------------------------------ synthetic gate

def synth(seed, ntw=2):
    """tournament.corpus with the plaintext kept (same RNG consumption)."""
    from tournament import PA, WEIGHTS, dense_under
    from tournament5 import make_lossy
    enc = make_lossy(seed, ntw=ntw)
    rng = random.Random(1000 + seed)
    L = 120
    T = rng.choices(PA, weights=WEIGHTS, k=L)
    names = [f'sim-{i}' for i in range(9)]
    plain = {}
    for nm in names:
        s = list(T); s[0] = rng.choice(PA)
        for _ in range(3):
            p = rng.randrange(12, L - 8)
            s[p:p + 7] = rng.choices(PA, weights=WEIGHTS, k=7)
        plain[nm] = s
    phs = [dense_under(enc, rng, n) for n in (30, 28, 26)]

    def ins(nm, ph, pos):
        if ph:
            plain[nm][pos:pos + len(ph)] = list(ph)
    ins('sim-0', phs[0], 30); ins('sim-1', phs[0], 61); ins('sim-5', phs[0], 34)
    ins('sim-2', phs[1], 40); ins('sim-4', phs[1], 67); ins('sim-3', phs[1], 44)
    ins('sim-6', phs[2], 46); ins('sim-8', phs[2], 60); ins('sim-7', phs[2], 30)
    S = {nm: enc(plain[nm]) for nm in names}
    pit = [0] * P
    for i, c in enumerate(enc.D):
        pit[c] = i
    return S, names, plain, pit, enc


def thin(blocks, target_pairs, rng):
    """Random subset of blocks with about the real corpus's aligned-pair count."""
    idx = list(range(len(blocks))); rng.shuffle(idx)
    out, n = [], 0
    for i in idx:
        if n >= target_pairs:
            break
        out.append(blocks[i]); n += len(set(blocks[i]))
    return out


def gate_jobs(seeds, crib_len=20, n_wrong=6):
    """(seed, kind, symbols, text, blocks) per crib: the true opening plus
    shuffled / shifted / random wrong cribs of the same length."""
    from tournament import PA, WEIGHTS
    from tournament import N9
    real_blocks, _ = blocks_from(load_msgs(), N9)
    target = len({pr for b in real_blocks for pr in b})
    jobs = []
    for sd in seeds:
        S, names, plain, pit, enc = synth(sd)
        rng = random.Random(sd)
        blocks = thin(blocks_from(S, names)[0], target, rng)
        nm = names[0]
        L0 = min(next(i for i in range(1, len(S[nm])) if S[nm][i] != S[n][i]) for n in names[1:4])
        L = min(crib_len, L0 - 1)
        syms = S[nm][1:1 + L]
        true_text = ''.join(plain[nm][1:1 + L])
        jobs.append((sd, 'true', syms, true_text, blocks))
        for w in range(n_wrong):
            kind = w % 3
            if kind == 0:
                t = list(true_text); rng.shuffle(t); t = ''.join(t)
            elif kind == 1:
                t = ''.join(plain[nm][1 + 3 + w:1 + 3 + w + L])
            else:
                t = ''.join(rng.choices(PA, weights=WEIGHTS, k=L))
            jobs.append((sd, ['shuffled', 'shifted', 'random'][kind], syms, t, blocks))
    return jobs


def gate_one(job):
    sd, kind, syms, text, blocks = job
    t0 = time.time()
    status, hyps, n = test_crib(blocks, [(syms, text)], node_cap=GATE_NODE_CAP)
    rec = {'seed': sd, 'kind': kind, 'len': len(syms), 'text': text, 'status': status,
           'hyps_tried': n, 'hyps_kept': len(hyps), 's': round(time.time() - t0, 1)}
    acc = [h for h in hyps if h['status'] == 'accept']
    if acc:
        rec['accepted_as'] = {'tau': acc[0]['tau'], 'm': acc[0]['m']}
        if kind == 'true':
            _, _, _, pit, enc = synth(sd)
            info = next(lf['info'] for lf in acc[0]['leaves'] if lf['info'])
            ks = [x for x in info if x != '_k']
            ref = min(ks)
            best = max(sum(1 for x in ks if (u * (pit[x] - pit[ref])) % P == info[x])
                       for u in range(1, P))
            rec['pi_agree'] = f'{best}/{len(ks)}'
    return rec


GATE_NODE_CAP = 20_000


def gate(seeds=range(1, 13)):
    from multiprocessing import Pool
    jobs = gate_jobs(seeds)
    rep = []
    with Pool(os.cpu_count()) as pool:
        for rec in pool.imap_unordered(gate_one, jobs):
            rep.append(rec)
            print(rec, flush=True)
    summ = {}
    for kind in ('true', 'shuffled', 'shifted', 'random'):
        rs = [r['status'] for r in rep if r['kind'] == kind]
        summ[kind] = {st: rs.count(st) for st in ('accept', 'undetermined', 'capped', 'reject')}
    print('SUMMARY', summ)
    return {'summary': summ, 'rows': sorted(rep, key=lambda r: (r['seed'], r['kind']))}


def FAMILIES_ALL():
    from tournament import N9
    return N9


def real_cribs(args):
    S = load_msgs()
    cribs = []
    for fam, text in zip(args[0::2], args[1::2]):
        nm, end = FAMILIES[fam]
        t = text
        syms = S[nm][1:1 + min(len(t), end - 1)]
        cribs.append((syms, t[:len(syms)]))
    return cribs


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'gate'
    if mode == 'gate':
        rep = gate()
        json.dump(rep, open(os.path.join(ROOT, 'data', 'cribfit_gate.json'), 'w'), indent=1)
        return
    from tournament import N9
    blocks, _ = blocks_from(load_msgs(), N9)
    if mode == 'real':
        cribs = real_cribs([a if i % 2 == 0 else norm_crib(a, False)
                            for i, a in enumerate(sys.argv[2:])])
        status, hyps, n = test_crib(blocks, cribs)
        print(status, n, json.dumps(hyps[:5], default=str)[:2000])
    elif mode == 'file':
        from multiprocessing import Pool
        jobs = []
        for line in open(sys.argv[2]):
            if not line.strip() or line.startswith('#'):
                continue
            fam, text = line.rstrip('\n').split('\t', 1)
            for spaces in (False, True):
                jobs.append((fam, text, spaces))
        out = []
        with Pool(os.cpu_count()) as pool:
            for rec in pool.imap_unordered(real_one, jobs):
                out.append(rec)
                print(rec['fam'], rec['status'], repr(rec['crib']), rec['s'], 's', flush=True)
        json.dump(out, open(os.path.join(ROOT, 'data', 'cribfit_real.json'), 'w'),
                  indent=1, default=str)


def real_one(job):
    """One real crib (family, text, spaces) against the real isomorph blocks."""
    from tournament import N9
    fam, text, spaces = job
    blocks, _ = blocks_from(load_msgs(), N9)
    t = norm_crib(text, spaces)
    t0 = time.time()
    status, hyps, n = test_crib(blocks, real_cribs([fam, t]))
    return {'fam': fam, 'text': text, 'spaces': spaces, 'crib': t[:FAMILIES[fam][1] - 1],
            'status': status, 'hyps_tried': n, 'hyps': hyps[:3], 's': round(time.time() - t0, 1)}


if __name__ == '__main__':
    main()
