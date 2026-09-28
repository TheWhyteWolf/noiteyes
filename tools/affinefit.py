#!/usr/bin/env python3
"""Step F: affine-pi solver -- is every sigma an affine map in some hidden order?

Every prior solve was translation-only (sigma_k: pi(b) = pi(a) + delta_k) and
the real core system is inconsistent under it (walk_dim == 1). The tournament's
best candidate, gak41 (Group Autokey over C83:C41, the one group the
community's GAK classification leaves standing), predicts instead

    pi(b) = beta_k * pi(a) + alpha_k      on every aligned pair (a, b) of block k,

with beta_k in the order-41 subgroup of Z83* -- which is exactly the set of
nonzero QUADRATIC RESIDUES mod 83 (index-2 subgroup of a cyclic group). The
alpha_k are free. For a fixed beta assignment the system is linear over
GF(83) in (pi, alpha); the solver branches over beta per block (41 choices)
and prunes on two exact, deterministic conditions:

  * consistency: the null space restricted to the touched letters must hold a
    NON-constant vector (the constant pi=1, alpha_k=1-beta_k always solves);
  * injectivity: two letters are forced to the same pi value exactly when
    their rows in the null-space basis are identical; pi is a bijection, so
    any forced merge kills the branch.

Blocks are the same witness-span core blocks sigma_web.py uses (jump-safe:
a mechanism change strictly inside a span breaks the mirrored repeat unless
two distinct affine maps agree on the witness letter, ~1/83 per witness).
Each block gets its OWN (beta, alpha): no assumption that blocks of one pair,
or occurrences of one web, share state.

Control gates (all must PASS before the real verdict is read):
  gak41  -- synthetic C83:C41 group-autokey corpus: the search must find the
            true beta assignment and recover pi up to affine gauge (u*pi+v)
            on every pinned letter.
  transl -- additive walk (all beta = 1): must be accepted with beta = 1
            on every block (translation is the beta=1 special case).
  fresh  -- unrelated random alphabet per stretch (no hidden order at all):
            must yield NO surviving assignment (the solver cannot invent
            affine structure that isn't there).

Gauge note: pi -> u*pi + v leaves every beta unchanged, so beta assignments
are gauge-invariant and directly comparable to ground truth.

Plaintext note: the game is Finnish; any decrypt-level check must score
Finnish as well as English (IoC ~0.07+ vs English ~0.066).
"""
import json
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
from pi_solver import P, find_isomorph_pairs, load_msgs
from sigma_web import MINLEN, PA, ROOT, WEIGHTS, core_blocks, witness_trim

SEED = 1
QR = sorted({pow(x, 2, P) for x in range(1, P)})
assert len(QR) == 41 and QR == sorted({pow(4, k, P) for k in range(41)})
NODE_CAP = 2_000_000


# ------------------------------------------------------------------ search

class Echelon:
    """Incrementally maintained reduced row-echelon form over GF(P).
    add(row) is O(rank * cols); null_basis() reads the null space off the
    RREF directly. Immutable-by-copy so the search can backtrack cheaply."""

    def __init__(self, ncols):
        self.ncols = ncols
        self.R = np.zeros((0, ncols), dtype=np.int64)
        self.piv = []                           # pivot column per row

    def copy(self):
        e = Echelon.__new__(Echelon)
        e.ncols, e.R, e.piv = self.ncols, self.R.copy(), list(self.piv)
        return e

    def add(self, row):
        r = row % P
        for i, c in enumerate(self.piv):
            if r[c]:
                r = (r - r[c] * self.R[i]) % P
        nz = np.nonzero(r)[0]
        if nz.size == 0:
            return
        c = int(nz[0])
        r = (r * pow(int(r[c]), P - 2, P)) % P
        if self.R.shape[0]:
            col = self.R[:, c].copy()
            hit = np.nonzero(col)[0]
            if hit.size:
                self.R[hit] = (self.R[hit] - np.outer(col[hit], r)) % P
        self.R = np.vstack([self.R, r])
        self.piv.append(c)

    def null_basis(self):
        pset = set(self.piv)
        free = [c for c in range(self.ncols) if c not in pset]
        B = np.zeros((self.ncols, len(free)), dtype=np.int64)
        for j, fc in enumerate(free):
            B[fc, j] = 1
        if self.piv:
            B[self.piv] = (-self.R[:, free]) % P
        return B


class Problem:
    """blocks: list of lists of (a, b) aligned letter pairs."""

    def __init__(self, blocks):
        self.blocks = [sorted(set(b)) for b in blocks]
        letters = sorted({x for b in self.blocks for pr in b for x in pr})
        self.letters = letters
        self.lidx = {v: i for i, v in enumerate(letters)}
        self.nlet = len(letters)
        self.ncols = self.nlet + len(self.blocks)
        order, touched = [], set()
        rest = set(range(len(self.blocks)))
        while rest:                             # biggest first, then overlap
            k = max(rest, key=lambda k: (len(self.bl_letters(k) & touched),
                                         len(self.blocks[k]), -k))
            order.append(k)
            touched |= self.bl_letters(k)
            rest.discard(k)
        self.order = order

    def bl_letters(self, k):
        return {x for pr in self.blocks[k] for x in pr}

    def rows_for(self, k, beta):
        out = []
        for a, b in self.blocks[k]:
            r = np.zeros(self.ncols, dtype=np.int64)
            r[self.lidx[b]] += 1
            r[self.lidx[a]] -= beta
            r[self.nlet + k] -= 1
            out.append(r % P)
        return out

    def ok(self, ech, touched):
        """Non-constant solution on the touched letters AND no two touched
        letters forced equal (identical null-basis rows <=> equal in every
        solution). The constant solution alone collapses all rows to one."""
        sub = ech.null_basis()[sorted(self.lidx[x] for x in touched)]
        return len({row.tobytes() for row in sub}) == len(touched)

    def solve(self, max_skip=0, node_cap=NODE_CAP, want=50, first=None):
        """All full assignments (block -> beta, or None = skipped) using at
        most max_skip skips. `first` pins the first block's beta (used to
        split the search across processes). Returns (sols, nodes, capped)."""
        sols, nodes = [], [0]
        order = self.order

        def rec(depth, ech, touched, assign, skips):
            if nodes[0] >= node_cap or len(sols) >= want:
                return
            if depth == len(order):
                sols.append(dict(assign))
                return
            k = order[depth]
            t2 = touched | self.bl_letters(k)
            for beta in (QR if depth or first is None else [first]):
                nodes[0] += 1
                e2 = ech.copy()
                for r in self.rows_for(k, beta):
                    e2.add(r)
                if self.ok(e2, t2):
                    assign[k] = beta
                    rec(depth + 1, e2, t2, assign, skips)
                    del assign[k]
            if skips < max_skip:
                assign[k] = None
                rec(depth + 1, ech, touched, assign, skips + 1)
                del assign[k]

        rec(0, Echelon(self.ncols), set(), {}, 0)
        return sols, nodes[0], nodes[0] >= node_cap

    def psolve(self, node_cap=NODE_CAP, workers=None):
        """solve(max_skip=0) split over the first block's 41 betas."""
        from multiprocessing import Pool
        with Pool(workers or os.cpu_count()) as pool:
            parts = pool.starmap(_solve_part, [(self, b, node_cap) for b in QR])
        sols = [s for p in parts for s in p[0]]
        return sols, sum(p[1] for p in parts), any(p[2] for p in parts)

    def min_skip_solve(self, max_budget=4, node_cap=NODE_CAP):
        """Iterative deepening on the skip budget: the smallest number of
        blocks that must be dropped for an affine-QR assignment to exist."""
        log = []
        for s in range(max_budget + 1):
            sols, nodes, capped = self.solve(max_skip=s, node_cap=node_cap)
            log.append({'skip': s, 'n_sols': len(sols), 'nodes': nodes,
                        'capped': capped})
            if sols or capped:
                return s, sols, log
        return None, [], log

    def pi_for(self, assign, rng):
        """Generic solution for a full assignment; pinned letters are those
        whose value is fixed up to affine gauge (identical across random
        null-space samples after normalizing two reference letters)."""
        ech = Echelon(self.ncols)
        used = set()
        for k, beta in assign.items():
            if beta is not None:
                for r in self.rows_for(k, beta):
                    ech.add(r)
                used |= self.bl_letters(k)
        L = ech.null_basis()[:self.nlet]
        idx = sorted(self.lidx[x] for x in used)
        samples = []
        for _ in range(8):
            c = np.array([rng.randrange(P) for _ in range(L.shape[1])])
            samples.append((L @ c) % P)
        a = idx[0]
        b = next(i for i in idx[1:] if all((s[i] - s[a]) % P for s in samples))
        norm = [((s - s[a]) * pow(int((s[b] - s[a]) % P), P - 2, P)) % P
                for s in samples]
        pinned = [i for i in idx if len({int(n[i]) for n in norm}) == 1]
        return {self.letters[i]: int(norm[0][i]) for i in pinned}, L.shape[1]


def _solve_part(prob, beta, node_cap):
    return prob.solve(max_skip=0, node_cap=node_cap, first=beta)


# -------------------------------------------------------- corpus -> blocks

WITNESS_DEPTH = 2


def excursions(S, names, min_run=4, max_gap=6):
    """Short EXCURSIONS: two messages whose ciphertexts are identical (same
    state) on runs of >= min_run letters at one relative offset, but differ on
    a gap of <= max_gap letters between two such runs. The state left and came
    back -- e.g. east-1/west-1 at 25-28 and 33-36. Witness spans cannot see a
    temporary change that returns (repeats on either side still mirror), so
    these positions are excluded from every block. Chance false positives:
    1 in 30 synthetic corpora (gak41/additive/quagA, 10 seeds each)."""
    out, found = set(), []
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            s1, s2 = S[names[a]], S[names[b]]
            for d in range(-len(s2) + 1, len(s1)):
                runs, cur, i = [], None, max(0, d)
                while i < len(s1) and i - d < len(s2):
                    if s1[i] == s2[i - d]:
                        cur = [i, i + 1] if cur is None else [cur[0], i + 1]
                    else:
                        if cur and cur[1] - cur[0] >= min_run:
                            runs.append(cur)
                        cur = None
                    i += 1
                if cur and cur[1] - cur[0] >= min_run:
                    runs.append(cur)
                for r1, r2 in zip(runs, runs[1:]):
                    if 0 < r2[0] - r1[1] <= max_gap:
                        found.append((names[a], names[b], d, r1[1], r2[0]))
                        for p in range(r1[1], r2[0]):
                            out.add((names[a], p))
                            out.add((names[b], p - d))
    return out, found


def blocks_from(S, names, need=WITNESS_DEPTH, cut_excursions=True):
    """Aligned-pair blocks: maximal runs of window positions covered by at
    least `need` distinct witness spans (a witness = a letter repeating inside
    the window, mirrored on the other side).

    Why depth 2, not sigma_web's depth 1: the isomorph finder extends a window
    until a repeat MISmatches, so windows reach back into the region where the
    two plaintexts have not yet converged whenever one repeat there happens to
    match by chance -- and that lone chance witness is then the only thing
    "protecting" the pre-convergence positions. Measured on synthetic gak41
    (20 seeds, known pi): depth 1 leaves 9/377 blocks non-affine under the
    TRUE pi (all prefix/suffix contamination of exactly this kind); depth 2
    leaves 0/420. Two independent chance matches over the same positions are
    ~1/83^2. Excursion positions (see excursions()) are cut out first."""
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    bad = excursions(S, names)[0] if cut_excursions else set()
    out, meta = [], []
    for n1, l1, r1, n2, l2, r2 in sorted(raw):
        w = S[n1][l1:r1]
        pos = {}
        for j, a in enumerate(w):
            pos.setdefault(a, []).append(j)
        depth = [0] * len(w)
        for p in pos.values():
            if len(p) >= 2:
                for j in range(p[0], p[-1] + 1):
                    depth[j] += 1
        for j in range(len(w)):
            if (n1, l1 + j) in bad or (n2, l2 + j) in bad:
                depth[j] = 0
        j = 0
        while j < len(w):
            if depth[j] < need:
                j += 1
                continue
            k = j
            while k < len(w) and depth[k] >= need:
                k += 1
            if k - j >= 2:
                out.append([(S[n1][l1 + i], S[n2][l2 + i]) for i in range(j, k)])
                meta.append((n1, l1 + j, n2, l2 + j, k - j))
            j = k
    return out, meta


def make_fresh(seed):
    """Control: a brand-new random alphabet per plaintext-notch stretch (no
    hidden order, no reuse). Same corpus harness as gak41, so its constraint
    density is comparable -- the power check the sparse alphacount corpus
    could not provide."""
    rng = random.Random(seed + 11)
    notch = set(rng.sample(PA, 5))
    def enc(s):
        alph = dict(zip(PA, rng.sample(range(P), len(PA))))
        seq = []
        for ch in s:
            seq.append(alph[ch])
            if ch in notch:
                alph = dict(zip(PA, rng.sample(range(P), len(PA))))
        return seq
    return enc


# ---------------------------------------------------------- ground truths

def gak41_truth(seed):
    """Replays tournament.make_gak41's RNG to recover its deck D (pi = D^-1)."""
    import tournament as T
    rng = random.Random(seed + 9)
    sub = sorted({pow(4, k, P) for k in range(41)})
    twist = set(rng.sample(T.PA, 5))
    {c: rng.randrange(1, P) for c in T.PA}
    [rng.choice([x for x in sub if x != 1]) if c in twist else 1 for c in T.PA]
    D = list(range(P))
    rng.shuffle(D)
    return {c: i for i, c in enumerate(D)}


def block_betas(prob, pi):
    """Ground-truth beta per block under a known pi (None if the block is not
    affine-with-QR-beta under it)."""
    out = {}
    for k, prs in enumerate(prob.blocks):
        rs = set()
        a1, b1 = prs[0]
        for a, b in prs[1:]:
            da = (pi[a] - pi[a1]) % P
            if da:
                rs.add(((pi[b] - pi[b1]) * pow(da, P - 2, P)) % P)
        out[k] = rs.pop() if len(rs) == 1 and rs <= set(QR) else None
    return out


def gauge_match(pinned, pi_true):
    """Fraction of pinned letters on which recovered pi = u*pi_true + v for
    the (u, v) fitted on the first two distinct pinned letters."""
    ks = [k for k in pinned if k in pi_true]
    if len(ks) < 3:
        return 0.0
    x0, x1 = ks[0], next((k for k in ks[1:] if (pi_true[k] - pi_true[ks[0]]) % P), None)
    if x1 is None:
        return 0.0
    u = ((pinned[x1] - pinned[x0]) * pow((pi_true[x1] - pi_true[x0]) % P, P - 2, P)) % P
    v = (pinned[x0] - u * pi_true[x0]) % P
    return sum(1 for k in ks if (u * pi_true[k] + v) % P == pinned[k]) / len(ks)


# ------------------------------------------------------------------ driver

def gate(name, ok):
    print(f'  {"PASS" if ok else "FAIL":4s}  {name}', flush=True)
    return ok


def run(label, S, names, node_cap=NODE_CAP, cut=True):
    blocks, meta = blocks_from(S, names, cut_excursions=cut)
    prob = Problem(blocks)
    neq = len({pr for b in prob.blocks for pr in b})
    t0 = time.time()
    sols, nodes, capped = prob.psolve(node_cap=node_cap)
    dt = time.time() - t0
    print(f'  {label}: {len(blocks)} blocks, {prob.nlet} letters, {neq} aligned '
          f'pairs, {prob.nlet + 2 * len(blocks)} unknowns -> {len(sols)} '
          f'assignment(s){" (CAPPED)" if capped else ""}, {nodes} nodes, {dt:.1f}s',
          flush=True)
    return prob, meta, sols, nodes, capped


def determination(prob, sols, rng):
    """How far a set of surviving assignments pins pi down: null-space
    dimension and pinned-letter count of the first solution (2 pinned = only
    the gauge reference letters, i.e. nothing determined)."""
    if not sols:
        return {'null_dim': None, 'n_pinned': 0}
    pinned, nd = prob.pi_for(sols[0], rng)
    return {'null_dim': nd, 'n_pinned': len(pinned), 'pi': pinned}


def main():
    import tournament as T
    rng = random.Random(SEED)
    report, gates = {'qr': QR, 'witness_depth': WITNESS_DEPTH}, []

    print('== controls (all on the tournament corpus harness) ==', flush=True)
    # gak41: must find the truth, uniquely, and recover pi up to affine gauge
    S, names = T.corpus(T.make_gak41(SEED), random.Random(SEED))
    prob, meta, sols, nodes, capped = run('gak41', S, names)
    pi_t = gak41_truth(SEED)
    tb = block_betas(prob, pi_t)
    hit = [s for s in sols if all(s[k] == tb[k] for k in s)]
    gm = gauge_match(prob.pi_for(hit[0], rng)[0], pi_t) if hit else 0.0
    report['gak41'] = {'n_sols': len(sols), 'nodes': nodes, 'capped': capped,
                       'truth_clean': all(v is not None for v in tb.values()),
                       'truth_found': bool(hit), 'gauge_match': gm}
    gates.append(gate('gak41: every block affine-QR under the true pi (blocks clean)',
                      all(v is not None for v in tb.values())))
    gates.append(gate('gak41: true beta assignment found, search not capped',
                      bool(hit) and not capped))
    gates.append(gate(f'gak41: pi recovered up to affine gauge (match {gm:.2f} == 1.0)',
                      gm == 1.0))

    # transl: the additive walk (translates of one deck) is the beta=1 case.
    # A direct check, not a search: many other assignments also fit a walk,
    # so enumerating them is pointless; what matters is that beta=1 passes.
    S, names = T.corpus(T.make_additive(SEED), random.Random(SEED + 1))
    blocks, _meta = blocks_from(S, names)
    prob = Problem(blocks)
    ech = Echelon(prob.ncols)
    for k in range(len(blocks)):
        for r in prob.rows_for(k, 1):
            ech.add(r)
    ok1 = prob.ok(ech, set(prob.letters))
    print(f'  transl: {len(blocks)} blocks, {prob.nlet} letters; all-beta=1 '
          f'assignment {"accepted" if ok1 else "REJECTED"}', flush=True)
    report['transl'] = {'n_blocks': len(blocks), 'beta1_accepted': ok1}
    gates.append(gate('transl: all-beta=1 (pure translation) assignment accepted', ok1))

    # fresh: unrelated alphabets at comparable density -> nothing may survive
    S, names = T.corpus(make_fresh(SEED), random.Random(SEED + 2))
    prob, meta, sols, nodes, capped = run('fresh', S, names, node_cap=100_000)
    report['fresh'] = {'n_sols': len(sols), 'nodes': nodes, 'capped': capped,
                       'n_blocks': len(prob.blocks),
                       'n_pairs': len({pr for b in prob.blocks for pr in b}),
                       'rejection_proven': not sols and not capped}
    # Exhausting this dense search takes hours; within budget the gate only
    # asserts no spurious fit is FOUND. A REFUTED real verdict additionally
    # requires the rejection to be proven (not capped) -- see below.
    gates.append(gate('fresh: no affine-QR assignment found for unrelated alphabets'
                      + ('' if not capped else ' (budget-capped: rejection NOT proven)'),
                      not sols))

    allok = all(gates)
    print('\n== control gates ==' + ('  all PASS' if allok else
          '\n  gate failure -> real verdict NOT trustworthy'), flush=True)

    print('\n== REAL corpus ==', flush=True)
    S = load_msgs()
    print(f'  excursions cut: {excursions(S, pi_solver.N9)[1]}')
    real = {'excursions': excursions(S, pi_solver.N9)[1]}
    for tag, cut in (('uncut', False), ('cut', True)):
        prob, meta, sols, nodes, capped = run(f'real ({tag})', S, pi_solver.N9, cut=cut)
        det = determination(prob, sols, rng)
        print(f'    null_dim {det["null_dim"]}, pinned letters {det["n_pinned"]}'
              f'/{prob.nlet}', flush=True)
        real[tag] = {'n_blocks': len(prob.blocks), 'n_letters': prob.nlet,
                     'n_pairs': len({pr for b in prob.blocks for pr in b}),
                     'blocks': [list(m) for m in meta], 'n_sols': len(sols),
                     'nodes': nodes, 'capped': capped,
                     'null_dim': det['null_dim'], 'n_pinned': det['n_pinned']}
    c = real['cut']
    if c['capped']:
        verdict = 'UNDETERMINED (search capped)'
    elif c['n_sols'] == 0 and not report['fresh']['rejection_proven']:
        verdict = ('NO FIT, but the solver rejection power is unproven at '
                   'this budget (fresh control capped) -- not a refutation')
    elif c['n_sols'] == 0:
        verdict = ('REFUTED: no affine map with square multiplier fits the real '
                   'sigma web -- C83:C41 group autokey (gak41) is ruled out')
    elif c['n_pinned'] >= 0.9 * c['n_letters']:
        verdict = ('DETERMINED: an affine-QR pi is pinned -- proceed to the '
                   'STOP-and-decrypt checks (score Finnish AND English)')
    else:
        verdict = (f'INCONCLUSIVE: the real web fits the affine model but is far '
                   f'too sparse to pin it ({c["n_sols"]}+ assignments, null dim '
                   f'{c["null_dim"]}, {c["n_pinned"]}/{c["n_letters"]} letters '
                   f'pinned). The uncut web conflicts ONLY through the '
                   f'east-1/west-1 excursions, which break the witness-span '
                   f'assumption. gak41 is neither refuted nor recoverable from '
                   f'the isomorphs alone.')
    real['verdict'] = verdict
    report['real'] = real
    print(f'\n  verdict: {verdict}')

    json.dump(report, open(f'{ROOT}/data/affinefit_report.json', 'w'),
              indent=1, default=str)
    print('\nwrote data/affinefit_report.json')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
