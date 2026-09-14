#!/usr/bin/env python3
"""Sigma-algebra diagnostics on the isomorph web + witness-core walk refit.

PLAN.md next-phase item 1, plus a corrected redo of the walk refutation
(item 2's foundation). Model-free: extracts the letter-mappings sigma between
isomorph occurrences and asks what algebraic family they belong to.
Reproducible from data/messages.json alone; synthetic controls with known
ground truth validate every diagnostic before the real-data verdict is read.

METHODOLOGICAL REPAIR ESTABLISHED BY THE CONTROLS HERE: the maximal isomorph
windows overrun the true shared phrase through repeat-free stretches (window
extension only stops at an UNMATCHED repeat), so full-window equation systems
go inconsistent even when the strict walk model is TRUE (the strict-walk
control reproduces pi_solver.py's "dim 1 refuted" on synthetic data where the
model holds by construction). pi_solver.py's refutation is therefore unsafe.

The safe system: a letter repeating inside a window at positions p<q (a
"witness") has its repeat mirrored on the other side, which forces equal
state offsets at p and q. With jumps that cannot cancel exactly, no jump can
sit strictly inside a witness span without breaking the window, so witness
spans — merged only when STRICTLY overlapping; a merge at a touch point
would leave the boundary gap unprotected — are jump-free segments for ANY
piecewise-translation mechanism, at any jump density. The core-block system
must therefore be consistent, with an injective pi on its PINNED letters
(letters whose value is forced up to affine gauge, i.e. stable across random
nullspace samples) and a plaintext-level difference-stream IoC, for every
bijective-deck walk, strict or jumped. A dim-1 core system on the real
corpus refutes the whole piecewise-translation family (up to exactly
cancelling jump pairs inside a span), not just the strict walk.

Pre-registered signatures (validated on the controls each run):
                    strict      sparse-jump  row-rekey     affine   shuffle
  pairs L>=15       many        many         aligned only  ~0       ~0
  gauge             consistent  consistent   consistent    -        -
  fixed points      0           ~0           ~0            -        -
  commutation       ~100% (a)   high         high          -        -
  closure zones     0 (a)       at jumps     at row ends   -        -
  full-window GF    dim 1 (a)   dim 1        dim 1         -        -
  core GF pinned    inj, IoC>   inj, IoC>    inj, IoC>     -        -
                    .05, truth   .05, truth   .05, truth
                    match ~1     match ~1     match ~1
  (a) up to overrun contamination, quantified by the strict control.
"""
import json
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
from pi_solver import P, find_isomorph_pairs, gauss_nullspace, load_msgs

ROOT = pi_solver.D
SEED = 1
MINLEN = 15


# ------------------------------------------------------ witness trimming

def witness_trim(S, raw_pairs, min_span=6):
    """Per pair: witness spans (merged only on STRICT overlap) and the trimmed
    window. Span [first,last+1) of a repeated letter protects the gaps
    first+1..last; strict overlap keeps protection contiguous."""
    joint = []
    for (n1, l1, r1, n2, l2, r2) in raw_pairs:
        w = S[n1][l1:r1]
        pos = defaultdict(list)
        for j, a in enumerate(w):
            pos[a].append(j)
        ivs = sorted((p[0], p[-1] + 1) for p in pos.values() if len(p) >= 2)
        if not ivs:
            continue
        spans = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s < spans[-1][1]:
                spans[-1][1] = max(spans[-1][1], e)
            else:
                spans.append([s, e])
        t0, t1 = spans[0][0], spans[-1][1]
        if t1 - t0 < min_span:
            continue
        joint.append({'raw': (n1, l1, r1, n2, l2, r2),
                      'trim': (n1, l1 + t0, l1 + t1, n2, l2 + t0, l2 + t1),
                      'spans': [tuple(sp) for sp in spans],
                      'covered': sum(e - s for s, e in spans) / (r1 - l1)})
    return joint


def core_blocks(joint):
    out, tags = [], []
    for k, ent in enumerate(joint):
        n1, l1, r1, n2, l2, r2 = ent['raw']
        for s, e in ent['spans']:
            if e - s >= 2:
                out.append((n1, l1 + s, n2, l2 + s, e - s))
                tags.append(k)
    return out, tags


# ---------------------------------------------------------------- webs

def build_webs(pairs):
    """Cluster windows into occurrence nodes; assign frame offsets (potentials)
    per connected component; detect alignment slippage (gauge violations)."""
    occs = []                                   # [msg, lo, hi]

    def occ_of(msg, l, r):
        for k, (m, lo, hi) in enumerate(occs):
            if m == msg and min(hi, r) - max(lo, l) >= 0.5 * min(r - l, hi - lo):
                occs[k] = [m, min(lo, l), max(hi, r)]
                return k
        occs.append([msg, l, r])
        return len(occs) - 1

    edges = [(occ_of(n1, l1, r1), occ_of(n2, l2, r2), l2 - l1, k)
             for k, (n1, l1, r1, n2, l2, r2) in enumerate(pairs)]
    while True:                                 # stabilize clustering
        merge = next(((i, j) for i in range(len(occs)) for j in range(i + 1, len(occs))
                      if occs[i][0] == occs[j][0]
                      and min(occs[i][2], occs[j][2]) - max(occs[i][1], occs[j][1])
                      >= 0.5 * min(occs[i][2] - occs[i][1], occs[j][2] - occs[j][1])), None)
        if merge is None:
            break
        i, j = merge
        occs[i] = [occs[i][0], min(occs[i][1], occs[j][1]), max(occs[i][2], occs[j][2])]
        occs.pop(j)
        edges = [(a - (a > j), b - (b > j), o, k) for a, b, o, k in
                 ((i if a == j else a, i if b == j else b, o, k) for a, b, o, k in edges)]

    pot, web_id, violations = {}, {}, []
    nweb = 0
    adj = defaultdict(list)
    for a, b, off, k in edges:
        adj[a].append((b, off, k))
        adj[b].append((a, -off, k))
    for s in range(len(occs)):
        if s in pot:
            continue
        pot[s], web_id[s] = 0, nweb
        queue = [s]
        while queue:
            u = queue.pop(0)
            for vtx, off, k in adj[u]:
                if vtx not in pot:
                    pot[vtx], web_id[vtx] = pot[u] + off, nweb
                    queue.append(vtx)
                elif pot[vtx] != pot[u] + off:
                    violations.append((u, vtx, k, pot[vtx], pot[u] + off))
        nweb += 1
    webs = defaultdict(list)
    for o, w in web_id.items():
        webs[w].append(o)
    webs = [sorted(v) for v in webs.values() if len(v) >= 2]
    return occs, pot, webs, violations, edges


# ------------------------------------------------------- sigma statistics

def extract_sigmas(S, pairs):
    return [{S[n1][l1 + j]: S[n2][l2 + j] for j in range(r1 - l1)}
            for n1, l1, r1, n2, l2, r2 in pairs]


def cycle_stats(m):
    fixed = sum(1 for a, b in m.items() if a == b)
    cycles, seen = 0, set()
    for a in m:
        if a in seen:
            continue
        path, cur = {a}, a
        while True:
            cur = m.get(cur)
            if cur is None or cur in seen:
                break
            if cur == a:
                if len(path) > 1:
                    cycles += 1
                break
            if cur in path:
                break
            path.add(cur)
        seen |= path
    return fixed, cycles


def commutation_stats(sigmas):
    agree = tot = 0
    per = []
    for i in range(len(sigmas)):
        for j in range(i + 1, len(sigmas)):
            s1, s2 = sigmas[i], sigmas[j]
            a = t = 0
            for x in set(s1) & set(s2):
                p, q = s2.get(s1[x]), s1.get(s2[x])
                if p is not None and q is not None:
                    t += 1
                    a += (p == q)
            if t:
                agree, tot = agree + a, tot + t
                per.append([i, j, a, t])
    return agree, tot, per


# ------------------------------------------------------------- closure

def closure_analysis(S, occs, ext, pot, webs, trims_by_occ):
    """Mutual-isomorphy test over full (raw-extent) web coverage. A violating
    position pair (p,q) — repeat on one side, none on the other — brackets a
    net state-offset jump, or the end of the true shared phrase, inside
    (p, q]. Zones intersecting witness-trimmed coverage on both sides are
    tagged evidence-backed (phrase-end noise cannot live there)."""
    stats = {'cells': 0, 'clean': 0, 'viol': 0}
    zones = []
    for web in webs:
        for ii in range(len(web)):
            for jj in range(ii + 1, len(web)):
                i, j = web[ii], web[jj]
                m1, m2 = occs[i][0], occs[j][0]
                (s1, e1), (s2, e2) = ext[i], ext[j]
                p1, p2 = pot[i], pot[j]
                lo, hi = max(s1 - p1, s2 - p2), min(e1 - p1, e2 - p2)
                if hi - lo < 4:
                    continue
                A = S[m1][lo + p1:hi + p1]
                B = S[m2][lo + p2:hi + p2]
                L = hi - lo
                viols, cleans = [], []
                for p in range(L):
                    for q in range(p + 1, L):
                        ea, eb = A[p] == A[q], B[p] == B[q]
                        if not (ea or eb):
                            continue
                        stats['cells'] += 1
                        (cleans if ea == eb else viols).append((p, q))
                        stats['clean' if ea == eb else 'viol'] += 1
                ivs = sorted((p + 1, q) for p, q in viols)
                group = []
                for iv in ivs + [(None, None)]:
                    if group and (iv[0] is None or iv[0] > min(g[1] for g in group)):
                        z0 = max(g[0] for g in group)
                        z1 = min(g[1] for g in group)
                        if z1 < z0:
                            z0, z1 = min(g[0] for g in group), max(g[1] for g in group)
                        fz0, fz1 = lo + z0, lo + z1
                        evid = all(any(t0 - pot[o] <= fz1 and fz0 <= t1 - pot[o]
                                       for t0, t1 in trims_by_occ.get(o, []))
                                   for o in (i, j))
                        zones.append({'occ_a': i, 'occ_b': j, 'frame': (fz0, fz1),
                                      'n_viol': len(group), 'evidence': evid})
                        group = []
                    if iv[0] is not None:
                        group.append(iv)
    return stats, zones


# ------------------------------------------------- GF(p) block system

def gauss_p(rows, ncols, p):
    """gauss_nullspace over an arbitrary prime modulus."""
    rows = [dict(r) for r in rows if r]
    reduced, pivots = [], {}
    for r in rows:
        r = dict(r)
        for col, ri in list(pivots.items()):
            if col in r:
                factor = r[col]
                for c2, v2 in reduced[ri].items():
                    r[c2] = (r.get(c2, 0) - factor * v2) % p
                r = {c: v for c, v in r.items() if v}
        if not r:
            continue
        piv = min(r)
        inv = pow(r[piv], p - 2, p)
        r = {c: (v * inv) % p for c, v in r.items()}
        for i, rr in enumerate(reduced):
            if piv in rr:
                factor = rr[piv]
                for c2, v2 in r.items():
                    rr[c2] = (rr.get(c2, 0) - factor * v2) % p
                reduced[i] = {c: v for c, v in rr.items() if v}
        reduced.append(r)
        reduced = [rr for rr in reduced if rr]
        pivots = {min(rr): i for i, rr in enumerate(reduced)}
    rank = len(reduced)
    free_cols = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free_cols:
        vec = [0] * ncols
        vec[fc] = 1
        for pc, ri in pivots.items():
            vec[pc] = (-reduced[ri].get(fc, 0)) % p
        basis.append(vec)
    return rank, free_cols, basis


def solve_blocks(S, blocks, rng, tries=60, p=P):
    letters = sorted({S[m][o + j] for n1, l1, n2, l2, L in blocks
                      for m, o in ((n1, l1), (n2, l2)) for j in range(L)})
    lidx = {v: i for i, v in enumerate(letters)}
    nlet, ncols = len(letters), len(letters) + len(blocks)
    rows = []
    for k, (n1, l1, n2, l2, L) in enumerate(blocks):
        for a, b in {(S[n1][l1 + j], S[n2][l2 + j]) for j in range(L)}:
            row = defaultdict(int)
            row[lidx[b]] += 1
            row[lidx[a]] -= 1
            row[nlet + k] -= 1
            r = {c: v % p for c, v in row.items() if v % p}
            if r:
                rows.append(r)
    rank, _, basis = gauss_p(rows, ncols, p)
    dim = ncols - rank
    best = None
    samples = []
    for _ in range(tries if dim >= 2 else 0):
        coef = [rng.randrange(p) for _ in basis]
        vec = [0] * ncols
        for c, bv in zip(coef, basis):
            if c:
                for idx, bvi in enumerate(bv):
                    if bvi:
                        vec[idx] = (vec[idx] + c * bvi) % p
        lp = vec[:nlet]
        if len(set(lp)) > 1:
            samples.append(lp)
            if best is None or len(set(lp)) > best[0]:
                best = (len(set(lp)), vec)
    return {'dim': dim, 'rank': rank, 'nlet': nlet, 'nblocks': len(blocks),
            'letters': letters, 'neqs': len(rows), 'samples': samples[:16],
            'distinct': best[0] if best else 1, 'p': p,
            'pi': ({letters[i]: best[1][i] for i in range(nlet)} if best else None)}


def forced_merges(res):
    """Letters identical in EVERY nullspace sample are forced equal by the
    system (chance identity across >=6 samples is ~p^-6). Count letters
    merged beyond one representative per class."""
    smp = res['samples']
    if not smp:
        return res['nlet'] - 1
    keys = list(zip(*smp))
    return len(keys) - len(set(keys))


def pinned_eval(S, names, res):
    """Letters whose value is forced by the system up to affine gauge: stable
    across random nullspace samples after gauge normalization. The decisive
    metrics live on this set only — floaters in underdetermined systems would
    otherwise dilute IoC/injectivity with arbitrary values."""
    samps, letters, p = res['samples'], res['letters'], res['p']
    if len(samps) < 4:
        return None
    n = len(letters)
    best = None
    for a in range(min(6, n)):
        for b in range(a + 1, min(a + 7, n)):
            if any((s[b] - s[a]) % p == 0 for s in samps):
                continue
            norm = [[((s[i] - s[a]) * pow(s[b] - s[a], p - 2, p)) % p for i in range(n)]
                    for s in samps]
            pinned = [i for i in range(n) if len({nm[i] for nm in norm}) == 1]
            if best is None or len(pinned) > len(best[0]):
                best = (pinned, norm[0])
    if best is None:
        return None
    pinned, norm0 = best
    pi = {letters[i]: norm0[i] for i in pinned}
    ioc, nt = diff_ioc(S, names, pi, p)
    return {'n_pinned': len(pinned), 'nlet': n,
            'distinct': len({norm0[i] for i in pinned}),
            'ioc': [round(ioc, 5), nt], 'pi': pi}


def diff_ioc(S, names, pi, p=P):
    hist = defaultdict(int)
    for n in names:
        s = S[n]
        for i in range(2, len(s)):
            if s[i] in pi and s[i - 1] in pi:
                hist[(pi[s[i]] - pi[s[i - 1]]) % p] += 1
    n = sum(hist.values())
    return (sum(c * (c - 1) for c in hist.values()) / (n * (n - 1)) if n > 1 else 0.0), n


def pairwise_matrix(S, blocks, rng, p=P):
    """Mutual consistency of core segments: a pair conflicts when its joint
    system is inconsistent (dim<2) or forces letter merges."""
    conflicts = Counter()
    ncons = ntot = 0
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            r = solve_blocks(S, [blocks[i], blocks[j]], rng, tries=10, p=p)
            ok = r['dim'] >= 2 and forced_merges(r) == 0
            ntot += 1
            ncons += ok
            if not ok:
                conflicts[i] += 1
                conflicts[j] += 1
    return ncons, ntot, conflicts.most_common(5)


def consensus_prune(S, names, blocks, rng, min_pinned=25):
    """RANSAC-style consensus: jumps in BOTH streams of a pair can cancel in
    the relative offset, silently corrupting a witness span (the sparse
    control exhibits this), but a corrupted segment betrays itself by
    conflicting with clean segments sharing its letters. Iteratively drop the
    most-conflicted segment until the joint system pins an injective solution
    or nothing usable remains. Pruning only removes constraints, so surviving
    pinned values are forced by clean evidence alone."""
    cur, dropped = list(blocks), []
    history = []
    while len(cur) >= 3:
        res = solve_blocks(S, cur, rng, tries=60)
        pin = pinned_eval(S, names, res) if res['dim'] >= 2 else None
        history.append((len(cur), res['dim'], pin['n_pinned'] if pin else 0))
        if pin and pin['distinct'] >= 0.95 * pin['n_pinned'] and pin['n_pinned'] >= min_pinned:
            return {'verdict': 'CONSISTENT', 'kept': cur, 'dropped': dropped,
                    'res': res, 'pin': pin, 'history': history}
        ncons, ntot, top = pairwise_matrix(S, cur, rng)
        if not top:
            return {'verdict': ('UNDERDETERMINED' if res['dim'] >= 2 else 'INCONSISTENT'),
                    'kept': cur, 'dropped': dropped, 'res': res, 'pin': pin,
                    'history': history, 'pairwise': [ncons, ntot]}
        dropped.append(cur.pop(top[0][0]))
    return {'verdict': 'EXHAUSTED', 'kept': cur, 'dropped': dropped,
            'res': None, 'pin': None, 'history': history}


def sol_str(r):
    return (f"dim {r['dim']:2d}  distinct {r['distinct']}/{r['nlet']}"
            + (' DEGENERATE' if r['pi'] and r['distinct'] < 0.7 * r['nlet'] else ''))


# --------------------------------------------------------- row analysis

def rowpair_boundaries(S):
    """Letter positions where a new physical row-pair (triangle row) starts,
    from raw_eye_rows in messages.json."""
    raw = json.load(open(f'{ROOT}/data/messages.json'))['raw_eye_rows']
    bounds = {}
    for n, rows in raw.items():
        counts = [(len(rows[k]) + len(rows[k + 1])) // 3 for k in range(0, len(rows) - 1, 2)]
        if sum(counts) != len(S[n]):
            print(f'  WARNING: row-pair letters {sum(counts)} != {len(S[n])} for {n}')
        cum, acc = [], 0
        for c in counts[:-1]:
            acc += c
            cum.append(acc)
        bounds[n] = cum
    return bounds


def boundary_distance_stat(events, bounds, S):
    """events: (msg, pos). Observed mean distance to nearest row boundary vs
    the uniform-position null for the same messages."""
    obs, null_mu, null_var = [], 0.0, 0.0
    for n, p in events:
        B = bounds.get(n)
        if not B:
            continue
        obs.append(min(abs(p - b) for b in B))
        dists = [min(abs(t - b) for b in B) for t in range(len(S[n]))]
        m = sum(dists) / len(dists)
        null_mu += m
        null_var += sum((d - m) ** 2 for d in dists) / len(dists)
    if not obs:
        return None
    n = len(obs)
    mu, sd = null_mu / n, (null_var / n / n) ** 0.5
    om = sum(obs) / n
    return {'n': n, 'obs_mean': round(om, 2), 'null_mean': round(mu, 2),
            'z': round((om - mu) / sd, 2) if sd else None}


# ------------------------------------------------------------ simulators

PA = 'abcdefghijklmnopqrstuvwxyz_'
WEIGHTS = [8.2, 1.5, 2.8, 4.3, 12.7, 2.2, 2.0, 6.1, 7.0, 0.2, 0.8, 4.0, 2.4,
           6.7, 7.5, 1.9, 0.1, 6.0, 6.3, 9.1, 2.8, 1.0, 2.4, 0.2, 2.0, 0.1, 18.0]


def rich_phrase(rng, v, L, jump_offsets=()):
    """Random phrase whose walk image is repeat-dense enough to seed the
    finder (needs a 12-window with >=3 internally repeated letters; repeats
    depend only on v-sums, plus any fixed +8 jumps at the given offsets)."""
    while True:
        ph = ''.join(rng.choices(PA, weights=WEIGHTS, k=L))
        pre = [0]
        for k, ch in enumerate(ph):
            pre.append((pre[-1] + v[ch] + (8 if k in jump_offsets else 0)) % P)
        pos_by_state = defaultdict(list)
        for k in range(1, L + 1):
            pos_by_state[pre[k]].append(k - 1)
        reps = [ix for ix in pos_by_state.values() if len(ix) >= 2]
        if sum(len(ix) * (len(ix) - 1) // 2 for ix in reps) < 5:
            continue
        for w0 in range(0, L - 11):
            if sum(1 for ix in reps
                   if len([i for i in ix if w0 <= i < w0 + 12]) >= 2) >= 3:
                return ph


def synth(mode, seed):
    """9-message corpus mimicking the real structure. Modes:
    strict — pure additive walk; sparse — additive walk with rare random +1
    state jumps (rate .03); rowkey — additive walk with a fixed +8 state jump
    at every 13-letter row boundary (phrase copies at row-aligned positions,
    one misaligned demo copy); affine — x -> 5x + v (destroys exact
    isomorphy; expect no pairs)."""
    rng = random.Random(seed)
    deck = list(range(P))
    rng.shuffle(deck)
    v = {ch: rng.randrange(1, P) for ch in PA}
    text = lambda L: ''.join(rng.choices(PA, weights=WEIGHTS, k=L))
    names = [f'sim-{i}' for i in range(9)]
    plain = {n: list(text(rng.randrange(100, 138))) for n in names}
    ins = lambda n, ph, pos: plain[n].__setitem__(slice(pos, pos + len(ph)), list(ph))
    if mode == 'rowkey':
        # boundary hits letter j with (j+1) % 13 == 0; insert at 4 (mod 13)
        # puts them at phrase offsets {8, 21} for all phrase lengths
        phrase1 = rich_phrase(rng, v, 30, jump_offsets={8, 21})
        phrase2 = rich_phrase(rng, v, 34, jump_offsets={8, 21})
        phrase3 = rich_phrase(rng, v, 28, jump_offsets={8, 21})
        ins('sim-0', phrase1, 30); ins('sim-0', phrase1, 56)
        ins('sim-1', phrase1, 43); ins('sim-1', phrase1, 69)
        ins('sim-5', phrase1, 30); ins('sim-5', phrase1, 63)   # 63: misaligned
        ins('sim-4', phrase2, 30); ins('sim-8', phrase2, 43)
        ins('sim-3', phrase2[:22], 56)
        ins('sim-2', phrase3, 30); ins('sim-2', phrase3, 69); ins('sim-6', phrase3, 43)
    else:
        phrase1, phrase2 = rich_phrase(rng, v, 30), rich_phrase(rng, v, 34)
        phrase3 = rich_phrase(rng, v, 28)
        ins('sim-0', phrase1, 30); ins('sim-0', phrase1, 61)
        ins('sim-1', phrase1, 35); ins('sim-1', phrase1, 72)
        ins('sim-5', phrase1, 30)
        ins('sim-4', phrase2, 35); ins('sim-8', phrase2, 36); ins('sim-3', phrase2[:22], 46)
        ins('sim-2', phrase3, 40); ins('sim-2', phrase3, 75); ins('sim-6', phrase3, 33)
    S, jumps = {}, defaultdict(list)
    for i, n in enumerate(names):
        plain[n][0:2] = list('qk')
        x, out = 17, [(7 * i + 3) % P]
        for j, ch in enumerate(plain[n]):
            if mode == 'affine':
                x = (5 * x + v[ch]) % P
            else:
                x = (x + v[ch]) % P
                if mode == 'sparse' and rng.random() < 0.02:
                    x = (x + 1) % P
                    jumps[n].append(j + 1)
                if mode == 'dense':                 # signed, cancellation-prone
                    u = rng.random()
                    if u < 0.08:
                        x = (x + 8) % P
                        jumps[n].append(j + 1)
                    elif u < 0.16:
                        x = (x - 8) % P
                        jumps[n].append(j + 1)
                if mode == 'rowkey' and (j + 1) % 13 == 0:
                    x = (x + 8) % P
                    jumps[n].append(j + 1)
            out.append(deck[x])
        S[n] = out
    return S, {'jumps': dict(jumps), 'deck': deck}


def pi_matches_truth(pi, deck):
    """Fraction of recovered pi consistent with one affine map of true state."""
    pinv = {deck[s]: s for s in range(P)}
    items = [(pinv[c], val) for c, val in pi.items()]
    (s0, v0) = items[0]
    best = 0
    for s1, v1 in items[1:]:
        if (s1 - s0) % P == 0:
            continue
        a = ((v1 - v0) * pow(s1 - s0, P - 2, P)) % P
        b = (v0 - a * s0) % P
        best = max(best, sum(1 for s, val in items if (a * s + b) % P == val))
    return best / len(items)


# -------------------------------------------------------------- battery

def battery(S, names, label, truth=None, rng=None, bounds=None):
    rng = rng or random.Random(SEED)
    raw_pairs = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    M = {'label': label, 'pairs': len(raw_pairs)}
    print(f'\n== {label} ==')
    print(f'isomorph pairs L>={MINLEN}: {len(raw_pairs)}')
    if not raw_pairs:
        return M
    joint = witness_trim(S, raw_pairs)
    trimmed = [e['trim'] for e in joint]
    occs, pot, webs, gauge_viol, edges = build_webs(trimmed)

    ext = {k: (lo, hi) for k, (m, lo, hi) in enumerate(occs)}
    trims_by_occ = defaultdict(list)
    for a, b, off, k in edges:
        rw, tw = joint[k]['raw'], joint[k]['trim']
        ext[a] = (min(ext[a][0], rw[1]), max(ext[a][1], rw[2]))
        ext[b] = (min(ext[b][0], rw[4]), max(ext[b][1], rw[5]))
        trims_by_occ[a].append((tw[1], tw[2]))
        trims_by_occ[b].append((tw[4], tw[5]))

    sigmas = extract_sigmas(S, trimmed)
    fixed = sum(cycle_stats(m)[0] for m in sigmas)
    cycles = sum(cycle_stats(m)[1] for m in sigmas)
    agree, tot, commut_per = commutation_stats(sigmas)
    # commutation split by message family (position-3 letter, per PLAN.md)
    fam = {n: S[n][3] for n in names}
    def cat(i, j):
        msgs = {trimmed[i][0], trimmed[i][3], trimmed[j][0], trimmed[j][3]}
        fams = {fam[m] for m in msgs}
        return f'pure-{sorted(fams)[0]}' if len(fams) == 1 else 'mixed'
    csplit = defaultdict(lambda: [0, 0])
    for i, j, a, t in commut_per:
        c = cat(i, j)
        csplit[c][0] += a
        csplit[c][1] += t
    cstats, zones = closure_analysis(S, occs, ext, pot, webs, trims_by_occ)

    blocks_raw = [(n1, l1, n2, l2, r1 - l1) for n1, l1, r1, n2, l2, r2 in raw_pairs]
    cores, core_tags = core_blocks(joint)
    full = solve_blocks(S, blocks_raw, rng)
    core = solve_blocks(S, cores, rng)
    pinned = pinned_eval(S, names, core) if core['dim'] >= 2 else None

    cover = sum(e - s for ent in joint for s, e in ent['spans'])
    total = sum(p[2] - p[1] for p in raw_pairs)
    M.update(webs=[len(w) for w in webs], gauge_violations=len(gauge_viol),
             fixed_points=fixed, closed_cycles=cycles,
             commut=[agree, tot], commut_per=commut_per, closure=cstats,
             zones=len(zones),
             zones_evidence=sum(1 for z in zones if z['evidence']),
             zone_detail=[{**z, 'msg_a': occs[z['occ_a']][0],
                           'pos_a': [z['frame'][0] + pot[z['occ_a']],
                                     z['frame'][1] + pot[z['occ_a']]],
                           'msg_b': occs[z['occ_b']][0],
                           'pos_b': [z['frame'][0] + pot[z['occ_b']],
                                     z['frame'][1] + pot[z['occ_b']]]} for z in zones],
             occurrences=[{'msg': m, 'span': [lo, hi], 'pot': pot[k]}
                          for k, (m, lo, hi) in enumerate(occs)],
             witness_cover=round(cover / total, 3), n_cores=len(cores),
             full={k: full[k] for k in ('dim', 'distinct', 'nlet', 'neqs')},
             core={k: core[k] for k in ('dim', 'distinct', 'nlet', 'neqs')})
    if pinned:
        M['pinned'] = {k: pinned[k] for k in ('n_pinned', 'nlet', 'distinct', 'ioc')}

    print(f"webs: {M['webs']}  gauge violations: {M['gauge_violations']}  "
          f"witness cover {M['witness_cover']}  cores: {len(cores)}")
    M['commut_split'] = {k: list(v) for k, v in csplit.items()}
    print(f"sigma fixed points: {fixed}  closed cycles: {cycles}  "
          f"commutation: {agree}/{tot} "
          + ' '.join(f'[{k} {v[0]}/{v[1]}]' for k, v in sorted(csplit.items())))
    print(f"closure: {cstats['viol']} violations / {cstats['clean']} clean; "
          f"{M['zones']} zones, {M['zones_evidence']} evidence-backed")
    print(f"GF(83) raw windows:  {sol_str(full)}")
    print(f"GF(83) CORES:        {sol_str(core)}  ({core['neqs']} eqs, "
          f"{len(cores)} segments)")
    if pinned:
        print(f"pinned letters: {pinned['n_pinned']}/{pinned['nlet']}  "
              f"distinct {pinned['distinct']}  IoC {pinned['ioc']}")

    cons = consensus_prune(S, names, cores, rng)
    ncons0, ntot0, top0 = pairwise_matrix(S, cores, rng)
    M['pairwise_consistent'] = [ncons0, ntot0]
    M['consensus'] = {'verdict': cons['verdict'], 'n_kept': len(cons['kept']),
                      'n_dropped': len(cons['dropped']),
                      'dropped': [list(map(str, b)) for b in cons['dropped']],
                      'history': cons['history']}
    if cons.get('pin'):
        M['consensus'].update(n_pinned=cons['pin']['n_pinned'],
                              distinct=cons['pin']['distinct'],
                              ioc=cons['pin']['ioc'])
        M['consensus_pi'] = {str(k): val for k, val in cons['pin']['pi'].items()}
    print(f"pairwise segment consistency: {ncons0}/{ntot0}")
    print(f"consensus: {cons['verdict']} after dropping "
          f"{len(cons['dropped'])}/{len(cores)} segments"
          + (f"; pinned {cons['pin']['n_pinned']} distinct {cons['pin']['distinct']}"
             f" IoC {cons['pin']['ioc']}" if cons.get('pin') else ''))

    if cons['verdict'] != 'CONSISTENT':
        # deep probe 1: subsystems — per-message/per-web deck hypotheses
        occ_web = {o: wi for wi, w in enumerate(webs) for o in w}
        pair_web = {k: occ_web.get(a) for a, b, off, k in edges}
        subs = {'within-msg': [blk for blk, k in zip(cores, core_tags)
                               if joint[k]['raw'][0] == joint[k]['raw'][3]],
                'cross-msg': [blk for blk, k in zip(cores, core_tags)
                              if joint[k]['raw'][0] != joint[k]['raw'][3]]}
        for wi in sorted(set(pair_web.values())):
            subs[f'web{wi}'] = [blk for blk, k in zip(cores, core_tags)
                                if pair_web.get(k) == wi]
        M['subsystems'] = {}
        for nm, bl in subs.items():
            if len(bl) < 2:
                continue
            r = solve_blocks(S, bl, rng)
            pin = pinned_eval(S, names, r) if r['dim'] >= 2 else None
            M['subsystems'][nm] = {'segs': len(bl), 'dim': r['dim'],
                                   'merges': forced_merges(r) if r['dim'] >= 2 else None,
                                   'pinned': pin['n_pinned'] if pin else 0,
                                   'distinct': pin['distinct'] if pin else 0,
                                   'ioc': pin['ioc'] if pin else None}
            print(f"  subsystem {nm}: {len(bl)} segs  dim {r['dim']}"
                  + (f"  merges {forced_merges(r)}  pinned "
                     f"{pin['n_pinned'] if pin else 0}"
                     + (f" distinct {pin['distinct']} IoC {pin['ioc']}" if pin else '')
                     if r['dim'] >= 2 else ''))
        # deep probe 2: state group Z_N for N != 83 (bigger deck)
        sweep = {}
        for pr in (59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127):
            r = solve_blocks(S, cores, rng, tries=25, p=pr)
            sweep[pr] = r['dim']
            if r['dim'] >= 2 and pr != P:
                pin = pinned_eval(S, names, r)
                if pin:
                    print(f"  modulus {pr}: dim {r['dim']}  pinned {pin['n_pinned']} "
                          f"distinct {pin['distinct']}  IoC {pin['ioc']}")
        M['modulus_sweep'] = sweep
        print(f"  modulus sweep (dim per prime): {sweep}")

    if truth:
        pi = (cons.get('pin') or pinned or {}).get('pi') or core.get('pi')
        if pi and len(pi) >= 2:
            M['pi_truth_match'] = round(pi_matches_truth(pi, truth['deck']), 3)
        tj = {(n, t) for n, ts in truth['jumps'].items() for t in ts}
        wins = {(q[0], q[1], q[2]) for q in raw_pairs} | {(q[3], q[4], q[5]) for q in raw_pairs}
        det = {(n, t) for n, t in tj if any(pn == n and l < t < r for pn, l, r in wins)}
        zpos = defaultdict(list)
        for z in M['zone_detail']:
            zpos[z['msg_a']].append(z['pos_a'])
            zpos[z['msg_b']].append(z['pos_b'])
        hit = {(n, t) for n, t in det
               if any(a - 1 <= t <= b + 1 for a, b in zpos.get(n, []))}
        M['jump_recall'] = f'{len(hit)}/{len(det)} detectable (of {len(tj)} total)'
        print(f"ground truth: jump recall {M['jump_recall']}"
              f"  pi-vs-truth match {M.get('pi_truth_match')}")

    if bounds:
        ev_win = [(n, p) for (n1, l1, r1, n2, l2, r2) in raw_pairs
                  for n, p in ((n1, l1), (n1, r1), (n2, l2), (n2, r2))
                  if 0 < p < len(S[n]) - 1]
        ev_zone = [(z['msg_a'], pp) for z in M['zone_detail'] for pp in z['pos_a']] + \
                  [(z['msg_b'], pp) for z in M['zone_detail'] for pp in z['pos_b']]
        for tag, ev in (('window ends', ev_win), ('zones', ev_zone)):
            st = boundary_distance_stat(ev, bounds, S)
            M[f'rowdist_{tag.replace(" ", "_")}'] = st
            if st:
                print(f"row-boundary distance ({tag}): obs {st['obs_mean']} vs "
                      f"null {st['null_mean']} (n={st['n']}, z={st['z']})")
    return M


def main():
    S = load_msgs()
    names = pi_solver.N9
    rng = random.Random(SEED)
    report = {}

    gates = []
    for mode in ('strict', 'sparse', 'rowkey', 'dense', 'affine'):
        Ssim, truth = synth(mode, SEED)
        simb = {n: list(range(13, len(Ssim[n]) - 1, 13)) for n in Ssim} \
            if mode == 'rowkey' else None
        r = battery(Ssim, sorted(Ssim), f'control: {mode} walk',
                    truth=truth if mode != 'affine' else None, rng=rng, bounds=simb)
        report[mode] = r
        if mode == 'affine':
            gates.append(('affine: no exact isomorphs', r['pairs'] == 0))
        elif mode == 'dense':
            pass                    # informational: catalog richness comparison
        else:
            cn = r.get('consensus', {})
            ok = (cn.get('verdict') == 'CONSISTENT'
                  and cn.get('ioc', [0])[0] > 0.05
                  and r.get('pi_truth_match', 0) >= 0.95)
            gates.append((f'{mode}: consensus cores crack layer 1', ok))

    shuf = {}
    for n in names:
        s = list(S[n])
        rng.shuffle(s)
        shuf[n] = s
    npairs = len([p for p in find_isomorph_pairs(shuf, names=names)
                  if p[2] - p[1] >= MINLEN])
    print(f'\n== control: shuffled real corpus ==\nisomorph pairs L>={MINLEN}: {npairs}')
    report['shuffle'] = {'pairs': npairs}
    gates.append(('shuffle: no pairs', npairs == 0))

    print('\n== control gates ==')
    for name, ok in gates:
        print(f'  {"PASS" if ok else "FAIL"}  {name}')

    bounds = rowpair_boundaries(S)
    report['real'] = battery(S, names, 'REAL CORPUS', rng=rng, bounds=bounds)

    raw_pairs = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    joint = witness_trim(S, raw_pairs)
    sig = [{'pair': list(e['trim']), 'spans': [list(s) for s in e['spans']],
            'sigma': {str(a): b for a, b in m.items()}}
           for e, m in zip(joint, extract_sigmas(S, [e['trim'] for e in joint]))]
    json.dump(sig, open(f'{ROOT}/data/sigmas.json', 'w'), indent=1)
    for r in report.values():
        r.pop('commut_per', None)
    json.dump(report, open(f'{ROOT}/data/sigma_report.json', 'w'), indent=1, default=str)
    print(f'\nwrote data/sigmas.json ({len(sig)} sigma maps) and data/sigma_report.json')
    if not all(ok for _, ok in gates):
        print('NOTE: control gate failure — real-corpus verdicts not trustworthy yet')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
