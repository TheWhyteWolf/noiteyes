#!/usr/bin/env python3
"""Test the additive-walk ("cut deck") model against the isomorph catalog.

Model: c_i = D[x_i],  x_i = x_{i-1} + v(p_i)  (mod 83)
  D = secret deck order (bijection Z83 -> letters), pi = D^-1,
  v = secret per-plaintext-letter cut amounts (nonzero).

A true isomorph pair (same plaintext phrase at two places, walk offset delta)
forces  pi(B_j) = pi(A_j) + delta_pair  for every aligned position j.

All such constraints form a homogeneous linear system over GF(83) in
unknowns pi(0..82) and delta_pair.  Solution space always contains the
translation vector (pi=const, delta=0).  If the model is right and pairs
are true, the null space is exactly 2-dimensional: translation + the true
solution (up to scale).  Dimension 1 = model refuted (only trivial
solutions).  Dimension >2 = underdetermined (disconnected constraints).

Recovering pi up to affine gauge (a*pi+b) is enough: the difference stream
d_i = pi(c_i) - pi(c_{i-1}) is then an affine relabeling of v(p_i), i.e.
a monoalphabetic image of the plaintext -> attackable by frequency methods.
"""
import json
import os
import sys
from collections import defaultdict

P = 83
D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N9 = ['east-1','east-2','east-3','east-4','east-5','west-1','west-2','west-3','west-4']


def load_msgs():
    return {n: v for n, v in json.load(open(f'{D}/data/messages.json'))['messages'].items()}


def mapping_ok(w1, w2):
    f, g = {}, {}
    for a, b in zip(w1, w2):
        if f.setdefault(a, b) != b: return False
        if g.setdefault(b, a) != a: return False
    return True


def find_isomorph_pairs(S, L0=12, min_rep=3):
    def repeat_sig(w):
        first = {}; cnt = defaultdict(int)
        for v in w: cnt[v] += 1
        return tuple(first.setdefault(v, i) if cnt[v] >= 2 else -1 for i, v in enumerate(w))
    seeds = defaultdict(list)
    for n in N9:
        s = S[n]
        for i in range(len(s) - L0 + 1):
            w = s[i:i+L0]
            cnt = defaultdict(int)
            for v in w: cnt[v] += 1
            if sum(1 for v, k in cnt.items() if k >= 2) >= min_rep:
                seeds[repeat_sig(w)].append((n, i))
    raw = set()
    for sig, locs in seeds.items():
        for a in range(len(locs)):
            for b in range(a + 1, len(locs)):
                (n1, i1), (n2, i2) = locs[a], locs[b]
                if n1 == n2 and abs(i1 - i2) < L0: continue
                w1, w2 = S[n1][i1:i1+L0], S[n2][i2:i2+L0]
                if n1 != n2 and i1 == i2 and w1 == w2: continue
                if not mapping_ok(w1, w2): continue
                l1, l2, r1, r2 = i1, i2, i1 + L0, i2 + L0
                while l1 > 0 and l2 > 0 and mapping_ok(S[n1][l1-1:r1], S[n2][l2-1:r2]):
                    l1 -= 1; l2 -= 1
                while r1 < len(S[n1]) and r2 < len(S[n2]) and mapping_ok(S[n1][l1:r1+1], S[n2][l2:r2+1]):
                    r1 += 1; r2 += 1
                raw.add((n1, l1, r1, n2, l2, r2))
    maximal = []
    for p in sorted(raw, key=lambda p: -(p[2] - p[1])):
        n1, l1, r1, n2, l2, r2 = p
        if any(m[0] == n1 and m[3] == n2 and m[1] <= l1 and r1 <= m[2] and m[4] <= l2 and r2 <= m[5]
               for m in maximal):
            continue
        maximal.append(p)
    return maximal


def gauss_nullspace(rows, ncols):
    """Null space basis of homogeneous system over GF(83). rows: dict col->coef."""
    rows = [dict(r) for r in rows if r]
    pivots = {}          # col -> row index in reduced list
    reduced = []
    for r in rows:
        r = dict(r)
        for col, ri in list(pivots.items()):
            if col in r:
                factor = r[col]
                for c2, v2 in reduced[ri].items():
                    r[c2] = (r.get(c2, 0) - factor * v2) % P
                r = {c: v for c, v in r.items() if v}
        if not r:
            continue
        piv = min(r)
        inv = pow(r[piv], P - 2, P)
        r = {c: (v * inv) % P for c, v in r.items()}
        # back-substitute into existing rows
        for i, rr in enumerate(reduced):
            if piv in rr:
                factor = rr[piv]
                for c2, v2 in r.items():
                    rr[c2] = (rr.get(c2, 0) - factor * v2) % P
                reduced[i] = {c: v for c, v in rr.items() if v}
        pivots = {}
        reduced.append(r)
        for i, rr in enumerate(reduced):
            if rr:
                pivots[min(rr)] = i
        reduced = [rr for rr in reduced if rr]
        pivots = {min(rr): i for i, rr in enumerate(reduced)}
    rank = len(reduced)
    free_cols = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free_cols:
        vec = [0] * ncols
        vec[fc] = 1
        for pc, ri in pivots.items():
            vec[pc] = (-reduced[ri].get(fc, 0)) % P
        basis.append(vec)
    return rank, free_cols, basis


def main():
    S = load_msgs()
    pairs = find_isomorph_pairs(S)
    min_len = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    strong = [p for p in pairs if p[2] - p[1] >= min_len]
    print(f'isomorph pairs: {len(pairs)} maximal, {len(strong)} with L>={min_len}')

    letters = sorted({v for (n1, l1, r1, n2, l2, r2) in strong
                      for v in list(S[n1][l1:r1]) + list(S[n2][l2:r2])})
    lidx = {v: i for i, v in enumerate(letters)}
    ncols = len(letters) + len(strong)
    rows = []
    for pi_, (n1, l1, r1, n2, l2, r2) in enumerate(strong):
        for a, b in set(zip(S[n1][l1:r1], S[n2][l2:r2])):
            row = defaultdict(int)
            row[lidx[b]] = (row[lidx[b]] + 1) % P
            row[lidx[a]] = (row[lidx[a]] - 1) % P
            row[len(letters) + pi_] = P - 1        # -delta_pair
            rows.append({c: v for c, v in row.items() if v})
    print(f'letters involved: {len(letters)}/83, equations: {len(rows)}, unknowns: {ncols}')

    rank, free_cols, basis = gauss_nullspace(rows, ncols)
    print(f'rank: {rank}, null space dimension: {ncols - rank}')
    print(f'free columns: {free_cols}')

    # translation vector = (pi all-ones, delta all-zero); find a basis vector
    # not proportional to it (restricted to letter columns)
    def letters_part(vec): return vec[:len(letters)]
    indep = None
    for vec in basis:
        lp = letters_part(vec)
        if len(set(lp)) > 1:
            indep = vec
            break
    if indep is None:
        print('MODEL REFUTED at this constraint set: only constant pi solutions')
        return 1
    if ncols - rank == 2:
        print('NULL SPACE DIM == 2: unique non-trivial solution (up to affine gauge)')
    pi = {letters[i]: indep[i] for i in range(len(letters))}
    json.dump({'pi': {str(k): v for k, v in pi.items()},
               'letters_covered': letters,
               'pairs_used': strong},
              open(f'{D}/data/pi_partial.json', 'w'), indent=1)
    print(f'saved data/pi_partial.json with pi for {len(pi)} letters')

    # difference stream where both endpoints known
    known = set(pi)
    hist = defaultdict(int)
    total = have = 0
    for n in N9:
        s = S[n]
        for i in range(2, len(s)):          # skip header transition 0->1
            total += 1
            if s[i] in known and s[i-1] in known:
                have += 1
                hist[(pi[s[i]] - pi[s[i-1]]) % P] += 1
    print(f'\ndifference stream coverage: {have}/{total}')
    top = sorted(hist.items(), key=lambda kv: -kv[1])
    print('top differences:', top[:15])
    n = sum(hist.values())
    ic = sum(c * (c - 1) for c in hist.values()) / (n * (n - 1)) if n > 1 else 0
    print(f'differences observed: {n}, distinct: {len(hist)}, IoC: {ic:.5f} '
          f'(uniform={1/83:.5f}, English~0.066, Finnish~0.074)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
