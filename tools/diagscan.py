#!/usr/bin/env python3
"""Diagonal isomorph scan + deduplicated recurrence: period-4 is a duplication artifact.

Two things this does that the rest of the toolchain did not:

1. EXHAUSTIVE isomorph scan. pi_solver.find_isomorph_pairs seeds only on
   12-letter windows holding >=3 repeated letters, so sparse isomorphs (long
   stretches whose repeats are spread out) are never seeded. Here every
   message pair and every relative offset is walked as a diagonal
   x[i] ~ y[i+o]; for each start the maximal window whose letters map through
   ONE consistent bijection is taken. Evidence = number of repeats matched on
   both sides. By chance consistency breaks at the first unmatched repeat, so
   evidence >=3 is essentially never reached -- confirmed by two nulls:
     * shuffled   -- each message independently permuted
     * reversed   -- one message reversed (keeps local structure, kills
                     alignment)
   Windows that are mostly letter-identical (shared-prefix copies) are
   reported separately from genuine non-identity isomorphs.

2. DEDUPLICATED recurrence by distance. The corpus repeats material heavily:
   family A (east-1,east-2,west-1) and family B (east-4,east-5,west-4) share
   prefixes letter-for-letter, and isomorph windows repeat the same plaintext
   under a relabelling. Counting x[i]==x[i-d] naively counts one plaintext
   event once per copy. Union-find over (message, position) merges same-
   position identical letters across messages plus every aligned position of
   every evidence>=3 window; each recurrence class (class(i-d), class(i)) is
   then counted once and tested against Poisson(slots/83).

Result (see data/dedup_recurrence.json): the "period-4" excess (raw d4-ratio
2.6, 26 vs ~11.6) collapses to 11 vs 7.7, p=0.16 -- level with d5, d9, d13.
What survives is near-repeat DEPLETION: d1=0, d2 and d3 low.

Writes data/isomorphs_diag.json (all windows, evidence>=3, with a NEW flag
against data/isomorphs.json) and data/dedup_recurrence.json.
Does not modify data/isomorphs.json (other tools depend on it).
"""
import json
import math
import os
import random
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load():
    M = json.load(open(os.path.join(ROOT, 'data', 'messages.json')))['messages']
    return {k: (v if isinstance(v, list) else v['letters']) for k, v in M.items()}


def windows(x, y, o):
    """Maximal consistent-bijection windows on diagonal x[i] ~ y[i+o]."""
    out = []; last_end = -1
    lo = max(0, -o); hi = min(len(x), len(y) - o)
    for s in range(lo, hi):
        m = {}; r = {}; ev = 0; e = s
        while e < hi:
            a, b = x[e], y[e + o]
            if m.get(a, b) != b or r.get(b, a) != a:
                break
            if a in m:
                ev += 1
            m[a] = b; r[b] = a; e += 1
        if e > last_end:
            out.append((s, e, ev)); last_end = e
    return out


def scan(S, minev=5, skip_identical=True):
    """All pairs, all offsets. Returns (a,s,e,b,s2,e2,evidence,identical)."""
    names = list(S); res = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            for o in range(-len(S[a]) + 1, len(S[b])):
                for s, e, ev in windows(S[a], S[b], o):
                    if ev < minev:
                        continue
                    ident = sum(1 for k in range(s, e) if S[a][k] == S[b][k + o])
                    if skip_identical and ident >= (e - s) // 2:
                        continue
                    res.append((a, s, e, b, s + o, e + o, ev, ident))
    return res


def dedup_classes(S, minev=3):
    """Union-find over (msg,pos): identical same-position letters + evidence>=minev windows."""
    names = list(S); par = {}

    def f(x):
        while par.get(x, x) != x:
            par[x] = par.get(par[x], par[x]); x = par[x]
        return x

    def u(a, b):
        a, b = f(a), f(b)
        if a != b:
            par[a] = b
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            for p in range(min(len(S[a]), len(S[b]))):
                if S[a][p] == S[b][p]:
                    u((a, p), (b, p))
    for a, s, e, b, s2, e2, ev, _ in scan(S, minev=minev):
        for k in range(e - s):
            u((a, s + k), (b, s2 + k))
    return f


def pois_sf(k, lam):
    return 1 - sum(math.exp(-lam) * lam ** j / math.factorial(j) for j in range(k))


def dedup_recurrence(S, dmax=15, minev=3):
    """Per distance d: raw count, deduplicated count, dedup slots, Poisson expectation and p(>=)."""
    names = list(S); f = dedup_classes(S, minev)
    out = {}
    for d in range(1, dmax + 1):
        raw = sum(1 for n in names for i in range(d, len(S[n])) if S[n][i] == S[n][i - d])
        ev = {(f((n, i - d)), f((n, i))) for n in names for i in range(d, len(S[n]))
              if S[n][i] == S[n][i - d]}
        slots = {(f((n, i - d)), f((n, i))) for n in names for i in range(d, len(S[n]))}
        lam = len(slots) / 83
        out[d] = {'raw': raw, 'dedup': len(ev), 'slots': len(slots),
                  'expected': round(lam, 2), 'p_ge': round(pois_sf(len(ev), lam), 4)}
    return out


def dedup_d4_ratio(S, minev=3):
    """Deduplicated analogue of tournament's d4_ratio: d4 / mean(d3, d5), class-counted."""
    r = dedup_recurrence(S, dmax=5, minev=minev)
    ref = (r[3]['dedup'] + r[5]['dedup']) / 2 or 1
    return round(r[4]['dedup'] / ref, 2)


def catalog_flag(r, cat):
    a, s, e, b, s2, e2 = r[:6]
    for c in cat:
        for (x, xs, xe, y, ys, ye) in (c, (c[3], c[4], c[5], c[0], c[1], c[2])):
            if x == a and y == b and s2 - s == ys - xs and max(s, xs) < min(e, xe):
                return 'catalog'
    return 'NEW'


def main():
    S = load(); names = list(S)
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    real = scan(S, minev=3)
    P = json.load(open(os.path.join(ROOT, 'data', 'isomorphs.json')))['pairs']
    cat = [(p['a']['msg'], p['a']['start'], p['a']['end'],
            p['b']['msg'], p['b']['start'], p['b']['end']) for p in P]
    rows = []
    for r in sorted(real, key=lambda r: (-r[6], r[0], r[1])):
        flag = catalog_flag(r, cat)
        rows.append({'a': {'msg': r[0], 'start': r[1], 'end': r[2]},
                     'b': {'msg': r[3], 'start': r[4], 'end': r[5]},
                     'length': r[2] - r[1], 'evidence': r[6], 'identical': r[7],
                     'status': flag})
        print(f'{r[0]}[{r[1]}:{r[2]}] ~ {r[3]}[{r[4]}:{r[5]}]  len {r[2]-r[1]}  ev {r[6]}  {flag}')

    # nulls
    rng = random.Random(1); sh = []
    for _ in range(trials):
        T = {k: rng.sample(v, len(v)) for k, v in S.items()}
        sh.append(Counter(r[6] for r in scan(T, minev=3)))
    rev = []
    for n in names:
        T = dict(S); T[n] = S[n][::-1]
        rev.append(sum(1 for r in scan(T, minev=3) if n in (r[0], r[3])))
    nulls = {}
    for thr in (3, 5, 7):
        rv = sum(1 for r in real if r[6] >= thr)
        nv = [sum(v for k, v in c.items() if k >= thr) for c in sh]
        nulls[thr] = {'real': rv, 'shuffled_mean': round(sum(nv) / len(nv), 2),
                      'shuffled_max': max(nv)}
        print(f'ev>={thr}: real {rv}  shuffled mean {nulls[thr]["shuffled_mean"]} max {max(nv)}')
    print(f'reversed-message null: windows ev>=3 involving the reversed message: {sum(rev)}')

    rec = dedup_recurrence(S)
    print('\n d  raw  dedup  expected  p(>=)')
    for d, v in rec.items():
        print(f'{d:2d} {v["raw"]:4d} {v["dedup"]:6d} {v["expected"]:9.1f}  {v["p_ge"]:.3f}')
    print('dedup d4_ratio', dedup_d4_ratio(S))

    json.dump({'finder': 'diagscan.scan (all pairs, all offsets, minev=3, non-identical)',
               'nulls': {'shuffled_trials': trials, 'by_threshold': nulls,
                         'reversed_total_ev3': sum(rev)},
               'pairs': rows},
              open(os.path.join(ROOT, 'data', 'isomorphs_diag.json'), 'w'), indent=1)
    json.dump({'method': 'union-find dedup (same-position identical + diagscan ev>=3 windows); '
                         'Poisson(slots/83)',
               'd4_ratio_raw': round(rec[4]['raw'] / ((rec[3]['raw'] + rec[5]['raw']) / 2), 2),
               'd4_ratio_dedup': dedup_d4_ratio(S), 'by_distance': rec},
              open(os.path.join(ROOT, 'data', 'dedup_recurrence.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
