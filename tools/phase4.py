#!/usr/bin/env python3
"""Step A: is the confirmed period-4 keystream cycle a KEY INDEX that reduces
the cipher, in ADDITIVE form?

The only period-4 form compatible with the observed L=33 exact isomorphs (a
full re-key every 4 is impossible) is a single deck D with a period-4 additive
offset on the state walk:
    c_i = D[x_i],  x_i = x_{i-1} + v(p_i) + k[i mod 4]   (mod 83)
For two occurrences of a phrase at relative offset delta = (l2 - l1):
  * delta == 0 mod 4 -> the k-terms cancel termwise, so the pair is a PURE
    translation (sigma = D-translation).
  * delta != 0 mod 4 -> the k-difference oscillates with period 4 and never
    cancels; exact isomorphy then REQUIRES every internal repeat to sit at a
    position-pair that is 0 mod 4 (only same-phase positions can co-repeat).

This yields three predictions of the additive-period-4 model, each checked
against synthetic corpora of known construction before the real read:

  S1 delta-distribution  additive-p4 SUPPRESSES delta!=0 isomorphs (they can
     only survive with fully mod-4-aligned repeats), so pairs concentrate at
     delta==0. A plain walk spreads pairs over all delta.
  S2 repeat alignment    additive-p4 forces delta!=0 pairs to ~1.0 mod-4
     repeat alignment; a plain walk sits at chance ~0.25.
  S3 stratified solve     under additive-p4 the delta==0 segments form a
     consistent, injective translation system recovering D, even when the
     POOLED system is dim-1.

Gates (must pass every run, else the real verdict is not read):
  period4 walk -> S1 concentrates at delta==0; S3 delta==0 recovers the deck
  plain walk   -> S1 spread; S2 chance ~0.25; S3 pooled cracks whole
  shuffle      -> S3 no delta-group cracks (noise / over-subset guard)
"""
import json
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
from pi_solver import P, find_isomorph_pairs, load_msgs
from sigma_web import (MINLEN, PA, ROOT, WEIGHTS, core_blocks, pi_matches_truth,
                       pinned_eval, solve_blocks, witness_trim)

SEED = 1


def delta4(block):
    return (block[3] - block[1]) % 4


def get_cores(S, names):
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    joint = witness_trim(S, raw)
    cores, tags = core_blocks(joint)
    return raw, cores, tags


# ------------------------------------------------------------ signatures

def delta_distribution(raw):
    """S1: how isomorph pairs spread over delta mod 4. additive-p4 piles them
    at 0; a plain walk spreads them."""
    d = Counter((l2 - l1) % 4 for (n1, l1, r1, n2, l2, r2) in raw)
    n = sum(d.values())
    frac0 = d.get(0, 0) / n if n else 0
    return {'counts': dict(sorted(d.items())), 'n': n, 'frac_delta0': round(frac0, 3)}


def repeat_alignment(S, raw):
    """S2: mean fraction of a pair's internal repeats whose position-gap is
    0 mod 4, grouped by delta. additive-p4 requires ~1.0 for delta!=0."""
    def frac_aligned(seq):
        pos = defaultdict(list)
        for j, x in enumerate(seq):
            pos[x].append(j)
        prs = [(a, b) for ix in pos.values() if len(ix) >= 2
               for i, a in enumerate(ix) for b in ix[i + 1:]]
        if not prs:
            return None
        return sum((b - a) % 4 == 0 for a, b in prs), len(prs)
    byd = defaultdict(lambda: [0, 0])
    for n1, l1, r1, n2, l2, r2 in raw:
        fa = frac_aligned(S[n1][l1:r1])
        if fa:
            byd[(l2 - l1) % 4][0] += fa[0]
            byd[(l2 - l1) % 4][1] += fa[1]
    return {d: {'aligned': a, 'total': t, 'frac': round(a / t, 3)}
            for d, (a, t) in sorted(byd.items()) if t}


# ------------------------------------------------------- stratified solve

def crack_ok(pin):
    return bool(pin and pin['n_pinned'] >= 12
                and pin['distinct'] >= 0.95 * pin['n_pinned']
                and pin['ioc'][0] > 0.05)


def solve_report(S, names, blocks, rng, truth=None):
    if len(blocks) < 2:
        return {'nseg': len(blocks), 'dim': None, 'distinct': 0, 'pinned': 0,
                'ioc': None, 'crack': False, 'truth': None}
    r = solve_blocks(S, blocks, rng, tries=60)
    pin = pinned_eval(S, names, r) if r['dim'] >= 2 else None
    out = {'nseg': len(blocks), 'dim': r['dim'],
           'pinned': pin['n_pinned'] if pin else 0,
           'distinct': pin['distinct'] if pin else 0,
           'ioc': pin['ioc'] if pin else None, 'crack': crack_ok(pin),
           'truth': round(pi_matches_truth(pin['pi'], truth['deck']), 3)
           if (truth and pin) else None}
    return out


def stratify(S, names, raw, cores, rng, truth=None, label=''):
    print(f'\n== {label} ==')
    s1 = delta_distribution(raw)
    s2 = repeat_alignment(S, raw)
    print(f'S1 pair delta mod4: {s1["counts"]}  (frac@0 = {s1["frac_delta0"]})')
    print('S2 repeat alignment by delta: '
          + '  '.join(f'd{d}:{v["frac"]}({v["total"]})' for d, v in s2.items()))
    out = {'label': label, 's1_delta_dist': s1, 's2_alignment': s2}
    pooled = solve_report(S, names, cores, rng, truth)
    out['pooled'] = pooled
    print(f"S3 POOLED : segs {pooled['nseg']} dim {pooled['dim']} "
          f"pinned {pooled['pinned']} distinct {pooled['distinct']} "
          f"IoC {pooled['ioc']} crack={pooled['crack']} truth={pooled['truth']}")
    out['strata'] = {}
    for d in range(4):
        bl = [b for b in cores if delta4(b) == d]
        rep = solve_report(S, names, bl, rng, truth)
        out['strata'][d] = rep
        print(f"   delta={d}: segs {rep['nseg']} dim {rep['dim']} "
              f"pinned {rep['pinned']} distinct {rep['distinct']} "
              f"IoC {rep['ioc']} crack={rep['crack']} truth={rep['truth']}")
    return out


# ------------------------------------------------------ synthetic corpora

REAL_LENS = [99, 118, 137, 119, 114, 103, 102, 124, 120]


def rich_phrase_k(rng, v, k, L):
    """A phrase whose walk image -- INCLUDING the period-4 offset k, with the
    phrase inserted at an offset == 0 mod 4 -- is repeat-dense enough to seed
    the finder (a 12-window with >=3 internally repeated states)."""
    while True:
        ph = ''.join(rng.choices(PA, weights=WEIGHTS, k=L))
        pre = [0]
        for idx, ch in enumerate(ph):
            pre.append((pre[-1] + v[ch] + k[(idx + 1) % 4]) % P)
        by = defaultdict(list)
        for idx in range(1, L + 1):
            by[pre[idx]].append(idx - 1)
        reps = [ix for ix in by.values() if len(ix) >= 2]
        if sum(len(ix) * (len(ix) - 1) // 2 for ix in reps) < 5:
            continue
        for w0 in range(0, L - 11):
            if sum(1 for ix in reps
                   if len([i for i in ix if w0 <= i < w0 + 12]) >= 2) >= 3:
                return ph


def gen(mode, seed):
    rng = random.Random(seed)
    deck = list(range(P))
    rng.shuffle(deck)
    v = {ch: rng.randrange(1, P) for ch in PA}
    k = ([0, rng.randrange(1, P), rng.randrange(1, P), rng.randrange(1, P)]
         if mode == 'period4' else [0, 0, 0, 0])
    names = [f'sim-{i}' for i in range(9)]
    plain = {n: list(''.join(rng.choices(PA, weights=WEIGHTS, k=L)))
             for n, L in zip(names, REAL_LENS)}
    rp = (lambda L: rich_phrase_k(rng, v, k, L))
    ph1, ph2, ph3 = rp(30), rp(34), rp(28)

    def ins(n, ph, pos):
        plain[n][pos:pos + len(ph)] = list(ph)
    if mode == 'period4':                     # every offset == 0 mod 4
        ins('sim-0', ph1, 32); ins('sim-1', ph1, 60); ins('sim-5', ph1, 28)
        ins('sim-2', ph2, 40); ins('sim-4', ph2, 72); ins('sim-3', ph2, 36)
        ins('sim-6', ph3, 44); ins('sim-8', ph3, 60); ins('sim-7', ph3, 32)
    else:                                     # offsets spread over mod 4
        ins('sim-0', ph1, 30); ins('sim-1', ph1, 35); ins('sim-5', ph1, 61)
        ins('sim-2', ph2, 40); ins('sim-4', ph2, 33); ins('sim-3', ph2, 46)
        ins('sim-6', ph3, 44); ins('sim-8', ph3, 33); ins('sim-7', ph3, 51)
    S = {}
    for n in names:
        plain[n][0:2] = list('qk')
        x = 17
        out = [deck[x]]
        for i, ch in enumerate(plain[n]):
            x = (x + v[ch] + k[(i + 1) % 4]) % P
            out.append(deck[x])
        S[n] = out
    return S, {'deck': deck, 'k': k}


def gate(name, ok):
    print(f'  {"PASS" if ok else "FAIL":4s}  {name}')
    return ok


def main():
    rng = random.Random(SEED)
    report, gates = {}, []

    Sp, tp = gen('period4', SEED)
    raw, cores, _ = get_cores(Sp, sorted(Sp))
    r = stratify(Sp, sorted(Sp), raw, cores, rng, truth=tp,
                 label='control: period4 walk')
    report['period4'] = r
    gates.append(gate('period4: S1 pairs concentrate at delta==0 (frac>=0.8)',
                      r['s1_delta_dist']['frac_delta0'] >= 0.8))
    d0 = r['strata'][0]
    gates.append(gate('period4: S3 delta==0 recovers the deck (truth>=0.9)',
                      d0['truth'] is not None and d0['truth'] >= 0.9))

    Spl, tpl = gen('plain', SEED)
    raw, cores, _ = get_cores(Spl, sorted(Spl))
    r = stratify(Spl, sorted(Spl), raw, cores, rng, truth=tpl,
                 label='control: plain walk')
    report['plain'] = r
    gates.append(gate('plain: S1 spread (frac@0 < 0.6)',
                      r['s1_delta_dist']['frac_delta0'] < 0.6))
    a2 = r['s2_alignment']
    ali_ne = [v['frac'] for d, v in a2.items() if d != 0]
    gates.append(gate('plain: S2 delta!=0 alignment ~ chance (all < 0.55)',
                      all(f < 0.55 for f in ali_ne) if ali_ne else False))
    gates.append(gate('plain: S3 pooled cracks whole (truth>=0.9)',
                      r['pooled']['crack'] and r['pooled']['truth']
                      and r['pooled']['truth'] >= 0.9))

    Ssh = {}
    for n in sorted(Sp):
        s = list(Sp[n]); rng.shuffle(s); Ssh[n] = s
    raw, cores, _ = get_cores(Ssh, sorted(Ssh))
    anycrack = False
    if len(cores) >= 2:
        r = stratify(Ssh, sorted(Ssh), raw, cores, rng, label='control: shuffle')
        report['shuffle'] = r
        anycrack = r['pooled']['crack'] or any(s['crack'] for s in r['strata'].values())
    else:
        report['shuffle'] = {'cores': len(cores)}
        print('\n== control: shuffle ==\ncores < 2 (as expected)')
    gates.append(gate('shuffle: no delta-group cracks (noise guard)', not anycrack))

    print('\n== control gates ==')
    allok = all(gates)
    if not allok:
        print('  control gate failure -> real verdict NOT trustworthy')

    S = load_msgs()
    names = pi_solver.N9
    raw, cores, _ = get_cores(S, names)
    report['real'] = stratify(S, names, raw, cores, rng, label='REAL CORPUS')

    json.dump(report, open(f'{ROOT}/data/phase4_report.json', 'w'),
              indent=1, default=str)
    print('\nwrote data/phase4_report.json')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
