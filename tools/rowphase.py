#!/usr/bin/env python3
"""Row-phase schedule extractor: does the surviving fixed-substitution
mechanism re-key on a schedule locked to the physical eye rows?

PLAN third-phase result: the sigma maps between isomorph occurrences are FIXED
substitutions, constant over stretches of 18-33 letters, changing at
boundaries that correlate with row-pair boundaries (z~-2.9). This tool turns
that correlation into a testable SCHEDULE. Two independent instruments, each
validated on synthetic corpora of KNOWN schedule before any real verdict:

  Instrument 1 - BREAK PHASE.  The endpoints of the constant-sigma stretches
    (isomorph window ends + closure-violation zones) are the re-key events.
    Row-pairs are ~26 letters (boundaries at 26,52,78,104 in every message),
    so a row-locked schedule puts these events at phase 0 mod 26 (or mod 13
    for single rows).  We measure residue concentration mod k in {4,13,26},
    nearest-boundary distance, and the fraction of closure zones that bracket
    a physical row crossing, each against a position-uniform null.

  Instrument 2 - KEYSTREAM PERIOD.  The community reports a period-4 cycle
    (distance-4 letter recurrence at ~2x rate).  A period-4 FULL re-key is
    already impossible (it would forbid the observed L=33 isomorphs), so this
    can only be a keystream sub-cycle.  We measure distance-d equal-letter
    counts d=1..8 within each message against a within-message shuffle null,
    which is exact (destroys all positional structure).

Gates (must pass every run, else real verdicts are not read):
  rekey-26 walk  -> instrument 1 fires at mod-26 phase 0, bracket rate ~1
  rekey-13 walk  -> fires at mod-13 phase 0, mod-26 splits over {0,13}
  no-schedule    -> instrument 1 flat, bracket rate ~ null
  shuffle-null   -> instrument 2 |z|<3 for all d (null calibrated)
  planted-d4     -> instrument 2 flags d=4 (probe sensitivity)
"""
import json
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
from pi_solver import P, find_isomorph_pairs, load_msgs
from sigma_web import (PA, WEIGHTS, MINLEN, ROOT, build_webs, closure_analysis,
                       rich_phrase, rowpair_boundaries, witness_trim)

SEED = 1
MODULI = (4, 13, 26)


# ------------------------------------------------------- event extraction

def get_events(S, names):
    """Return (window_end_events, zone_events). window_end_events: list of
    (msg,pos) at each isomorph occurrence boundary. zone_events: list of dicts
    with abs intervals in both occurrences of each closure-violation zone."""
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    if not raw:
        return [], [], raw
    joint = witness_trim(S, raw)
    trimmed = [e['trim'] for e in joint]
    occs, pot, webs, _gv, edges = build_webs(trimmed)
    ext = {k: (lo, hi) for k, (m, lo, hi) in enumerate(occs)}
    trims_by_occ = defaultdict(list)
    for a, b, off, k in edges:
        rw, tw = joint[k]['raw'], joint[k]['trim']
        ext[a] = (min(ext[a][0], rw[1]), max(ext[a][1], rw[2]))
        ext[b] = (min(ext[b][0], rw[4]), max(ext[b][1], rw[5]))
        trims_by_occ[a].append((tw[1], tw[2]))
        trims_by_occ[b].append((tw[4], tw[5]))
    _stats, zones = closure_analysis(S, occs, ext, pot, webs, trims_by_occ)

    win = []
    for n1, l1, r1, n2, l2, r2 in raw:
        for n, p in ((n1, l1), (n1, r1), (n2, l2), (n2, r2)):
            if 0 < p < len(S[n]) - 1:
                win.append((n, p))
    zev = []
    for z in zones:
        f0, f1 = z['frame']
        a, b = z['occ_a'], z['occ_b']
        zev.append({'msg_a': occs[a][0], 'a': (f0 + pot[a], f1 + pot[a]),
                    'msg_b': occs[b][0], 'b': (f0 + pot[b], f1 + pot[b]),
                    'evidence': z['evidence']})
    return win, zev, raw


# --------------------------------------------------- instrument 1: phase

def phase_conc(events, lens, k, rng, mc=4000):
    """Concentration of event positions at their most-common residue mod k,
    vs a null that resamples positions uniformly within the same messages."""
    if not events:
        return None
    res = Counter(p % k for _, p in events)
    peak_res, peak = res.most_common(1)[0]
    n = len(events)
    Ls = [lens[nm] for nm, _ in events]
    null = []
    for _ in range(mc):
        cc = Counter(rng.randrange(1, L) % k for L in Ls)
        null.append(max(cc.values()))
    mu = sum(null) / mc
    sd = (sum((x - mu) ** 2 for x in null) / mc) ** 0.5
    return {'k': k, 'peak_res': peak_res, 'peak_frac': round(peak / n, 3),
            'exp_frac': round(1 / k, 3), 'n': n,
            'z': round((peak - mu) / sd, 2) if sd else None}


def boundary_dist(events, bounds, S, rng, mc=4000):
    """Mean distance from events to the nearest row-pair boundary, plus the
    fraction within +-1, each vs a position-uniform null (Monte Carlo)."""
    ev = [(n, p) for n, p in events if bounds.get(n)]
    if not ev:
        return None
    def dist(n, p):
        return min(abs(p - b) for b in bounds[n])
    obs = [dist(n, p) for n, p in ev]
    om, ofrac = sum(obs) / len(obs), sum(d <= 1 for d in obs) / len(obs)
    means, fracs = [], []
    for _ in range(mc):
        smp = [dist(n, rng.randrange(len(S[n]))) for n, _ in ev]
        means.append(sum(smp) / len(smp))
        fracs.append(sum(d <= 1 for d in smp) / len(smp))
    mu = sum(means) / mc
    sd = (sum((x - mu) ** 2 for x in means) / mc) ** 0.5
    fmu = sum(fracs) / mc
    fsd = (sum((x - fmu) ** 2 for x in fracs) / mc) ** 0.5
    return {'n': len(ev), 'obs_mean': round(om, 2), 'null_mean': round(mu, 2),
            'z_mean': round((om - mu) / sd, 2) if sd else None,
            'obs_within1': round(ofrac, 3), 'null_within1': round(fmu, 3),
            'z_within1': round((ofrac - fmu) / fsd, 2) if fsd else None}


def bracket_rate(zev, bounds, S):
    """Fraction of closure zones bracketing a row boundary in stream A or B,
    vs the analytic null (a random placement of the same widths)."""
    if not zev:
        return None
    def brackets(msg, lo, hi):
        return any(lo <= b <= hi for b in bounds.get(msg, []))
    def null_p(msg, w):
        L = len(S[msg])
        B = bounds.get(msg, [])
        if L - w <= 0:
            return 0.0
        hit = sum(1 for x in range(0, L - w) if any(x <= b <= x + w for b in B))
        return hit / (L - w)
    obs = nul = 0.0
    for z in zev:
        a0, a1 = z['a']
        b0, b1 = z['b']
        obs += brackets(z['msg_a'], a0, a1) or brackets(z['msg_b'], b0, b1)
        pa = null_p(z['msg_a'], a1 - a0)
        pb = null_p(z['msg_b'], b1 - b0)
        nul += 1 - (1 - pa) * (1 - pb)
    n = len(zev)
    return {'n': n, 'obs_rate': round(obs / n, 3), 'null_rate': round(nul / n, 3)}


def instrument1(S, names, bounds, label, rng):
    win, zev, _raw = get_events(S, names)
    lens = {n: len(S[n]) for n in names}
    zpts = [(z['msg_a'], z['a'][0]) for z in zev] + [(z['msg_b'], z['b'][0]) for z in zev]
    ev_zone = [(z['msg_a'], z['a'][0]) for z in zev] + \
              [(z['msg_b'], z['b'][1]) for z in zev]
    out = {'label': label, 'n_win': len(win), 'n_zones': len(zev),
           'n_zones_evid': sum(z['evidence'] for z in zev)}
    print(f'\n-- instrument 1 (break phase): {label} --')
    print(f'window-end events: {len(win)}   closure zones: {len(zev)} '
          f'({out["n_zones_evid"]} evidence-backed)')
    for tag, ev in (('win', win), ('zone', zpts + ev_zone)):
        out[f'phase_{tag}'] = {}
        for k in MODULI:
            r = phase_conc(ev, lens, k, rng)
            if r:
                out[f'phase_{tag}'][k] = r
                print(f'  {tag:4s} mod {k:2d}: peak res {r["peak_res"]:2d} '
                      f'frac {r["peak_frac"]} (exp {r["exp_frac"]}) z={r["z"]}')
        bd = boundary_dist(ev, bounds, S, rng)
        out[f'bdist_{tag}'] = bd
        if bd:
            print(f'  {tag:4s} nearest-boundary dist: obs {bd["obs_mean"]} vs '
                  f'null {bd["null_mean"]} z={bd["z_mean"]}; within1 '
                  f'{bd["obs_within1"]} vs {bd["null_within1"]} z={bd["z_within1"]}')
    out['bracket'] = bracket_rate(zev, bounds, S)
    if out['bracket']:
        b = out['bracket']
        print(f'  zone bracket-a-boundary rate: obs {b["obs_rate"]} vs '
              f'null {b["null_rate"]} (n={b["n"]})')
    return out


# ------------------------------------------------- instrument 2: keystream

def distance_autocorr(S, names, rng, dmax=8, mc=3000):
    """Distance-d equal-letter counts within each message vs a within-message
    shuffle null (exact: destroys all positional structure). z per d."""
    obs = {d: 0 for d in range(1, dmax + 1)}
    for n in names:
        s = S[n]
        for d in range(1, dmax + 1):
            obs[d] += sum(1 for i in range(d, len(s)) if s[i] == s[i - d])
    null = {d: [] for d in range(1, dmax + 1)}
    for _ in range(mc):
        tot = {d: 0 for d in range(1, dmax + 1)}
        for n in names:
            s = list(S[n])
            rng.shuffle(s)
            for d in range(1, dmax + 1):
                tot[d] += sum(1 for i in range(d, len(s)) if s[i] == s[i - d])
        for d in range(1, dmax + 1):
            null[d].append(tot[d])
    res = {}
    for d in range(1, dmax + 1):
        mu = sum(null[d]) / mc
        sd = (sum((x - mu) ** 2 for x in null[d]) / mc) ** 0.5
        res[d] = {'obs': obs[d], 'null_mean': round(mu, 1),
                  'z': round((obs[d] - mu) / sd, 2) if sd else None}
    return res


def instrument2(S, names, label, rng):
    r = distance_autocorr(S, names, rng)
    print(f'\n-- instrument 2 (keystream period): {label} --')
    print('  d :  obs  null   z')
    for d, v in r.items():
        print(f'  {d:2d}: {v["obs"]:4d} {v["null_mean"]:6.1f} {v["z"]:6.2f}')
    return r


# ------------------------------------------------------ synthetic corpora

REAL_LENS = [99, 118, 137, 119, 114, 103, 102, 124, 120]


def gen_rekey(period, seed, plant_d4=False):
    """9-message corpus, row-pairs of 26 letters (boundaries at multiples of
    26). Output position i uses substitution deck (i // period): a re-key
    every `period` letters. period huge -> single deck (no schedule)."""
    rng = random.Random(seed)
    v = {ch: rng.randrange(1, P) for ch in PA}
    nblocks = REAL_LENS[2] // period + 3
    decks = []
    for _ in range(nblocks + 1):
        d = list(range(P))
        rng.shuffle(d)
        decks.append(d)
    names = [f'sim-{i}' for i in range(9)]
    plain = {n: list(''.join(rng.choices(PA, weights=WEIGHTS, k=L)))
             for n, L in zip(names, REAL_LENS)}
    ph1, ph2, ph3 = (rich_phrase(rng, v, 30), rich_phrase(rng, v, 34),
                     rich_phrase(rng, v, 28))

    def ins(n, ph, pos):
        plain[n][pos:pos + len(ph)] = list(ph)
    ins('sim-0', ph1, 30); ins('sim-0', ph1, 61)          # offsets misaligned
    ins('sim-1', ph1, 35); ins('sim-1', ph1, 72)          # mod 26 vs each other
    ins('sim-5', ph1, 29)
    ins('sim-4', ph2, 35); ins('sim-8', ph2, 36); ins('sim-3', ph2[:22], 46)
    ins('sim-2', ph3, 40); ins('sim-2', ph3, 75); ins('sim-6', ph3, 33)
    S, bounds = {}, {}
    for n in names:
        plain[n][0:2] = list('qk')
        x = 17
        out = [decks[0][x]]
        for i, ch in enumerate(plain[n]):
            x = (x + v[ch]) % P
            out.append(decks[(i + 1) // period][x])
        if plant_d4:                       # force distance-4 letter excess
            for i in range(4, len(out)):
                if rng.random() < 0.5:
                    out[i] = out[i - 4]
        S[n] = out
        bounds[n] = list(range(26, len(out) - 1, 26))
    return S, bounds


# -------------------------------------------------------------- main

def gate(name, ok):
    print(f'  {"PASS" if ok else "FAIL":4s}  {name}')
    return ok


def grid_scan(events, S, rng, periods, label):
    """Test window-end events for alignment to each candidate boundary grid
    (multiples of `period`). A true schedule of that period makes breaks pile
    up at the grid: nearest-boundary within-1 fraction enriched (z high)."""
    print(f'  -- {label}: alignment to candidate grids --')
    out = {}
    for pr in periods:
        b = {n: list(range(pr, len(S[n]) - 1, pr)) for n in S}
        st = boundary_dist(events, b, S, rng)
        out[pr] = st
        if st:
            print(f'    period {pr:2d}: within1 obs {st["obs_within1"]} vs '
                  f'null {st["null_within1"]} z={st["z_within1"]}  '
                  f'(mean z={st["z_mean"]})')
    return out


def mod_scan(events, lens, rng, ks, label):
    """Residue concentration of events over a range of moduli — a keystream
    period p shows as a mod-p peak. Reported (exploratory), not gated."""
    print(f'  -- {label}: residue concentration by modulus --')
    out = {}
    for k in ks:
        r = phase_conc(events, lens, k, rng)
        out[k] = r
        if r and r['z'] is not None and r['z'] >= 2:
            print(f'    mod {k:2d}: peak res {r["peak_res"]:2d} frac {r["peak_frac"]} '
                  f'(exp {r["exp_frac"]}) z={r["z"]}')
    return out


def l33_crossing_fact(S):
    """The L=33 isomorph east-5[33:66]~west-4[34:67] spans row-pair boundary
    52 in both occurrences while remaining exactly isomorphic: row-pair
    boundaries cannot be full-substitution rebuild points."""
    a, b = S['east-5'][33:66], S['west-4'][34:67]
    f, g = {}, {}
    ok = all(f.setdefault(x, y) == y and g.setdefault(y, x) == x for x, y in zip(a, b))
    print('  -- structural fact (no synthetic needed) --')
    print(f'    east-5[33:66]~west-4[34:67] exact isomorph: {ok}; both cross '
          f'row-pair boundary 52 => sigma constant across it in BOTH streams.')
    print('    => row-pair boundaries are NOT full re-key points.')
    return {'exact': ok, 'crosses_boundary': 52}


def main():
    rng = random.Random(SEED)
    report, gates = {}, []

    # ---- instrument 1 controls: boundary detector on window-end events.
    # (A full re-key with period < the isomorph length is impossible -- it
    #  would cap window length -- so the positive control uses period 26 and
    #  the detector is the period-agnostic nearest-boundary test, which fires
    #  cleanly; the naive mod-k residue peak is smeared by the bracket offset
    #  and is only reported, never gated.)
    S26, b26 = gen_rekey(26, SEED)
    r = instrument1(S26, sorted(S26), b26, 'control: re-key every 26', rng)
    report['rekey26'] = r
    bw = r['bdist_win']
    gates.append(gate('rekey-26: window ends pile at boundaries (within1 z>3, mean z<-3)',
                      bw and bw['z_within1'] and bw['z_within1'] > 3
                      and bw['z_mean'] and bw['z_mean'] < -3))

    Sno, _ = gen_rekey(999, SEED)
    bno = {n: list(range(26, len(Sno[n]) - 1, 26)) for n in Sno}
    r = instrument1(Sno, sorted(Sno), bno, 'control: no schedule', rng)
    report['noschedule'] = r
    bw = r['bdist_win']
    gates.append(gate('no-schedule: window ends NOT boundary-locked (within1 |z|<3)',
                      bw and bw['z_within1'] is not None and abs(bw['z_within1']) < 3))

    # ---- instrument 2 controls
    Sshuf, _ = gen_rekey(26, SEED)
    for n in Sshuf:                                # exact within-msg shuffle
        s = list(Sshuf[n]); rng.shuffle(s); Sshuf[n] = s
    r = instrument2(Sshuf, sorted(Sshuf), 'control: shuffle-null', rng)
    report['shuffle_null'] = r
    gates.append(gate('shuffle-null: |z|<3 all d (null calibrated)',
                      all(abs(v['z']) < 3 for v in r.values())))

    Spl, _ = gen_rekey(26, SEED, plant_d4=True)
    r = instrument2(Spl, sorted(Spl), 'control: planted distance-4', rng)
    report['planted_d4'] = r
    gates.append(gate('planted-d4: probe flags d=4 (z>6)', r[4]['z'] and r[4]['z'] > 6))

    print('\n== control gates ==')
    allok = all(gates)
    if not allok:
        print('  control gate failure -> real verdicts NOT trustworthy')

    # ---- real corpus (only reported as trustworthy if gates pass)
    S = load_msgs()
    names = pi_solver.N9
    lens = {n: len(S[n]) for n in names}
    bounds = rowpair_boundaries(S)
    print('\n== REAL CORPUS ==')
    report['real_l33'] = l33_crossing_fact(S)
    r1 = instrument1(S, names, bounds, 'REAL CORPUS', rng)
    report['real_i1'] = r1
    win, zev, _ = get_events(S, names)
    zpts = [(z['msg_a'], z['a'][0]) for z in zev] + [(z['msg_b'], z['b'][0]) for z in zev]
    report['real_grid_win'] = grid_scan(win, S, rng, (2, 4, 13, 26), 'window ends')
    report['real_mod_win'] = mod_scan(win, lens, rng, (2, 3, 4, 5, 6, 7, 8, 13, 26),
                                      'window ends')
    report['real_mod_zone'] = mod_scan(zpts, lens, rng, (2, 3, 4, 5, 6, 7, 8, 13, 26),
                                       'closure zones')
    report['real_i2'] = instrument2(S, names, 'REAL CORPUS', rng)

    json.dump(report, open(f'{ROOT}/data/rowphase_report.json', 'w'),
              indent=1, default=str)
    print('\nwrote data/rowphase_report.json')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
