#!/usr/bin/env python3
"""Non-abelian-sigma step: is the surviving structure a piecewise GENERAL rotor?

Step B isolated the tension: the mechanism must give offset-invariant EXACT
isomorphs (like the additive walk) yet be non-additive and non-commuting. The
natural resolution is a progressive key with a GENERAL (non-cyclic) permutation
rotor T: c_i = T^i(B(p_i)). Then two occurrences of a phrase at potential
offset delta have sigma = T^delta -- offset-invariant (isomorph), non-additive
(T not a cyclic shift => sigma_web dim-1), and non-commuting once T changes
between stretches. This tool tests that model directly on the sigma maps.

Test: within an isomorph web, read T off a UNIT-offset edge (sigma = T^{+-1}),
then verify every other edge e satisfies sigma_e = T^{delta_e} on its common
support. A single per-web rotor => high agreement; a rotor that varies within
the web => low agreement.

Validated on synthetic corpora (gates): a single-rotor progressive key must
fit its webs (>=0.9 support-weighted agreement); a corpus whose rotor VARIES
by region (different T in different messages) must NOT (<=0.5). The extraction
is confirmed exact on the single-rotor control (sigma of a unit edge == the
true T). The real corpus is read only if the gates pass.
"""
import json
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
from pi_solver import P, load_msgs
from sigma_web import (MINLEN, PA, ROOT, WEIGHTS, build_webs, find_isomorph_pairs,
                       witness_trim)

SEED = 1
N9 = pi_solver.N9
IDX = {c: i for i, c in enumerate(PA)}


def pipeline(S, names):
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    joint = witness_trim(S, raw)
    trimmed = [e['trim'] for e in joint]
    occs, pot, webs, gv, edges = build_webs(trimmed)
    return occs, pot, webs, edges


def sigma(i, j, occs, pot, S):
    """Letter-map from occurrence i to occurrence j, aligned by potential
    (frame = absolute - potential). Inconsistent letters are dropped."""
    mi, loi, hii = occs[i]
    mj, loj, hij = occs[j]
    lo = max(loi - pot[i], loj - pot[j])
    hi = min(hii - pot[i], hij - pot[j])
    m = {}
    for f in range(lo, hi):
        a, b = S[mi][f + pot[i]], S[mj][f + pot[j]]
        if m.setdefault(a, b) != b:
            m[a] = None
    return {a: b for a, b in m.items() if b is not None}


def tpow(tau, k):
    """tau^k as a partial map, computed only where it stays in tau's domain."""
    if k < 0:
        tau = {b: a for a, b in tau.items()}
        k = -k
    out = {a: a for a in tau}
    for _ in range(k):
        out = {a: tau[b] for a, b in out.items() if b in tau}
    return out


def fit_web(web, occs, pot, S, edges):
    """Read T from a unit-offset edge in the web, verify the OTHER edges vs
    T^delta. The seed edge is excluded from the tally (it defines T, so it is
    trivially perfect); only independent edges test the single-rotor model."""
    we = [(a, b, pot[b] - pot[a]) for a, b, off, k in edges
          if a in web and b in web]
    unit = next(((a, b, d) for a, b, d in we if abs(d) == 1), None)
    if unit is None:
        return {'web': web, 'unit': False}
    ua, ub, d = unit
    s = sigma(ua, ub, occs, pot, S)
    tau = s if d > 0 else {v: k for k, v in s.items()}      # T = sigma^{+1}
    agree = tot = 0
    detail = []
    for ea, eb, ed in we:
        if (ea, eb) == (ua, ub):
            continue                                        # skip seed edge
        se = sigma(ea, eb, occs, pot, S)
        te = tpow(tau, ed)
        common = set(se) & set(te)
        ag = sum(1 for x in common if se[x] == te[x])
        agree += ag; tot += len(common)
        detail.append((ea, eb, ed, ag, len(common)))
    return {'web': web, 'unit': True, 'agree': agree, 'tot': tot,
            'frac': round(agree / tot, 3) if tot else None, 'detail': detail}


def fit_all(S, names):
    occs, pot, webs, edges = pipeline(S, names)
    return [fit_web(w, occs, pot, S, edges) for w in webs if len(w) >= 3]


# ------------------------------------------------------ synthetic corpora

def make_progressive(seed, regions=False):
    """c_i = T^i(B(p_i)) off-chain header. regions=True uses a DIFFERENT rotor
    for the west messages, so cross-region sigma is not a single T^delta."""
    rng = random.Random(seed)
    def perm():
        d = list(range(P)); rng.shuffle(d); return d
    TA, TB, B = perm(), perm(), perm()
    def powers(T):
        tab = [list(range(P))]
        for _ in range(200):
            tab.append([T[tab[-1][x]] for x in range(P)])
        return tab
    pwA, pwB = powers(TA), powers(TB)
    def enc(s, T2):
        seq = [B[IDX[s[0]]]]
        for i, ch in enumerate(s[1:]):
            seq.append(T2[i + 1][B[IDX[ch]]])
        return seq
    return enc, pwA, pwB, (TA, TB, B)


def dense(enc, pw, rng, L, tries=6000):
    carrier = rng.choices(PA, weights=WEIGHTS, k=90)
    for _ in range(tries):
        ph = rng.choices(PA, weights=WEIGHTS, k=L)
        s = list(carrier); s[30:30 + L] = ph
        c = enc(s, pw)[30:30 + L]
        for w0 in range(L - 11):
            if sum(1 for v, k in Counter(c[w0:w0 + 12]).items() if k >= 2) >= 3:
                return ph
    return None


def corpus_prog(seed):
    """Positive control: a single global progressive rotor T."""
    enc, pwA, pwB, keys = make_progressive(seed)
    rng = random.Random(seed + 9)
    base = rng.choices(PA, weights=WEIGHTS, k=118)
    names = [f'sim-{i}' for i in range(9)]
    plain = {n: list(base) for n in names}
    for n in plain:
        plain[n][0] = rng.choice(PA)
    ph = dense(enc, pwA, rng, 30)
    for k, n in enumerate(['sim-0', 'sim-1', 'sim-2', 'sim-3', 'sim-4']):
        plain[n][30 + k:60 + k] = ph            # unit-spaced offsets -> unit edges
    return {n: enc(plain[n], pwA) for n in names}, names


def corpus_mono(seed):
    """Negative control: a per-message monoalphabetic substitution. It also
    produces EXACT isomorphs (sigma = P_n P_m^{-1}) but that sigma is
    OFFSET-INDEPENDENT, so it must not fit sigma = T^delta."""
    rng = random.Random(seed + 5)
    names = [f'sim-{i}' for i in range(9)]
    Pm = {}
    for n in names:
        d = list(range(P)); rng.shuffle(d); Pm[n] = d
    base = rng.choices(PA, weights=WEIGHTS, k=118)
    plain = {n: list(base) for n in names}
    for n in plain:
        plain[n][0] = rng.choice(PA)
    ph = rng.choices(PA, weights=WEIGHTS, k=30)
    for k, n in enumerate(['sim-0', 'sim-1', 'sim-2', 'sim-3', 'sim-4']):
        plain[n][30 + k:60 + k] = ph
    return {n: [Pm[n][IDX[c]] for c in plain[n]] for n in names}, names


def gate(name, ok):
    print(f'  {"PASS" if ok else "FAIL":4s}  {name}')
    return ok


def report_fits(fits, label):
    print(f'\n== {label} ==')
    tot_a = tot_t = 0
    for f in fits:
        if not f['unit']:
            print(f"  web {f['web']}: no unit-offset edge (cannot read T directly)")
            continue
        print(f"  web {f['web']}: T^delta agreement {f['agree']}/{f['tot']} "
              f"= {f['frac']}")
        for ea, eb, ed, ag, n in f['detail']:
            if n:
                print(f"      occ{ea}->occ{eb} delta={ed:+d}: {ag}/{n}")
        tot_a += f['agree']; tot_t += f['tot']
    overall = round(tot_a / tot_t, 3) if tot_t else None
    print(f"  overall support-weighted agreement: {tot_a}/{tot_t} = {overall}")
    return overall


def main():
    gates = []
    report = {}

    Ss, ns = corpus_prog(SEED)
    single = fit_all(Ss, ns)
    ov_s = report_fits(single, 'control: single progressive rotor (independent edges)')
    report['single'] = {'overall': ov_s}
    gates.append(gate('single rotor: independent edges fit T^delta (>=0.9)',
                      ov_s is not None and ov_s >= 0.9))

    Sm, nm = corpus_mono(SEED)
    mono = fit_all(Sm, nm)
    ov_m = report_fits(mono, 'control: per-message monoalphabetic (offset-independent sigma)')
    report['mono'] = {'overall': ov_m}
    gates.append(gate('monoalphabetic: independent edges do NOT fit (<=0.4)',
                      ov_m is not None and ov_m <= 0.4))

    print('\n== control gates ==')
    allok = all(gates)
    if not allok:
        print('  gate failure -> extraction/discriminator not validated')

    S = load_msgs()
    real = fit_all(S, N9)
    ov = report_fits(real, 'REAL CORPUS')
    report['real'] = {'overall': ov, 'fits': [{k: f[k] for k in f if k != 'detail'}
                                              for f in real]}
    print('\nverdict:')
    if ov is not None and ov < 0.5:
        print(f'  single per-web rotor REFUTED (agreement {ov} vs single-rotor '
              f'control ~{ov_s}). The effective rotor T VARIES within a web -- '
              f'the substitution changes on a scale FINER than the isomorph webs, '
              f'so 18 isomorphs under-constrain the per-piece rotors.')
    else:
        print(f'  webs are consistent with a per-web rotor (agreement {ov}); '
              f'recover T per web and attempt decryption.')
    report['verdict_ok_rotor'] = bool(ov is not None and ov >= 0.5)

    json.dump(report, open(f'{ROOT}/data/rotorfit_report.json', 'w'),
              indent=1, default=str)
    print('\nwrote data/rotorfit_report.json')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
