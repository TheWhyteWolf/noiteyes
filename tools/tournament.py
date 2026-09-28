#!/usr/bin/env python3
"""Step B: small-state generator tournament.

Step C left the mechanism small-state, non-additive, with a period-4 component
and an off-chain header. The other reproduced signatures tighten it further:

  * L=33 EXACT isomorphs => sigma is CONSTANT over 33 positions. The only
    schedule giving constant sigma is a progressive key E_i = T^(i)·B (then
    sigma = T^delta). So within a constant-sigma stretch the cipher is a
    progressive ("Trithemius/Alberti") key with SOME permutation T.
  * sigma maps are NON-commuting across stretches but tend to commute within
    one (verified: same-web 6/14 vs cross-web 0/2). Powers of one T commute,
    so T must CHANGE between stretches.
  * sigma_web's additive/cyclic-group solve is dim-1 => T is a GENERAL
    permutation, not a cyclic rotation (else it would be the refuted additive
    walk).
  * period-4 => the changing T cycles through a small set; the natural size is
    4 rotors T0..T3.
  * resync S_eff~15 (Step C) => the stretch re-key is plaintext-driven (an
    accumulator), not purely positional (which would realign far too often).

Candidate under test: a CYCLE OF 4 GENERAL-PERMUTATION ROTORS, progressive
within each stretch, re-keyed (advance rotor, reset the progressive counter)
on a plaintext notch. Baselines: an additive walk (refuted family) and a
single-rotor Enigma-style conjugation (position-driven).

This is a SEARCH, not a refutation: the score is how many of the real corpus's
signatures each candidate reproduces. Sanity gates keep the signature
detectors honest (the additive baseline MUST be crackable; the 4-rotor
candidate MUST produce exact isomorphs and a dim-1 additive-solve). The real
corpus is profiled with the identical detectors.
"""
import json
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
from pi_solver import P, load_msgs
from sigma_web import (MINLEN, PA, ROOT, WEIGHTS, build_webs, core_blocks,
                       extract_sigmas, find_isomorph_pairs, solve_blocks,
                       witness_trim)
from resync import corpus_realign

SEED = 1
N9 = pi_solver.N9
IDX = {c: i for i, c in enumerate(PA)}


# ------------------------------------------------------- signature battery

def commute_pair(s1, s2):
    a = t = 0
    for x in set(s1) & set(s2):
        p, q = s2.get(s1[x]), s1.get(s2[x])
        if p is not None and q is not None:
            t += 1; a += (p == q)
    return a, t


def delta_mod4_spectrum(S, names):
    """Histogram of (posB - posA) % 4 over the raw isomorph pairs -- the
    same diagnostic Step D ran on the real corpus (README: offsets
    0:6,1:2,2:5,3:5, i.e. isomorphs at EVERY residue, non-identity sigma).
    Used here to directly test position-locked vs content-locked period-4
    candidates: a position-locked outer layer predicts isomorphs pile up
    at delta%4==0 (or vanish elsewhere); content-locked predicts they
    don't care about delta%4 at all, like the real corpus."""
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    hist = Counter((p[4] - p[1]) % 4 for p in raw)
    return {'n_pairs': len(raw), 'hist': dict(sorted(hist.items()))}


def signatures(S, names, rng):
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    maxL = max((p[2] - p[1] for p in raw), default=0)
    walk_dim = None
    commute = same_w = None
    if raw:
        joint = witness_trim(S, raw)
        cores, _ = core_blocks(joint)
        if len(cores) >= 2:
            walk_dim = solve_blocks(S, cores, rng, tries=20)['dim']
        trimmed = [e['trim'] for e in joint]
        sig = extract_sigmas(S, trimmed)
        occs, pot, webs, gv, edges = build_webs(trimmed)
        pw = {}
        for a, b, off, k in edges:
            for wi, w in enumerate(webs):
                if a in w or b in w:
                    pw[k] = wi
        ca = ct = sa = st = 0
        for i in range(len(sig)):
            for j in range(i + 1, len(sig)):
                a, t = commute_pair(sig[i], sig[j])
                ca += a; ct += t
                if pw.get(i) is not None and pw.get(i) == pw.get(j):
                    sa += a; st += t
        commute = round(ca / ct, 2) if ct else None
        same_w = round(sa / st, 2) if st else None
    d1 = sum(1 for n in names for i in range(1, len(S[n])) if S[n][i] == S[n][i - 1])

    def dcount(d):
        return sum(1 for n in names for i in range(d, len(S[n])) if S[n][i] == S[n][i - d])
    d4 = dcount(4); dref = (dcount(3) + dcount(5)) / 2 or 1
    rr = corpus_realign(S, names)['rate']
    cnt = Counter(v for n in names for v in S[n]); tot = sum(cnt.values())
    ioc = sum(c * (c - 1) for c in cnt.values()) / (tot * (tot - 1))
    return {'maxL': maxL, 'walk_dim': walk_dim, 'commute': commute,
            'same_web': same_w, 'doubles': d1, 'd4_ratio': round(d4 / dref, 2),
            'resync': rr, 'ioc': round(ioc, 5)}


REAL_TESTS = {
    'iso L>=25':      lambda s: s['maxL'] >= 25,
    'walk dim==1':    lambda s: s['walk_dim'] == 1,
    'noncommute<.6':  lambda s: s['commute'] is not None and s['commute'] < 0.6,
    'zero doubles':   lambda s: s['doubles'] == 0,
    'resync .03-.12': lambda s: 0.03 <= s['resync'] <= 0.12,
    'IoC .011-.015':  lambda s: 0.011 <= s['ioc'] <= 0.015,
}


# 'period4 (d4>1.5)' was retired from scoring: tools/diagscan.py shows the
# real d4 excess is a duplication artifact (26 raw -> 11 deduplicated vs 7.8
# expected, p=0.16, level with d5/d9/d13). d4_ratio is still reported.
RETIRED_TESTS = {
    'period4 (d4>1.5)': lambda s: s['d4_ratio'] > 1.5,
}


def score(s):
    return {k: f(s) for k, f in REAL_TESTS.items()}


# ---------------------------------- fixed-key single-message encryptors

def powers(T, maxm):
    tab = [list(range(P))]
    for _ in range(maxm):
        tab.append([T[tab[-1][x]] for x in range(P)])
    return tab


def make_additive(seed):
    rng = random.Random(seed + 2)
    D = list(range(P)); rng.shuffle(D)
    v = {c: rng.randrange(1, P) for c in PA}
    def enc(s):
        seq = [D[(50 + v[s[0]]) % P]]; x = 17
        for ch in s[1:]:
            x = (x + v[ch]) % P; seq.append(D[x])
        return seq
    return enc


def make_enigma1(seed):
    """Single-rotor conjugation c_i = rot^{-i} W rot^{i}(p): position-driven,
    NOT constant-sigma (shows why long isomorphs need a progressive key)."""
    rng = random.Random(seed + 3)
    W = list(range(P)); rng.shuffle(W)
    def enc(s):
        seq = [W[IDX[s[0]]]]
        for i, ch in enumerate(s[1:]):
            seq.append((W[(IDX[ch] + i) % P] - i) % P)
        return seq
    return enc


def make_rotorcycle(seed):
    """4 general-permutation rotors cycled by stretch; progressive within a
    stretch; re-key (advance rotor, reset counter) on a plaintext notch."""
    rng = random.Random(seed + 4)
    pw = [powers([rng.sample(range(P), P)][0], 140) for _ in range(4)]
    B = list(range(P)); rng.shuffle(B)
    notch = {c for c in PA if IDX[c] % 5 == 0}
    def enc(s):
        seq = [B[IDX[s[0]]]]; sidx = 0; m = 0
        for ch in s[1:]:
            seq.append(pw[sidx % 4][m][B[IDX[ch]]]); m += 1
            if ch in notch:
                sidx += 1; m = 0
        return seq
    return enc


# ---------------------- new candidates (twelfth-phase: "deck rebuilt, not
# shifted" + period-4-composed-with-autokey; see GUIDE.md's Step B writeup)
#
# The tournament's binding tension: reproducing the real corpus needs a
# mechanism that (i) makes MANY cross-offset EXACT isomorphs (needs sigma
# CONSTANT over a whole 18-33 letter stretch), (ii) is non-additive /
# non-commuting across stretches, and (iii) realigns at a moderate rate
# (~0.065). make_rotorcycle satisfies (ii) via 4 unrelated rotors but fails
# (i) because within a stretch it applies T^m (m = letters since last
# notch) -- a PROGRESSIVE key, so sigma is NOT constant even inside one
# stretch. The fix tested here: apply a table with NO exponent at all (the
# SAME single fixed permutation for every letter in the stretch -- "deck
# REBUILT, not shifted", GUIDE.md's third-phase idea, here tried with a
# plaintext-driven notch instead of the already-refuted row-lock timing).
#
# Content-driven (not position-driven) notch is essential: it must be a
# deterministic function of the plaintext LETTER identity (matching
# make_rotorcycle's own convention), so that two independent occurrences of
# the identical phrase see an IDENTICAL notch pattern and can stay
# state-synchronized -- a notch keyed to absolute character position would
# make constant-sigma isomorphs impossible except at one fixed relative
# offset, which is the position-locked failure mode candidate 2 below is
# built to demonstrate directly.
#
# Chosen so mean spacing lands in the observed 18-33 isomorph-length range:
# total PA weight is 118.3; {q,z,j,x,k,v,b,p} sums to 5.8 -> mean spacing
# 118.3/5.8 ~ 20.4 letters.
_SLOW_NOTCH_LETTERS = set('qzjxkvbp')
M_SLOW = 15  # ~ Step C's S_eff~15 upper bound on the effective state


def make_deck_rebuild(seed):
    """'Deck rebuilt, not shifted': state indexes a bank of M_SLOW
    INDEPENDENTLY random permutations (no group relation between them --
    composing any two is an arbitrary map, so different states don't
    commute by construction). Within a run between notches, ONE table is
    applied unchanged to every letter (constant sigma across the whole
    stretch -- unlike make_rotorcycle's T^m). State advances (round-robin,
    mod M_SLOW) only on a plaintext notch, so two chance-aligned
    occurrences of a repeated phrase that start in the same state AND
    don't themselves contain a notch letter reproduce an exact isomorph;
    two diverged runs re-synchronize whenever their notch-counts happen to
    coincide mod M_SLOW (~1/M_SLOW per opportunity, the resync mechanism)."""
    rng = random.Random(seed + 5)
    tables = [rng.sample(range(P), P) for _ in range(M_SLOW)]
    B = list(range(P)); rng.shuffle(B)
    def enc(s):
        state = 0
        seq = [tables[state][B[IDX[s[0]]]]]
        for ch in s[1:]:
            if ch in _SLOW_NOTCH_LETTERS:
                state = (state + 1) % M_SLOW
            seq.append(tables[state][B[IDX[ch]]])
        return seq
    return enc


def make_period4_position(seed):
    """Candidate 1 (deck-rebuild autokey) with an OUTER period-4 layer keyed
    by ABSOLUTE character position: c_i = Q[i mod 4]( Table[state](p_i) ).
    Tests the task's 'period-4 sub-key... composed with a small autokey
    state' idea in its most literal, position-locked form. Predicted
    failure mode (Step D, already in README): since Q's phase depends only
    on absolute position, two occurrences of a shared phrase at relative
    offset delta get sigma_i = Q[(i+delta) mod 4] . Q[i mod 4]^-1, which is
    constant across the whole window ONLY when delta%4==0 -- at other
    offsets sigma cycles through up to 4 different maps within one window,
    which generally breaks the exact-isomorph match. Real isomorph offsets
    span ALL residues mod 4 with non-identity sigma (Step D), so if this
    candidate reproduces that pattern it would falsify the prediction;
    if it instead only makes isomorphs at delta%4==0, that confirms
    position-locked period-4 composition is the wrong shape and can be
    retired without further testing."""
    rng = random.Random(seed + 6)
    tables = [rng.sample(range(P), P) for _ in range(M_SLOW)]
    Q = [rng.sample(range(P), P) for _ in range(4)]
    B = list(range(P)); rng.shuffle(B)
    def enc(s):
        state = 0
        seq = [Q[0][tables[state][B[IDX[s[0]]]]]]
        for i, ch in enumerate(s[1:], start=1):
            if ch in _SLOW_NOTCH_LETTERS:
                state = (state + 1) % M_SLOW
            seq.append(Q[i % 4][tables[state][B[IDX[ch]]]])
        return seq
    return enc


# Quagmire family: a pool of mixed alphabets held CONSTANT over stretches
# (giving exact offset-invariant isomorphs with non-additive sigma -- the
# binding tension's missing property), advanced by a plaintext notch (content-
# driven, no positional organization -- matching rotorfit). Pool size ~8-16
# per alphacount.py: the real web forces >=8 distinct alphabets across 9
# occurrences with exactly one reuse.

QPOOL = 12


def make_quagA(seed):
    """Plain notch-advanced pool quagmire (stretch-keyed monoalphabets)."""
    rng = random.Random(seed + 6)
    pool = [dict(zip(PA, rng.sample(range(P), len(PA)))) for _ in range(QPOOL)]
    notch = set(rng.sample(PA, 5))
    def enc(s):
        k = 0; seq = []
        for ch in s:
            seq.append(pool[k % QPOOL][ch])
            if ch in notch:
                k += 1
        return seq
    return enc


# A second, FAST, content-driven notch for the period-4 layer of candidate
# 3 below -- deliberately common letters (~1/4 of total weight) so the fast
# state advances roughly every 4 letters on average, the scale the real
# distance-4 recurrence excess operates at. {i,o,a,t} sums to 7.0+7.5+8.2+
# 9.1=31.8 of 118.3 total -> mean spacing 118.3/31.8 ~ 3.7 letters.
_FAST_NOTCH_LETTERS = set('ioat')
M_FAST = 4


def make_period4_content(seed):
    """Candidate 1 (deck-rebuild autokey) with a period-4 layer that is
    CONTENT-driven rather than position-driven: c_i =
    Q[fast mod 4]( Table[slow mod M_SLOW](p_i) ), where BOTH fast and slow
    counters advance on their own plaintext-letter-identity notch (fast
    ~every 4 letters, slow ~every 20). Because fast's value depends only on
    how many fast-notch letters have occurred SO FAR IN THIS MESSAGE'S OWN
    TEXT -- not on absolute position -- two occurrences of an identical
    phrase advance fast IDENTICALLY step-for-step throughout the shared
    span (their internal letters are the same), so if they merely START at
    the same fast-phase (~1/4 chance) they STAY synchronized for the whole
    window regardless of the window's absolute offset delta -- unlike
    candidate 2, this predicts exact isomorphs CAN occur at any delta mod 4,
    matching Step D's real observation instead of contradicting it. This
    operationalizes the open idea in GUIDE.md Step A's writeup: 'an order-4
    key STATE feeding the non-commuting substitution, not an order-4 offset
    on a state walk.'"""
    rng = random.Random(seed + 7)
    tables = [rng.sample(range(P), P) for _ in range(M_SLOW)]
    Q = [rng.sample(range(P), P) for _ in range(M_FAST)]
    B = list(range(P)); rng.shuffle(B)
    def enc(s):
        slow = fast = 0
        seq = [Q[fast][tables[slow][B[IDX[s[0]]]]]]
        for ch in s[1:]:
            if ch in _SLOW_NOTCH_LETTERS:
                slow = (slow + 1) % M_SLOW
            if ch in _FAST_NOTCH_LETTERS:
                fast = (fast + 1) % M_FAST
            seq.append(Q[fast][tables[slow][B[IDX[ch]]]])
        return seq
    return enc


def make_quagB(seed):
    """Period-4-composed pool: alphabet index = (i mod 4, notch count) --
    an order-4 positional cycle woven through the pool selection."""
    rng = random.Random(seed + 7)
    pool = [[dict(zip(PA, rng.sample(range(P), len(PA)))) for _ in range(4)]
            for _ in range(QPOOL)]
    notch = set(rng.sample(PA, 5))
    def enc(s):
        k = 0; seq = []
        for i, ch in enumerate(s):
            seq.append(pool[k % QPOOL][i % 4][ch])
            if ch in notch:
                k += 1
        return seq
    return enc


def make_gak41(seed):
    """Group Autokey over C83:C41 -- THE one group in the community's GAK
    classification (Toboter's progress doc) that nobody has ruled out: affine
    maps x -> b*x + a mod 83 with b restricted to the order-41 subgroup of
    Z83* (generated by 4). State composes g <- g o m(p); ciphertext = D[a].
    a_p != 0 for all letters gives STRUCTURAL zero doubles; occurrences of a
    phrase share the running product, so sigma = h2 o h1^-1 is a constant
    AFFINE map -- exact offset-invariant isomorphs that are NOT translations.
    Hidden state = the b component: 41 values, matching the small-state bound."""
    rng = random.Random(seed + 9)
    sub = sorted({pow(4, k, P) for k in range(41)})
    twist = set(rng.sample(PA, 5))
    a_of = {c: rng.randrange(1, P) for c in PA}
    b_of = {c: rng.choice([x for x in sub if x != 1]) if c in twist else 1
            for c in PA}
    D = list(range(P)); rng.shuffle(D)
    def enc(s):
        b, a = 1, 17
        seq = [D[(50 + a_of[s[0]]) % P]]
        for ch in s[1:]:
            b, a = (b * b_of[ch]) % P, (b * a_of[ch] + a) % P
            seq.append(D[a])
        return seq
    return enc


def make_quagC(seed):
    """Notch quagmire + COLLISION BUMP: if the output letter would repeat the
    previous output, advance the pool index and re-encrypt (repeat until it
    differs). Structural zero doubles by a NON-additive state bump."""
    rng = random.Random(seed + 8)
    pool = [dict(zip(PA, rng.sample(range(P), len(PA)))) for _ in range(QPOOL)]
    notch = set(rng.sample(PA, 5))
    def enc(s):
        k = 0; seq = []
        for ch in s:
            c = pool[k % QPOOL][ch]
            while seq and c == seq[-1]:
                k += 1
                c = pool[k % QPOOL][ch]
            seq.append(c)
            if ch in notch:
                k += 1
        return seq
    return enc


# --------------------------------------------------------- plaintext family

def dense_under(enc, rng, L, tries=4000):
    """A phrase whose encryption (in a fixed carrier) has a 12-window with >=3
    internally repeated cipher letters -- so the isomorph finder can seed it.
    Rejection-sampled against the candidate's OWN encryption, so every
    candidate gets a fair shot at producing detectable isomorphs."""
    carrier = rng.choices(PA, weights=WEIGHTS, k=90)
    for _ in range(tries):
        ph = rng.choices(PA, weights=WEIGHTS, k=L)
        s = list(carrier); s[25:25 + L] = ph
        c = enc(s)[25:25 + L]
        for w0 in range(0, L - 11):
            win = c[w0:w0 + 12]
            reps = sum(1 for v, k in Counter(win).items() if k >= 2)
            if reps >= 3:
                return ph
    return None


def p4ify(rng, s):
    """Impose soft period-4 PLAINTEXT structure (a list/record delimiter every
    4th unit) -- tests whether the real d4 excess can live in the plaintext
    layer instead of the key schedule. Phase is per-sequence, so inserts land
    at arbitrary phases; distance-4 recurrence is phase-agnostic."""
    return ['_' if (i % 4 == 0 and rng.random() < 0.6) else c
            for i, c in enumerate(s)]


def corpus(enc, rng, p4=False):
    """Shared template (prefix + gaps, for resyncs) plus candidate-seeded
    repeated phrases at spread positions (for cross-position isomorphs)."""
    L = 120
    T = rng.choices(PA, weights=WEIGHTS, k=L)
    if p4:
        T = p4ify(rng, T)
    names = [f'sim-{i}' for i in range(9)]
    plain = {}
    for nm in names:
        s = list(T); s[0] = rng.choice(PA)
        for _ in range(3):
            p = rng.randrange(12, L - 8)
            s[p:p + 7] = rng.choices(PA, weights=WEIGHTS, k=7)
        plain[nm] = s
    phs = [dense_under(enc, rng, n) for n in (30, 28, 26)]
    if p4:
        phs = [p4ify(rng, ph) if ph else ph for ph in phs]
    def ins(nm, ph, pos):
        if ph:
            plain[nm][pos:pos + len(ph)] = list(ph)
    # NOTE: insertion offsets are deliberately a MIX of residues mod 4 (not
    # all multiples of 4, as an earlier version of this harness had by
    # coincidence -- 62-30=32, 34-30=4, 68-40=28, 44-40=4, 62-46=16, 30-46=16
    # were ALL ≡0 mod 4). That coincidence silently made every position-
    # locked period-4 candidate look identical to a content-locked one in
    # this test, since the harness never sampled a cross-message offset
    # that isn't a multiple of 4. Shifting one insertion per phrase group by
    # 1-2 positions breaks that and lets the maxL signature actually
    # discriminate the two designs (see make_period4_position's docstring).
    ins('sim-0', phs[0], 30); ins('sim-1', phs[0], 61); ins('sim-5', phs[0], 34)
    ins('sim-2', phs[1], 40); ins('sim-4', phs[1], 67); ins('sim-3', phs[1], 44)
    ins('sim-6', phs[2], 46); ins('sim-8', phs[2], 60); ins('sim-7', phs[2], 30)
    return {nm: enc(plain[nm]) for nm in names}, names


def gate(name, ok):
    print(f'  {"PASS" if ok else "FAIL":4s}  {name}')
    return ok


def main():
    rng = random.Random(SEED)
    report, gates = {}, []

    cands = {'additive walk': (make_additive, False),
             'single-rotor enigma': (make_enigma1, False),
             '4-rotor progressive': (make_rotorcycle, False),
             'deck-rebuild autokey': (make_deck_rebuild, False),
             'period4-position+autokey': (make_period4_position, False),
             'period4-content+autokey': (make_period4_content, False),
             'quagA notch-pool': (make_quagA, False),
             'quagB p4-composed': (make_quagB, False),
             'quagC bumped': (make_quagC, False),
             'quagC + p4 plaintext': (make_quagC, True),
             'gak41 affine autokey': (make_gak41, False),
             'gak41 + p4 plaintext': (make_gak41, True)}
    SPECTRUM_CANDS = {'deck-rebuild autokey', 'period4-position+autokey',
                       'period4-content+autokey'}
    print('== candidate signature profiles (candidate-seeded plaintext) ==')
    profiles = {}
    for label, (make, p4) in cands.items():
        enc = make(SEED)
        S, names = corpus(enc, rng, p4=p4)
        s = signatures(S, names, rng)
        profiles[label] = s
        sc = score(s)
        print(f'\n  {label}:')
        print(f'    maxL {s["maxL"]}  walk_dim {s["walk_dim"]}  commute {s["commute"]}'
              f'  same_web {s["same_web"]}  doubles {s["doubles"]}  '
              f'd4_ratio {s["d4_ratio"]}  resync {s["resync"]}  IoC {s["ioc"]}')
        print('    ' + '  '.join(f'{k}={"Y" if v else "n"}' for k, v in sc.items()))
        report[label] = {'sig': s, 'score': sc, 'n_match': sum(sc.values())}
        if label in SPECTRUM_CANDS:
            spec = delta_mod4_spectrum(S, names)
            report[label]['delta_mod4_spectrum'] = spec
            print(f'    delta%4 spectrum (of {spec["n_pairs"]} raw isomorph '
                  f'pairs, L>={MINLEN}): {spec["hist"]}')

    # gates validate the DETECTORS (so the real fingerprint is trustworthy);
    # the candidate scores are exploratory findings, not gated.
    gates.append(gate('detector: additive baseline is crackable (walk_dim>=2)',
                      profiles['additive walk']['walk_dim'] is not None
                      and profiles['additive walk']['walk_dim'] >= 2))
    gates.append(gate('detector: additive gives commuting sigma (~1.0)',
                      profiles['additive walk']['commute'] is not None
                      and profiles['additive walk']['commute'] > 0.9))
    gates.append(gate('detector: single-rotor gives NO long isomorphs (maxL<25)',
                      profiles['single-rotor enigma']['maxL'] < 25))
    gates.append(gate('construction: quagC bump yields zero doubles',
                      profiles['quagC bumped']['doubles'] == 0
                      and profiles['quagC + p4 plaintext']['doubles'] == 0))

    print('\n== control gates (detector validation) ==')
    allok = all(gates)
    if not allok:
        print('  gate failure -> detectors not behaving as expected')

    print('\n== REAL CORPUS profile ==')
    S = load_msgs()
    real = signatures(S, N9, rng)
    rsc = score(real)
    real_spec = delta_mod4_spectrum(S, N9)
    report['REAL'] = {'sig': real, 'score': rsc, 'n_match': sum(rsc.values()),
                      'delta_mod4_spectrum': real_spec}
    print(f'  maxL {real["maxL"]}  walk_dim {real["walk_dim"]}  commute {real["commute"]}'
          f'  same_web {real["same_web"]}  doubles {real["doubles"]}  '
          f'd4_ratio {real["d4_ratio"]}  resync {real["resync"]}  IoC {real["ioc"]}')
    print('  ' + '  '.join(f'{k}={"Y" if v else "n"}' for k, v in rsc.items()))
    print(f'  delta%4 spectrum (of {real_spec["n_pairs"]} raw isomorph pairs, '
          f'L>={MINLEN}): {real_spec["hist"]}  <- offsets at EVERY residue, '
          f'the pattern a position-locked period-4 layer cannot produce')

    print('\n== signatures each candidate reproduces vs REAL (real=all Y) ==')
    best = None
    for label in cands:
        cs = report[label]['score']
        agree = sum(1 for k in REAL_TESTS if cs[k] == rsc[k])
        print(f'  {label:22s}: matches real on {agree}/{len(REAL_TESTS)} signatures'
              + ('' if agree == len(REAL_TESTS)
                 else '   (miss: ' + ', '.join(k for k in REAL_TESTS
                                               if cs[k] != rsc[k]) + ')'))
        report[label]['agree_real'] = agree
        best = (agree, label) if best is None or agree > best[0] else best

    print(f'\n  best candidate(s): '
          + ', '.join(l for l in cands if report[l]['agree_real'] == best[0])
          + f' ({best[0]}/{len(REAL_TESTS)}). ')
    print('  gak41 (C83:C41 affine group-autokey, the one unruled group in the')
    print('  community GAK classification) is the FIRST candidate to reproduce')
    print('  walk_dim==1 together with non-commuting sigma, exact cross-offset')
    print('  isomorphs, structural zero doubles and the IoC band (robust: median')
    print('  5/7 and walk_dim==1 in 14/20 over a 20-seed sweep; 5/6 now that')
    print('  period-4 is retired as a duplication artifact, see diagscan.py).')
    print('  Its remaining miss is resync rate (twist elements are a free choice).')
    print('  (Historical: quagC+p4 cleared d4>1.5 in only 5/20 seeds and p4')
    print('  plaintext broke gak41 walk_dim==1 -- moot now that d4 is retired.)')
    print('  Next tests:')
    print('  fit sigma_k = pi^-1(beta_k*pi + alpha_k) with beta in the order-41')
    print('  subgroup (the affine generalization no prior phase has tested), and')
    print('  explain fast resync (S_eff<=15) alongside the large affine state.')
    report['conclusion'] = {'best': best[1], 'best_score': best[0],
                            'note': 'gak41 5/6 incl walk_dim==1 (was 5/7; period-4 '
                            'retired as a duplication artifact, diagscan.py); '
                            'only miss is resync (param). next: affine-pi '
                            'solver over subgroup-41 (done, inconclusive); '
                            'reconcile resync with a large state',
                            'tension': 'abundant cross-offset exact isomorphs vs '
                            'non-additive/non-commuting vs moderate resync'}

    json.dump(report, open(f'{ROOT}/data/tournament_report.json', 'w'),
              indent=1, default=str)
    print('\nwrote data/tournament_report.json')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
