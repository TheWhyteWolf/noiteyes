#!/usr/bin/env python3
"""Tournament round 5: a LOSSY state update for gak41 that resyncs.

Round 4 showed gak41 + short additive steps (1..27) fixes near-repeat
depletion, but resync stays ~0: the update g <- g o m(p) is a group action,
hence invertible, so two diverged messages re-agree only if their differing
plaintext has the same group product (a coincidence). A plaintext reset
(rare trigger letters) reaches the real mean resync but is all-or-nothing:
whether a trigger follows each divergence gap is fixed by the SHARED template
plaintext, so the per-seed rate swings between 0 and far too high (resync
band hit in only 15% of 48 seeds), and resets cost walk dim / non-commutation.

Round 5 makes the TWIST state-dependent, so the update map becomes
non-injective and two diverged states can MERGE by chance (a per-pair event,
not a per-template one):
  ptwist      -- round-4 short K=27 baseline: twist b_of[p] on 5 plaintext letters
  bset        -- plaintext twist letters OVERWRITE b (b <- b_of[p]) instead of
                 multiplying: b forgets its history, so two diverged states
                 merge whenever their a's agree at a twist letter
  stw-mul m   -- twist driven by the STATE: b <- b * t(a) when the current a is
                 in a set of m values (plaintext twists off)
  stw-mul m+p -- the same, plus the plaintext twists
  stw-set m+p -- plaintext twists multiply; state twists overwrite b <- t(a)
Each m in (2, 4, 8, 16). Plus:
  ptwist tw=t -- the lossless baseline with only t of its 5 twist letters
  coarse C    -- coarse-step autokey (make_coarse): update by letter CLASS
Two side tests per config: unaligned repeated ciphertext runs >= 4 (real: 0;
a reset predicts them) and resync under a SHORT-EDIT plaintext model (1-4
letter edits, like the real 1-4 letter mismatch islands) instead of the
harness's 7-letter blocks. Same corpus harness and 7 tests as rounds 3-4,
200 seeds by default; jobs run in parallel (meant for enter).

Result (200 seeds, data/tournament5_report.json):
  * state-driven twists destroy the isomorphs (iso L>=25 falls to 0.02-0.97
    as m grows, walk dim <= 0.16): two occurrences of a phrase start in
    different states, so they twist at different letters and sigma is no
    longer constant. Loss compatible with exact isomorphs must be
    PLAINTEXT-driven.
  * bset (plaintext overwrite of b) kills non-commutation (0.04).
  * resets stay bimodal and predict unaligned repeats the real corpus lacks
    (only 55% of seeds have none; resync band AND none: 9%).
  * the resync test is plaintext-model dependent: under short edits the
    LOSSLESS baseline with 1-3 twist letters lands in the band in 12-39% of
    seeds (5 twists: 5%), against 5% for resets; only bset (54%, but no
    non-commutation) does better. Round 4's "resync needs a lossy update"
    does not hold beyond the 7-letter-block harness.
  * coarse-step autokeys keep isomorphs but fail zero doubles (<= 6%).

Writes data/tournament5_report.json. Run from repo root:
    python3 tools/tournament5.py [n_seeds] [n_procs]
"""
import json
import os
import random
import statistics
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tournament import P, PA, WEIGHTS, load_msgs, N9, signatures  # noqa: E402
from tournament3 import depletion, tests, sigma_stat, summarize  # noqa: E402
from tournament import corpus  # noqa: E402
from resync import corpus_realign  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FREQ = {c: w / sum(WEIGHTS) for c, w in zip(PA, WEIGHTS)}
RARE4 = sorted(PA, key=lambda c: FREQ[c])[:4]


def make_lossy(seed, mode='ptwist', m=0, K=27, trig=(), ntw=5):
    """gak41 key material (tournament4.make_gak41x with short steps K) plus a
    state-dependent twist. mode: ptwist | bset | stw-mul | stw-mul+p | stw-set+p."""
    rng = random.Random(seed + 9)
    sub = sorted({pow(4, k, P) for k in range(41)})
    twist = set(rng.sample(PA, 5)[:ntw])            # ntw < 5: fewer twist letters
    a_of = {c: rng.randrange(1, P) for c in PA}   # keeps the rng stream aligned
    if K:
        a_of = dict(zip(PA, random.Random(seed * 7 + K).sample(range(1, K + 1), len(PA))))
    b_of = {c: rng.choice([x for x in sub if x != 1]) if c in twist else 1
            for c in PA}
    D = list(range(P)); rng.shuffle(D)
    srng = random.Random(seed * 101 + m)
    tset = set(srng.sample(range(P), m)) if m else set()
    t_of = {v: srng.choice([x for x in sub if x != 1]) for v in tset}
    ptw = mode in ('ptwist', 'stw-mul+p', 'stw-set+p', 'bset')
    trig = set(trig)

    def enc(s):
        b, a = 1, 17
        seq = [D[(50 + a_of[s[0]]) % P]]
        for ch in s[1:]:
            if ch in trig:
                b, a = 1, 17
            nb = b                            # as gak41: a uses the OLD b
            if a in tset:
                if mode.startswith('stw-mul'):
                    nb = (nb * t_of[a]) % P
                elif mode == 'stw-set+p':
                    nb = t_of[a]
            if ptw and b_of[ch] != 1:
                nb = b_of[ch] if mode == 'bset' else (nb * b_of[ch]) % P
            b, a = nb, (b * a_of[ch] + a) % P
            seq.append(D[a])
        return seq
    enc.D = D                                 # hidden order pi = D^-1 (stepfit gates)
    enc.twist = [(a_of[c], b_of[c]) for c in sorted(twist)]
    return enc


def make_coarse(seed, C=8, K=27, ntw=2):
    """Coarse-step autokey: the state g = (b, a) is APPLIED to a fixed letter
    code, c_t = D[b*e(p_t) + a], then updated by the letter's CLASS only,
    g <- g o m(kappa(p_t)), with kappa: 27 letters -> C classes (C=27 is
    injective). Letters in one class share a key step, so a plaintext
    difference within a class changes one ciphertext letter and nothing after
    it (a 1-letter island), and differing words with equal class products
    re-merge: loss that is PLAINTEXT-driven, as exact isomorphs require.
    ntw of the C classes twist b."""
    rng = random.Random(seed * 13 + C)
    sub = sorted({pow(4, k, P) for k in range(41)})
    cls = list(range(C)) + [rng.randrange(C) for _ in range(len(PA) - C)]
    rng.shuffle(cls)
    kappa = dict(zip(PA, cls))
    steps = rng.sample(range(1, K + 1), C)
    twc = set(rng.sample(range(C), ntw))
    bcl = {k: rng.choice([x for x in sub if x != 1]) if k in twc else 1 for k in range(C)}
    e = dict(zip(PA, rng.sample(range(P), len(PA))))
    D = list(range(P)); rng.shuffle(D)

    def enc(s):
        b, a = 1, 17
        seq = [D[(50 + e[s[0]]) % P]]
        for ch in s[1:]:
            seq.append(D[(b * e[ch] + a) % P])
            k = kappa[ch]
            b, a = (b * bcl[k]) % P, (b * steps[k] + a) % P
        return seq
    return enc


def offreps(S, names, L=4):
    """Maximal exact ciphertext runs >= L shared at DIFFERENT positions
    (unaligned). A reset makes the text after a trigger a fixed function of
    the plaintext, so repeated plaintext after a trigger repeats in the
    ciphertext at any offset. Real corpus: none."""
    out = []
    for i, a in enumerate(names):
        for b in names[i:]:
            A, B = S[a], S[b]
            for o in range(-len(A) + 1, len(B)):
                if o == 0 or (a == b and o < 0):
                    continue
                run = 0
                hi = min(len(A), len(B) - o)
                for x in range(max(0, -o), hi + 1):
                    if x < hi and A[x] == B[x + o]:
                        run += 1
                    else:
                        if run >= L:
                            out.append((a, x - run, b, x - run + o, run))
                        run = 0
    return out


def corpus_edits(enc, rng, glen=(1, 4), n_ed=4, L=120):
    """Alternative plaintext model: shared template with n_ed SHORT edits
    (1-4 letters) per message instead of the harness's 7-letter blocks."""
    T = rng.choices(PA, weights=WEIGHTS, k=L)
    names = [f'sim-{i}' for i in range(9)]
    out = {}
    for nm in names:
        s = list(T); s[0] = rng.choice(PA)
        for _ in range(n_ed):
            g = rng.randint(*glen); p = rng.randrange(12, L - g)
            s[p:p + g] = rng.choices(PA, weights=WEIGHTS, k=g)
        out[nm] = enc(s)
    return out, names


def build(kw, sd):
    return make_coarse(sd, **kw['coarse']) if 'coarse' in kw else make_lossy(sd, **kw)


def job(args):
    label, kw, sd = args
    enc = build(kw, sd)
    rng = random.Random(1000 + sd)
    S, names = corpus(enc, rng)
    sig = signatures(S, names, rng)
    dep = depletion(S, names)
    ss = sigma_stat(S, names)
    E, en = corpus_edits(enc, random.Random(500 + sd))
    return label, sd, {'sig': sig, 'dep': dep, 'n_offrep': len(offreps(S, names)),
                       'resync_edits': corpus_realign(E, en)['rate'],
                       'sig_groups': None if ss is None else
                       {k: ss[k] for k in ('n_sigmas', 'n_groups', 'n_links',
                                           'links_by_kind', 'mean_fixed_frac',
                                           'nonid_cross_web')}}


def configs():
    c = {'ptwist (short K=27)': {'mode': 'ptwist'},
         'ptwist + reset rare4': {'mode': 'ptwist', 'trig': RARE4},
         'bset': {'mode': 'bset'}}
    for t in (1, 2, 3):
        c[f'ptwist tw={t}'] = {'mode': 'ptwist', 'ntw': t}
    for m in (2, 4, 8, 16):
        c[f'stw-mul m={m}'] = {'mode': 'stw-mul', 'm': m}
        c[f'stw-mul m={m} + p'] = {'mode': 'stw-mul+p', 'm': m}
        c[f'stw-set m={m} + p'] = {'mode': 'stw-set+p', 'm': m}
    for C in (3, 5, 8, 12, 27):
        c[f'coarse C={C}'] = {'coarse': {'C': C}}
    return c


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    seeds = list(range(1, n + 1))
    S = load_msgs()
    real_sig = signatures(S, N9, random.Random(1))
    real_dep = depletion(S, N9)
    cf = configs()
    jobs = [(lab, kw, sd) for lab, kw in cf.items() for sd in seeds]
    rows = {lab: {} for lab in cf}
    with Pool(procs) as pool:
        for lab, sd, row in pool.imap_unordered(job, jobs, chunksize=4):
            rows[lab][sd] = row
    report = {'seeds': n, 'REAL': {'sig': real_sig, 'dep': real_dep,
                                   'offrep': offreps(S, N9),
                                   'realign': corpus_realign(S, N9)}, 'candidates': {}}
    print(f"REAL offrep {len(report['REAL']['offrep'])}  realign {report['REAL']['realign']}")
    for lab in cf:
        rr = [rows[lab][sd] for sd in seeds]
        sm = summarize(rr, real_dep)
        rs = [r['sig']['resync'] for r in rr]
        sm['resync_sd'] = round(statistics.pstdev(rs), 3)
        sm['resync_zero_frac'] = round(sum(1 for x in rs if x == 0) / len(rs), 2)
        band = [0.03 <= x <= 0.12 for x in rs]
        no_rep = [r['n_offrep'] == 0 for r in rr]
        sm['offrep_mean'] = round(statistics.mean(r['n_offrep'] for r in rr), 2)
        sm['offrep_zero_frac'] = round(sum(no_rep) / len(rr), 2)
        sm['resync_band_and_no_offrep'] = round(sum(x and y for x, y in zip(band, no_rep)) / len(rr), 3)
        re = [r['resync_edits'] for r in rr]
        sm['resync_edits_mean'] = round(statistics.mean(re), 3)
        sm['resync_edits_band'] = round(sum(0.03 <= x <= 0.12 for x in re) / len(re), 2)
        report['candidates'][lab] = sm
        pr = sm['pass_rate']
        print(f"{lab:26s} " + ' '.join(f'{v:.2f}' for v in pr.values()) +
              f" | mean {sm['mean_tests_passed']} all7 {sm['all7_rate']:.3f}"
              f" rs {sm['mean_resync']}±{sm['resync_sd']} (0:{sm['resync_zero_frac']})"
              f" r2 {sm['mean_r2']} r3 {sm['mean_r3']} comm {sm['mean_commute']}"
              f" fix {sm['mean_fixed_frac']}\n{'':26s} offrep {sm['offrep_mean']} (0:{sm['offrep_zero_frac']})"
              f" band&no-offrep {sm['resync_band_and_no_offrep']}"
              f" | short edits: rs {sm['resync_edits_mean']} band {sm['resync_edits_band']}", flush=True)
    print('columns: ' + ' | '.join(pr))
    out = os.path.join(ROOT, 'data', 'tournament5_report.json')
    with open(out, 'w') as f:
        json.dump(report, f, indent=1)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
