#!/usr/bin/env python3
"""Tournament round 4: can gak41 get its missing resync and near-repeat depletion?

Round 3 left gak41 (affine autokey over C83:C41) as the best candidate but
failing two tests structurally:
  * resync = 0. A group autokey's state update g <- g o m(p) is INVERTIBLE, so
    after two messages diverge their state difference is only ever conjugated,
    never cancelled: they can never re-agree. The real corpus re-agrees
    (rate ~0.065), so the real state update must LOSE information somewhere.
  * d2/d3 not depleted (r2 ~1.7 vs real 0.40).

Two minimal modifications, each testable against those failures:
  reset  -- on a trigger plaintext letter the state returns to the identity
            (b, a) = (1, a0). A reset is the simplest lossy update: after it,
            two diverged messages agree again. Triggers must be rare enough to
            leave >=25-letter isomorph windows unreset.
  avoid  -- if the output D[a] repeats one of the last w outputs, bump a by a
            fixed step until it doesn't ('state': the bump is kept in the
            state; 'out': only the emitted letter is bumped).
  short  -- additive steps drawn from 1..K (K < 83): a d-repeat needs the
            steps to sum to 0 mod 83, impossible on untwisted stretches while
            dK < 83, so d2/d3 are depleted softly instead of zeroed.
and their combination. Same harness, tests and 12 seeds as tournament3; the
leading configs are then re-scored over 48 seeds ('confirm_48').

Writes data/tournament4_report.json. Run from repo root: python3 tools/tournament4.py
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tournament import P, PA, WEIGHTS, make_gak41, load_msgs, N9, signatures  # noqa: E402
from tournament3 import depletion, tests, run_config, summarize, SEEDS  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FREQ = {c: w / sum(WEIGHTS) for c, w in zip(PA, WEIGHTS)}


def make_gak41x(seed, trig=(), w=0, mode='state', step=None, K=None):
    """gak41 (identical key material to tournament.make_gak41 for the same seed)
    plus optional reset triggers and an avoidance window of width w.
    K: draw the additive steps a_of as DISTINCT values from 1..K (short steps):
    on untwisted stretches a_t - a_(t-d) is a sum of d steps in [d, dK], so it
    cannot be 0 mod 83 while dK < 83 -- d-repeats vanish structurally there and
    only the twist letters (b != 1) let them back in (soft depletion)."""
    rng = random.Random(seed + 9)
    sub = sorted({pow(4, k, P) for k in range(41)})
    twist = set(rng.sample(PA, 5))
    a_of = {c: rng.randrange(1, P) for c in PA}
    if K:
        a_of = dict(zip(PA, random.Random(seed * 7 + K).sample(range(1, K + 1), len(PA))))
    b_of = {c: rng.choice([x for x in sub if x != 1]) if c in twist else 1
            for c in PA}
    D = list(range(P)); rng.shuffle(D)
    st = step or random.Random(seed * 31 + 5).randrange(1, P)
    trig = set(trig)

    def enc(s):
        b, a = 1, 17
        seq = [D[(50 + a_of[s[0]]) % P]]
        for ch in s[1:]:
            if ch in trig:
                b, a = 1, 17
            b, a = (b * b_of[ch]) % P, (b * a_of[ch] + a) % P
            o = a
            if w:
                recent = seq[-w:]
                k = 0
                while D[o] in recent and k < P:
                    o = (o + st) % P; k += 1
                if mode == 'state':
                    a = o
            seq.append(D[o])
        return seq
    return enc


def trig_sets():
    """Trigger sets spanning plaintext frequency mass ~0.2% .. ~15%."""
    by = sorted(PA, key=lambda c: FREQ[c])
    return {
        'qz': ['q', 'z'],
        'rare4': by[:4],                       # q z j x
        'kvb': ['k', 'v', 'b'],
        'y': ['y'],
        'm': ['m'],
        'n': ['n'],
        'e': ['e'],
        'space': ['_'],
    }


def main():
    S = load_msgs()
    real_sig = signatures(S, N9, random.Random(1))
    real_dep = depletion(S, N9)
    report = {'seeds': SEEDS, 'REAL': {'sig': real_sig, 'dep': real_dep,
                                       'score': tests(real_sig, real_dep, real_dep)}}
    print(f'REAL: {real_sig}\nREAL depletion: {real_dep}')

    configs = {'gak41 (baseline)': lambda sd: make_gak41(sd)}
    ts = trig_sets()
    for name, t in ts.items():
        configs[f'reset {name}'] = (lambda t: lambda sd: make_gak41x(sd, trig=t))(t)
    for w in (1, 2, 3):
        for mode in ('state', 'out'):
            configs[f'avoid w={w} {mode}'] = \
                (lambda w, mode: lambda sd: make_gak41x(sd, w=w, mode=mode))(w, mode)
    for name in ('kvb', 'y', 'm', 'n', 'space'):
        for w in (2, 3):
            for mode in ('state', 'out'):
                configs[f'reset {name} + avoid w={w} {mode}'] = \
                    (lambda t, w, mode: lambda sd: make_gak41x(sd, trig=t, w=w, mode=mode))(
                        ts[name], w, mode)

    for K in (27, 41, 55):
        configs[f'short K={K}'] = (lambda K: lambda sd: make_gak41x(sd, K=K))(K)
        for name in ('qz', 'rare4', 'y'):
            configs[f'short K={K} + reset {name}'] = \
                (lambda K, t: lambda sd: make_gak41x(sd, K=K, trig=t))(K, ts[name])
        configs[f'short K={K} + avoid w=1 out'] = \
            (lambda K: lambda sd: make_gak41x(sd, K=K, w=1, mode='out'))(K)
        configs[f'short K={K} + rare4 + avoid w=1 out'] = \
            (lambda K: lambda sd: make_gak41x(sd, K=K, w=1, mode='out', trig=ts['rare4']))(K)

    report['trigger_mass'] = {n: round(sum(FREQ[c] for c in t), 4) for n, t in ts.items()}
    print('trigger mass:', report['trigger_mass'])
    res = {}
    for label, mk in configs.items():
        sm = summarize(run_config(label, mk), real_dep)
        res[label] = sm
        pr = sm['pass_rate']
        print(f"  {label:34s} " + ' '.join(f'{v:.2f}' for v in pr.values()) +
              f"  | mean {sm['mean_tests_passed']} all7 {sm['all7_rate']:.2f}"
              f" maxL {sm['mean_maxL']} comm {sm['mean_commute']} rs {sm['mean_resync']}"
              f" r2 {sm['mean_r2']} r3 {sm['mean_r3']}", flush=True)
    print('  columns: ' + ' | '.join(next(iter(res.values()))['pass_rate']))
    report['candidates'] = res

    # re-score the leaders over 48 seeds: 2/12 all-7 passes is within noise
    import tournament3
    saved = list(tournament3.SEEDS)
    tournament3.SEEDS[:] = list(range(1, 49))
    lead = ['gak41 (baseline)', 'short K=27', 'short K=27 + reset qz',
            'short K=27 + reset rare4', 'reset m + avoid w=3 state']
    configs['reset rare4 + avoid w=3 out'] = \
        lambda sd: make_gak41x(sd, trig=ts['rare4'], w=3, mode='out')
    lead.append('reset rare4 + avoid w=3 out')
    report['confirm_48'] = {}
    print('\n48-seed confirmation:')
    for label in lead:
        sm = summarize(run_config(label, configs[label]), real_dep)
        report['confirm_48'][label] = sm
        print(f"  {label:34s} mean {sm['mean_tests_passed']} all7 {sm['all7_rate']:.2f}"
              f" rs {sm['mean_resync']} (pass {sm['pass_rate']['resync .03-.12']})"
              f" r2 {sm['mean_r2']} r3 {sm['mean_r3']}", flush=True)
    tournament3.SEEDS[:] = saved
    out = os.path.join(ROOT, 'data', 'tournament4_report.json')
    with open(out, 'w') as f:
        json.dump(report, f, indent=1)
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
