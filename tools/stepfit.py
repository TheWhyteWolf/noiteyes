#!/usr/bin/env python3
"""Step G: key fit for the round-5 leader (lossless short-step gak41).

Model: c_t = D[a_t], a_t = a_(t-1) + b_(t-1) * k(p_t) mod 83 with steps k in
1..27 and b in the order-41 subgroup, changed only at 2-3 twist letters by
that letter's fixed multiplier. With pi = D^-1 every consecutive difference
pi(c_t) - pi(c_(t-1)) must be b * (a step in 1..K), and b may change only
where that step is a twist letter's step. If pi and the twists were found,
each position's step identifies its plaintext letter (a monoalphabetic
readout).

tools/stepfit.c anneals pi (and the twist (step, multiplier) pairs) against a
Viterbi cost over b: 1 per position whose step falls outside 1..K, mu per
unexplained jump of b. This driver runs, in order:
  landscape -- a synthetic corpus of the model (tournament5.make_lossy, 2
               twists): cost of the TRUE pi, and of the truth after k random
               transpositions
  gate      -- annealing on that synthetic corpus, twists given and annealed;
               PASS only if it reaches the truth's cost and recovers pi
  real      -- the real corpus, the real corpus reversed (a short-step order
               cannot fit backwards: -1 is not a quadratic residue mod 83),
               and a full-step gak41 null; best reachable costs compared
The real verdict is read only if the gate passes.

Result: the gate FAILS -- see data/stepfit_report.json and README Step G.

Needs gcc; heavy (~1 h on 16 cores), so run it on enter:
    python3 tools/stepfit.py [iters] [restarts]
"""
import json
import os
import random
import subprocess
import sys
import time
from multiprocessing.pool import ThreadPool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tournament import corpus, load_msgs, N9, make_gak41  # noqa: E402
from tournament5 import make_lossy  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join('/tmp', f'stepfit-{os.getuid()}')
K, MU, T0, T1 = 27, 4, 8.0, 0.3


def build():
    subprocess.run(['gcc', '-O3', '-march=native', '-o', BIN,
                    os.path.join(ROOT, 'tools', 'stepfit.c'), '-lm'], check=True)


def feed(S, names, iters, seed, tw, pi=None):
    inp = (f'{K} {MU} {iters} {seed} {T0} {T1}\n{len(tw)} ' +
           ' '.join(f'{a} {b}' for a, b in tw) + f'\n{len(names)}\n' +
           ''.join(f'{len(S[n])} ' + ' '.join(map(str, S[n])) + '\n' for n in names))
    if pi is not None:
        inp += ' '.join(map(str, pi))
    out = subprocess.run([BIN], input=inp, capture_output=True, text=True,
                         check=True).stdout.split('\n')
    return int(out[0]), (list(map(int, out[1].split())) if pi is None else None), \
        (out[2].strip() if pi is None else None)


def agreement(pit, pi):
    """Symbols on which pi matches the truth up to affine gauge u*pi+v."""
    return max(sum(1 for x in range(83) if (u * pit[x] + v) % 83 == pi[x])
               for u in range(1, 83) for v in range(83))


def main():
    iters = int(sys.argv[1]) if len(sys.argv) > 1 else 3_000_000
    restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    build()
    enc = make_lossy(1, ntw=2)
    S, names = corpus(enc, random.Random(1001))
    pit = [0] * 83
    for i, c in enumerate(enc.D):
        pit[c] = i
    tw = enc.twist
    rep = {'params': {'K': K, 'mu': MU, 'T0': T0, 'T1': T1, 'iters': iters,
                      'restarts': restarts}, 'synthetic_twist': tw}

    # landscape around the truth
    rr = random.Random(3)
    land = {'truth': feed(S, names, 0, 1, tw, pit)[0],
            'random': feed(S, names, 0, 1, tw, rr.sample(range(83), 83))[0]}
    for k in (1, 2, 5, 10, 20, 40):
        cs = []
        for _ in range(5):
            pi = list(pit)
            for _ in range(k):
                x, y = rr.sample(range(83), 2); pi[x], pi[y] = pi[y], pi[x]
            cs.append(feed(S, names, 0, 1, tw, pi)[0])
        land[f'{k} swaps'] = cs
    rep['landscape'] = land
    print('landscape:', land, flush=True)

    R = load_msgs()
    rev = {n: [R[n][0]] + R[n][1:][::-1] for n in N9}
    full, fn = corpus(make_gak41(1), random.Random(1001))
    none = [(0, 0)] * 2
    jobs = []
    for s in range(restarts):
        jobs += [('gate, twists given', S, names, tw, s), ('gate, twists annealed', S, names, none, s),
                 ('real', R, N9, none, s), ('real reversed', rev, N9, none, s),
                 ('full-step gak41 null', full, fn, none, s)]

    def run(j):
        lab, SS, nn, t, s = j
        c, pi, twout = feed(SS, nn, iters, s, t)
        ag = agreement(pit, pi) if lab.startswith('gate') else None
        return lab, c, ag, twout

    t = time.time()
    with ThreadPool(os.cpu_count()) as pool:
        res = pool.map(run, jobs, chunksize=1)
    rep['anneal'] = {}
    for lab in dict.fromkeys(j[0] for j in jobs):
        rows = sorted((c, ag, tw_) for L, c, ag, tw_ in res if L == lab)
        SS, nn = next((j[1], j[2]) for j in jobs if j[0] == lab)
        npos = sum(len(SS[n]) - 2 for n in nn)
        rep['anneal'][lab] = {'positions': npos, 'costs': [r[0] for r in rows],
                              'best_per_position': round(rows[0][0] / npos, 3),
                              'agreement': [r[1] for r in rows] if rows[0][1] is not None else None,
                              'best_twists': rows[0][2]}
        print(lab, rep['anneal'][lab], flush=True)
    g = rep['anneal']['gate, twists given']
    rep['gate_pass'] = min(g['costs']) <= land['truth'] and max(g['agreement']) >= 75
    rep['seconds'] = round(time.time() - t)
    print('GATE', 'PASS' if rep['gate_pass'] else 'FAIL', f"({rep['seconds']} s)")
    with open(os.path.join(ROOT, 'data', 'stepfit_report.json'), 'w') as f:
        json.dump(rep, f, indent=1)


if __name__ == '__main__':
    main()
