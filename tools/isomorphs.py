#!/usr/bin/env python3
"""Dump the isomorph catalog to data/isomorphs.json.

Uses the finder from pi_solver.py (seed length 12, >=3 repeated letters,
maximal extension). Isomorphism (identical repeat pattern in both windows)
holds over each pair's full extent by construction; the evidential caveat is
that positions whose letter occurs only ONCE in the window contribute no
repeat-structure evidence, so maximal extension can overrun the true shared
phrase through such positions (see PLAN.md). Each entry therefore records
evidence_positions (positions with a letter repeated within the window) next
to the raw length — treat length as an upper bound on the shared phrase and
evidence_positions as its support.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter

from pi_solver import D, find_isomorph_pairs, load_msgs


def main():
    S = load_msgs()
    out = []
    for n1, l1, r1, n2, l2, r2 in sorted(find_isomorph_pairs(S),
                                         key=lambda p: -(p[2] - p[1])):
        w1, w2 = S[n1][l1:r1], S[n2][l2:r2]
        cnt = Counter(w1)
        out.append({'a': {'msg': n1, 'start': l1, 'end': r1},
                    'b': {'msg': n2, 'start': l2, 'end': r2},
                    'length': r1 - l1,
                    'evidence_positions': sum(1 for v in w1 if cnt[v] >= 2),
                    'ascii_a': ''.join(chr(32 + v) for v in w1),
                    'ascii_b': ''.join(chr(32 + v) for v in w2)})
    json.dump({'finder': 'pi_solver.find_isomorph_pairs (L0=12, min_rep=3, maximal)',
               'pairs': out},
              open(f'{D}/data/isomorphs.json', 'w'), indent=1)
    print(f'wrote data/isomorphs.json: {len(out)} maximal pairs')
    for p in out:
        print(f'  L={p["length"]:3d} ev={p["evidence_positions"]:2d} '
              f'{p["a"]["msg"]}[{p["a"]["start"]}:{p["a"]["end"]}] ~ '
              f'{p["b"]["msg"]}[{p["b"]["start"]}:{p["b"]["end"]}]')


if __name__ == '__main__':
    sys.exit(main())
