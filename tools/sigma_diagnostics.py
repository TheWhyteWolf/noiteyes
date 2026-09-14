#!/usr/bin/env python3
"""Sigma-composition diagnostics on the isomorph catalog (PLAN.md item 1).

For every maximal isomorph pair (msgA[l1:r1] ~ msgB[l2:r2]) the aligned
position offset delta = l2-l1 is constant across the whole window (by
construction: the finder extends both windows in lock-step). That lets us
collapse "is a real relation between msgA and msgB" down to a single key
(msgA, msgB, delta), and build sigma_{msgA,msgB,delta}: value-at-msgA ->
value-at-msgB, a partial function on the 83-letter alphabet, from the union
of every window sharing that key.

If the strict additive-walk model (c_i = D[x_i], x_i = x_{i-1}+v(p_i), D a
fixed bijection) were correct, every such sigma would equal D o (+delta) o
D^-1 for its own delta -- i.e. all sigmas are conjugates of shift maps by
the SAME D, so composition must hold exactly:

    sigma_{B,C,dBC} ( sigma_{A,B,dAB} (v) )  ==  sigma_{A,C,dAB+dBC} (v)

whenever both sides are defined. This script builds every pairwise sigma
from data/messages.json + the isomorph finder, checks each sigma is
internally well-defined (function AND injective) when multiple windows are
merged into one key, then exhaustively tests composition over every
A-B-C triple reachable through the catalog. It also decomposes same-message
self-relations (msgA == msgB, delta != 0) into cycles to characterize their
algebraic shape directly.

This is a strictly local/combinatorial test -- no linear algebra, no
assumption of a shared D beyond what composition itself requires -- so it
localizes *where* pi_solver.py's global refutation comes from: if particular
triangles/components compose cleanly while others don't, the break is
structural (resync/skip/etc.), not a global accounting error.
"""
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pi_solver import D, find_isomorph_pairs, load_msgs  # noqa: E402

GROUP_A = {'east-1', 'east-2', 'west-1'}
GROUP_B = {'east-4', 'east-5', 'west-4'}


def group_of(msg):
    if msg in GROUP_A:
        return 'A'
    if msg in GROUP_B:
        return 'B'
    return '?'


def build_windows(pairs):
    """key (n1, n2, delta) -> list of (l1, r1) windows in n1's coordinates."""
    windows = defaultdict(list)
    for n1, l1, r1, n2, l2, r2 in pairs:
        d = l2 - l1
        windows[(n1, n2, d)].append((l1, r1))
        windows[(n2, n1, -d)].append((l2, r2))
    return windows


def build_sigmas(S, windows):
    sigmas = {}
    report = []
    for key, wins in windows.items():
        n1, n2, d = key
        fwd, inv = {}, {}
        fn_conflicts = inj_conflicts = 0
        npos = 0
        for l1, r1 in wins:
            for p in range(l1, r1):
                p2 = p + d
                if not (0 <= p2 < len(S[n2])):
                    continue
                v, w = S[n1][p], S[n2][p2]
                npos += 1
                if v in fwd and fwd[v] != w:
                    fn_conflicts += 1
                else:
                    fwd[v] = w
                if w in inv and inv[w] != v:
                    inj_conflicts += 1
                else:
                    inv[w] = v
        sigmas[key] = fwd
        report.append((key, len(wins), npos, len(fwd), fn_conflicts, inj_conflicts))
    return sigmas, report


def test_composition(sigmas):
    results = []
    for kAB, sigAB in sigmas.items():
        nA, nB, dAB = kAB
        for kBC, sigBC in sigmas.items():
            nB2, nC, dBC = kBC
            if nB2 != nB:
                continue
            dAC = dAB + dBC
            if nA == nC and dAC == 0:
                continue  # trivial: inverse-pair round trip, guaranteed by construction
            kAC = (nA, nC, dAC)
            if kAC not in sigmas:
                continue
            sigAC = sigmas[kAC]
            agree, disagree, examples = 0, 0, []
            for v, w in sigAB.items():
                if w not in sigBC:
                    continue
                u = sigBC[w]
                if v not in sigAC:
                    continue
                if sigAC[v] == u:
                    agree += 1
                else:
                    disagree += 1
                    examples.append((v, sigAC[v], u))
            if agree + disagree > 0:
                results.append({
                    'A': nA, 'B': nB, 'C': nC, 'dAB': dAB, 'dBC': dBC, 'dAC': dAC,
                    'agree': agree, 'disagree': disagree, 'examples': examples[:5],
                })
    return results


def cycle_decompose(perm):
    """perm: dict value->value (partial, domain subset of range assumed for
    the cycle to close). Returns list of cycles (as tuples) using only
    elements whose forward chain returns to itself within the given dict."""
    seen = set()
    cycles = []
    for start in perm:
        if start in seen:
            continue
        chain = [start]
        cur = perm[start]
        ok = True
        while cur != start:
            if cur in seen or cur not in perm or cur in chain:
                ok = False
                break
            chain.append(cur)
            cur = perm[cur]
        for x in chain:
            seen.add(x)
        if ok:
            cycles.append(tuple(chain))
    return cycles


def main():
    S = load_msgs()
    pairs = find_isomorph_pairs(S)
    print(f'isomorph pairs (maximal, L0=12 seed): {len(pairs)}')

    windows = build_windows(pairs)
    sigmas, report = build_sigmas(S, windows)

    print(f'\ndistinct (msgA,msgB,delta) relations: {len(sigmas)} '
          f'({len(sigmas)//2} unordered pairs)')
    print('\nper-relation well-definedness (should all be 0 conflicts):')
    bad_relations = 0
    for key, nwin, npos, dom, fn_c, inj_c in sorted(report, key=lambda r: -r[3]):
        n1, n2, d = key
        flag = '' if fn_c == 0 and inj_c == 0 else '  <-- NOT WELL-DEFINED'
        if flag:
            bad_relations += 1
        print(f'  {n1:8s}->{n2:8s} delta={d:+4d}  windows={nwin}  positions={npos:3d} '
              f'domain={dom:2d}  fn_conflicts={fn_c}  inj_conflicts={inj_c}{flag}')
    print(f'\n{bad_relations}/{len(sigmas)} relations are internally inconsistent '
          f'as merged (multiple windows disagree on the same key).')

    print('\n--- self-relations (same message, delta != 0): cycle structure ---')
    for (n1, n2, d), sig in sigmas.items():
        if n1 != n2:
            continue
        cycles = cycle_decompose(sig)
        lens = sorted(len(c) for c in cycles)
        covered = sum(lens)
        print(f'  {n1} delta={d:+4d}: domain={len(sig)}, closed cycles cover {covered} '
              f'letters, lengths={lens}')
        for c in cycles:
            if len(c) > 1:
                print(f'      cycle len {len(c)}: {c}')

    print('\n--- composition test: sigma_BC(sigma_AB(v)) vs sigma_AC(v) ---')
    results = test_composition(sigmas)
    tot_agree = tot_disagree = 0
    clean = broken = 0
    for r in sorted(results, key=lambda r: -(r['agree'] + r['disagree'])):
        tot_agree += r['agree']
        tot_disagree += r['disagree']
        gA, gB, gC = group_of(r['A']), group_of(r['B']), group_of(r['C'])
        status = 'CLEAN' if r['disagree'] == 0 else f"BROKEN x{r['disagree']}"
        if r['disagree'] == 0:
            clean += 1
        else:
            broken += 1
        print(f"  [{gA}{gB}{gC}] {r['A']}--{r['dAB']:+d}-->{r['B']}--{r['dBC']:+d}-->{r['C']} "
              f"(direct delta {r['dAC']:+d}): agree={r['agree']:3d} disagree={r['disagree']:3d} "
              f"{status}")
        if r['disagree']:
            print(f"      example (v, direct_sigma_AC[v], composed): {r['examples']}")

    print(f'\ntriples tested: {len(results)}  (clean={clean}, broken={broken})')
    print(f'letter-instance comparisons: agree={tot_agree}, disagree={tot_disagree}')
    if tot_agree + tot_disagree:
        print(f'agreement rate: {tot_agree/(tot_agree+tot_disagree):.4f}')

    out = {
        'relations': [{'msgA': k[0], 'msgB': k[1], 'delta': k[2], 'domain_size': len(v)}
                      for k, v in sigmas.items()],
        'composition_tests': results,
    }
    path = f'{D}/data/sigma_diagnostics.json'
    json.dump(out, open(path, 'w'), indent=1)
    print(f'\nsaved {path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
