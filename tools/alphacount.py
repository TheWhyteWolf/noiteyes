#!/usr/bin/env python3
"""Minimum-alphabet factorization of the sigma web.

rotorfit left the surviving family as: sigma between two isomorph occurrences
is a fixed PER-OCCURRENCE substitution difference sigma = P_j o P_i^-1, with
no positional (rotor/progressive) organization. That leaves one question the
ciphertext can still answer on its own: HOW MANY distinct substitution
alphabets P are needed to explain the whole sigma web?

Model: each isomorph occurrence i carries one unknown bijection P_i
(plaintext->cipher over its window); an edge (i,j) with letter-map sigma_ij
constrains sigma_ij = P_j o P_i^-1 on its observed support. Occurrences may
be grouped into classes sharing one alphabet. A candidate partition is
CONSISTENT iff class-level relations rho_cd = Q_d o Q_c^-1, seeded from the
edges and closed under inverse and composition, stay (a) functional,
(b) injective, and (c) identity on the diagonal (rho_cc).

The tool computes, per connected web:
  * base ("all-singletons") closure -- violations here mean the
    one-alphabet-per-occurrence model itself fails inside some window;
  * the pairwise verdict for every occurrence pair: DISTINCT (some derived
    entry contradicts identity), SUPPORTED-EQUAL (nonempty derived relation,
    all identity), or FREE (no derived constraint at all);
  * m_min = minimum number of classes over all consistent partitions
    (exact branch-and-bound; the webs are tiny), plus every minimal
    partition (capped);
  * cross-checks: family membership and start-position mod 4 per class.

m_min is a LOWER bound on the alphabet count the data forces: FREE pairs can
always merge, so sparse webs under-determine the schedule -- quantifying
exactly that is the point.

Control gates (all must PASS before the real verdict is read):
  quag4  -- pool of 4 alphabets on a known schedule, phrases planted so every
            web edge is direct: recovered partition must EQUAL ground truth
            (accepts truth, m_min == true m, unique minimal partition).
  fresh  -- a fresh alphabet per stretch: truth (all-singleton) accepted,
            m_min == web size (no false merges).
  transl -- alphabets are translates of one deck (the refuted-walk shape):
            same requirements as fresh; validates that structured, commuting
            sigmas do not fool the closure.
"""
import json
import os
import random
import sys
from collections import Counter, deque

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
from pi_solver import P, find_isomorph_pairs, load_msgs
from sigma_web import MINLEN, PA, ROOT, WEIGHTS, build_webs, witness_trim

SEED = 1
FAM48 = {'east-1', 'east-2', 'west-1'}


# ------------------------------------------------- relation closure engine

def close_relations(seeds, labels):
    """seeds: list of (i, j, sigma_dict) occurrence-level edges.
    labels: occ -> class id. Returns (rels, violations); rels maps
    (c, d) -> {x: y} partial bijection, closed under inverse+composition.
    Violations are (kind, c, d, x, y1, y2) tuples; any violation means the
    labeling is inconsistent with per-class alphabets."""
    rels, rev, viol = {}, {}, []
    queue = deque()
    for i, j, sg in seeds:
        c, d = labels[i], labels[j]
        for x, y in sg.items():
            queue.append((c, d, x, y))
            queue.append((d, c, y, x))
    while queue:
        c, d, x, y = queue.popleft()
        if c == d and x != y:
            viol.append(('diag', c, d, x, y, None))
            continue
        m = rels.setdefault((c, d), {})
        r = rev.setdefault((c, d), {})
        if x in m:
            if m[x] != y:
                viol.append(('func', c, d, x, m[x], y))
            continue
        if y in r:
            viol.append(('inj', c, d, y, r[y], x))
            continue
        m[x] = y
        r[y] = x
        # right-compose: (c,d,x,y) + (d,e,y,z) -> (c,e,x,z)
        for (dd, e), m2 in list(rels.items()):
            if dd == d and y in m2:
                queue.append((c, e, x, m2[y]))
        # left-compose: (e,c,w,x) + (c,d,x,y) -> (e,d,w,y)
        for (e, cc), _m2 in list(rels.items()):
            if cc == c:
                w = rev[(e, cc)].get(x)
                if w is not None:
                    queue.append((e, d, w, y))
    return rels, viol


def consistent(seeds, labels):
    return not close_relations(seeds, labels)[1]


# ------------------------------------------------------ partition search

def min_partitions(nodes, seeds, cap=2000):
    """Exact minimum class count over consistent partitions of `nodes`,
    via restricted-growth branch-and-bound; returns (m_min, partitions)
    where partitions is every minimal partition (as tuples of frozensets),
    capped at `cap`."""
    n = len(nodes)
    if n == 0:
        return 0, []
    # order: most-constrained first (nodes on many edges)
    deg = {v: 0 for v in nodes}
    for i, j, _ in seeds:
        if i in deg:
            deg[i] += 1
        if j in deg:
            deg[j] += 1
    order = sorted(nodes, key=lambda v: -deg[v])

    best = [n + 1]
    found = []

    def rec(idx, classes):
        if len(classes) > best[0]:
            return
        if idx == n:
            m = len(classes)
            if m < best[0]:
                best[0] = m
                found.clear()
            if m == best[0] and len(found) < cap:
                found.append([frozenset(c) for c in classes])
            return
        v = order[idx]
        for c in classes:
            c.add(v)
            lab = {}
            for ci, cl in enumerate(classes):
                for u in cl:
                    lab[u] = ci
            sub = [(i, j, sg) for i, j, sg in seeds if i in lab and j in lab]
            if consistent(sub, lab):
                rec(idx + 1, classes)
            c.discard(v)
        if len(classes) + 1 <= best[0]:
            classes.append({v})
            rec(idx + 1, classes)
            classes.pop()

    rec(0, [])
    return best[0], found


def pairwise_verdicts(nodes, seeds):
    """Classify every occurrence pair from the base (all-singleton) closure:
    DISTINCT / SUPPORTED-EQUAL (with evidence count) / FREE."""
    labels = {v: v for v in nodes}
    rels, viol = close_relations(seeds, labels)
    out = {}
    for a in nodes:
        for b in nodes:
            if a >= b:
                continue
            m = rels.get((a, b), {})
            if not m:
                out[(a, b)] = ('FREE', 0)
            elif all(x == y for x, y in m.items()):
                out[(a, b)] = ('SUPPORTED-EQUAL', len(m))
            else:
                bad = sum(1 for x, y in m.items() if x != y)
                out[(a, b)] = ('DISTINCT', bad)
    return out, viol


# ------------------------------------------------------- corpus pipeline

def core_sigma(S, ent):
    """Sigma from the pair's LARGEST strictly-overlapping span block only,
    further restricted to repeat-evidenced positions inside that block.

    Rationale: maximal windows overrun the true shared phrase through
    singleton letters, and a mirrored carrier repeat in the overrun forms its
    own small span DISCONNECTED (strict-overlap-wise) from the phrase's span
    chain. witness_trim's [first span .. last span] extent keeps such junk;
    the largest contiguous block does not. Within one block every position is
    repeat-protected, so alphabet changes cannot hide inside it without
    breaking pattern equality (any residue is caught by the base-closure
    violation check)."""
    n1, l1, r1, n2, l2, r2 = ent['raw']
    s, e = max(ent['spans'], key=lambda sp: sp[1] - sp[0])
    w1, w2 = S[n1][l1 + s:l1 + e], S[n2][l2 + s:l2 + e]
    cnt = Counter(w1)
    sg = {a: b for a, b in zip(w1, w2) if cnt[a] >= 2}
    return sg, (n1, l1 + s, l1 + e, n2, l2 + s, l2 + e)


def web_structures(S, names):
    """Shared pipeline (identical for synthetic and real): isomorph catalog
    -> witness trim -> largest-block core sigmas -> occurrence webs."""
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    joint = witness_trim(S, raw)
    sigmas, cores = [], []
    for ent in joint:
        sg, core = core_sigma(S, ent)
        if core[2] - core[1] >= 6:
            sigmas.append(sg)
            cores.append(core)
    occs, pot, webs, gv, edges = build_webs(cores)
    seeds = [(a, b, sigmas[k]) for a, b, off, k in edges if a != b]
    return occs, webs, seeds


# --------------------------------------------------- synthetic controls

SL, NS, MLEN = 30, 4, 120           # stretch length, stretches, message length


def rep_phrase(rng, L):
    pool = rng.sample(PA, 6)
    return [pool[rng.randrange(6)] for _ in range(L)]


# (message, stretch, in-stretch offset) placements per phrase; designed so
# every web mixes same-alphabet and different-alphabet stretch pairs under
# the quag4 assignment asn[msg][stretch] = (msg + stretch) % 4.
PLACEMENTS = [
    [(0, 0, 4), (1, 1, 6), (2, 2, 5), (3, 0, 7)],   # quag4 alphs {0,2,0,3}
    [(4, 1, 4), (5, 2, 6), (6, 3, 5), (7, 1, 7)],   # quag4 alphs {1,3,1,0}
    [(8, 0, 4), (0, 2, 6), (1, 3, 5)],              # quag4 alphs {0,2,0}
]


def synth(mode, seed):
    """Returns (S, names, alph_of) where alph_of(msg_idx, pos) is the ground
    truth alphabet index at that position."""
    rng = random.Random(seed)
    names = [f'sim-{i}' for i in range(9)]
    if mode == 'quag4':
        asn = [[(i + j) % 4 for j in range(NS)] for i in range(9)]
        nal = 4
    else:
        asn = [[i * NS + j for j in range(NS)] for i in range(9)]
        nal = 9 * NS
    if mode == 'transl':
        deck = list(range(P))
        rng.shuffle(deck)
        base = dict(zip(PA, rng.sample(range(P), len(PA))))
        shifts = rng.sample(range(P), nal)
        alphs = [{c: deck[(base[c] + s) % P] for c in PA} for s in shifts]
    else:
        alphs = [dict(zip(PA, rng.sample(range(P), len(PA)))) for _ in range(nal)]
    plain = {nm: rng.choices(PA, weights=WEIGHTS, k=MLEN) for nm in names}
    for ph_i, places in enumerate(PLACEMENTS):
        ph = rep_phrase(rng, 22)
        for mi, st, off in places:
            pos = st * SL + off
            plain[names[mi]][pos:pos + len(ph)] = ph
    S = {}
    for i, nm in enumerate(names):
        S[nm] = [alphs[asn[i][p // SL]][c] for p, c in enumerate(plain[nm])]

    def alph_of(msg_idx, pos):
        return asn[msg_idx][pos // SL]
    return S, names, alph_of


def truth_labels(occs, names, alph_of):
    """Ground-truth alphabet index per occurrence; None if the (trimmed,
    merged) window straddles a stretch boundary."""
    out = []
    for msg, lo, hi in occs:
        mi = names.index(msg)
        ids = {alph_of(mi, p) for p in range(lo, hi)}
        out.append(ids.pop() if len(ids) == 1 else None)
    return out


def run_control(mode):
    S, names, alph_of = synth(mode, SEED)
    occs, webs, seeds = web_structures(S, names)
    truth = truth_labels(occs, names, alph_of)
    res = {'mode': mode, 'n_occs': len(occs), 'webs': [], 'ok': True,
           'n_straddling': sum(1 for t in truth if t is None)}
    checks = {'accept': True, 'm_min': True, 'exact': True}
    for web in webs:
        # occurrences whose window straddles a stretch boundary have no single
        # ground-truth alphabet; they are excluded from the gates (reported).
        nodes = [o for o in web if truth[o] is not None]
        sub = [(i, j, sg) for i, j, sg in seeds if i in nodes and j in nodes]
        # gate 1: ground-truth partition is consistent
        tl = {o: truth[o] for o in nodes}
        acc = consistent(sub, tl)
        checks['accept'] &= acc
        # gate 2: exact minimum matches true distinct count
        m_true = len({truth[o] for o in nodes})
        m_min, parts = min_partitions(nodes, sub)
        checks['m_min'] &= (m_min == m_true)
        # gate 3 (quag4): the minimal partition is unique and equals truth
        exact = True
        if mode == 'quag4':
            tp = {}
            for o in nodes:
                tp.setdefault(truth[o], set()).add(o)
            tset = {frozenset(v) for v in tp.values()}
            exact = len(parts) == 1 and set(parts[0]) == tset
            checks['exact'] &= exact
        res['webs'].append({'nodes': nodes, 'm_true': m_true, 'm_min': m_min,
                            'n_minimal': len(parts), 'accepts_truth': acc,
                            'exact': exact})
    res['checks'] = checks
    res['ok'] = all(checks.values())
    return res


def gate(name, ok):
    print(f'  {"PASS" if ok else "FAIL":4s}  {name}')
    return ok


# --------------------------------------------------------------- real run

def run_real():
    S = load_msgs()
    occs, webs, seeds = web_structures(S, pi_solver.N9)
    report = {'occurrences': [], 'webs': []}
    for k, (msg, lo, hi) in enumerate(occs):
        report['occurrences'].append(
            {'id': k, 'msg': msg, 'lo': lo, 'hi': hi,
             'family': 48 if msg in FAM48 else 49, 'start_mod4': lo % 4})
    print(f'  occurrences: {len(occs)}, webs: '
          f'{[len(w) for w in webs]}, edges: {len(seeds)}')
    for wi, web in enumerate(webs):
        sub = [(i, j, sg) for i, j, sg in seeds if i in web and j in web]
        pv, viol = pairwise_verdicts(web, sub)
        m_min, parts = min_partitions(web, sub)
        wr = {'web': wi, 'nodes': web,
              'node_info': [f"{occs[o][0]}[{occs[o][1]}:{occs[o][2]}]" for o in web],
              'base_violations': len(viol),
              'm_min': m_min, 'n_occs': len(web),
              'n_minimal_partitions': len(parts),
              'minimal_partitions': [sorted(sorted(c) for c in p) for p in parts[:50]],
              'pairs': {}}
        nd = sum(1 for v, _ in pv.values() if v == 'DISTINCT')
        ne = sum(1 for v, _ in pv.values() if v == 'SUPPORTED-EQUAL')
        nf = sum(1 for v, _ in pv.values() if v == 'FREE')
        for (a, b), (v, ev) in sorted(pv.items()):
            wr['pairs'][f'{a}-{b}'] = [v, ev]
        print(f'\n  web {wi}: {len(web)} occurrences '
              f'({", ".join(wr["node_info"])})')
        print(f'    base closure violations: {len(viol)}'
              + ('  (per-occurrence-alphabet model holds)' if not viol else
                 '  (!! some window is not monoalphabetic)'))
        print(f'    pairwise: {nd} DISTINCT, {ne} SUPPORTED-EQUAL, {nf} FREE')
        for (a, b), (v, ev) in sorted(pv.items()):
            if v == 'SUPPORTED-EQUAL':
                oa, ob = occs[a], occs[b]
                print(f'      = {oa[0]}[{oa[1]}:{oa[2]}] ~ {ob[0]}[{ob[1]}:{ob[2]}]'
                      f'  ({ev} identity letters)')
        print(f'    m_min = {m_min} of {len(web)} occurrences; '
              f'{len(parts)} minimal partition(s)')
        # cross-checks on the minimal partitions: family & mod-4 composition
        fams, mods = [], []
        for p in parts[:50]:
            for cl in p:
                if len(cl) >= 2:
                    fams.append(sorted({48 if occs[o][0] in FAM48 else 49 for o in cl}))
                    mods.append(sorted({occs[o][1] % 4 for o in cl}))
        wr['merged_class_families'] = fams
        wr['merged_class_mod4'] = mods
        report['webs'].append(wr)
    return report


def main():
    print('== control gates ==')
    gates, controls = [], {}
    for mode in ('quag4', 'fresh', 'transl'):
        r = run_control(mode)
        controls[mode] = r
        if r['n_straddling']:
            print(f'  info: {mode}: {r["n_straddling"]} straddling occurrence(s) '
                  'excluded from gates')
        c = r['checks']
        gates.append(gate(f'{mode}: ground-truth partition accepted', c['accept']))
        gates.append(gate(f'{mode}: m_min == true alphabet count '
                          f'({[w["m_min"] for w in r["webs"]]} vs '
                          f'{[w["m_true"] for w in r["webs"]]})', c['m_min']))
        if mode == 'quag4':
            gates.append(gate('quag4: minimal partition unique == truth', c['exact']))
    allok = all(gates)
    if not allok:
        print('\nGATE FAILURE -> real verdicts NOT trustworthy')

    print('\n== REAL CORPUS ==')
    real = run_real()

    json.dump({'controls': controls, 'real': real},
              open(f'{ROOT}/data/alphacount_report.json', 'w'), indent=1)
    print('\nwrote data/alphacount_report.json')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
