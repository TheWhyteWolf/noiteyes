#!/usr/bin/env python3
"""Tournament round 3: small-state machine vs hybrid register vs gak41.

Part 1 -- distinct-sigma discriminator on REAL data (data/sigmas.json):
  a K-state machine c_t = T[s_t][p_t] can only give a constant non-identity
  sigma across an isomorph window when the state PAIR (s_a, s_b) is held, so
  sigma = T[s_b] T[s_a]^-1 and there are at most K^2 of them, which must RECUR
  across unrelated pairs once #pairs is comparable to K^2. Count mutually
  compatible sigma groups (union is a partial bijection with >=MIN_OV shared
  agreeing keys), orientation-aware, and ask whether any recurrence links pairs
  from different webs / families (i.e. not explainable by duplication).

Part 2 -- multi-seed tournament on the tournament.py harness (corpus(),
  signatures()); d4 is REPORTED but RETIRED from scoring (duplication artifact,
  p=0.16 after dedup). Adds near-repeat depletion d1..d3 and a sigma-recurrence
  statistic computed identically on sims and real data.

Writes data/tournament3_report.json. Run from repo root: python3 tools/tournament3.py
"""
import json
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tournament import (P, PA, WEIGHTS, corpus, signatures, make_gak41,  # noqa: E402
                        load_msgs, N9, find_isomorph_pairs, MINLEN,
                        witness_trim, extract_sigmas, build_webs)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEEDS = list(range(1, 13))       # 12 seeds per config
MIN_OV = 3                        # shared agreeing keys to call two sigmas "the same"
FAM_A = {'east-1', 'east-2', 'west-1'}
FAM_B = {'east-4', 'east-5', 'west-4'}


# ------------------------------------------------------------- sigma algebra

def inv(m):
    return {v: k for k, v in m.items()}


def compat(m1, m2):
    """(agree, conflict): agree = shared keys with equal image; conflict = any
    key or value clash (union not a partial bijection)."""
    agree = conf = 0
    for k, v in m1.items():
        if k in m2:
            if m2[k] == v:
                agree += 1
            else:
                conf += 1
    r1, r2 = inv(m1), inv(m2)
    for v, k in r1.items():
        if v in r2 and r2[v] != k:
            conf += 1
    return agree, conf


def best_compat(m1, m2):
    """Orientation-aware: sigma(A->B) matches sigma(C->D) or its inverse."""
    a1, c1 = compat(m1, m2)
    a2, c2 = compat(m1, inv(m2))
    return max([(a1, c1, '+'), (a2, c2, '-')], key=lambda t: (t[1] == 0, t[0]))


def fixed_frac(m):
    return sum(1 for k, v in m.items() if k == v) / max(1, len(m))


def sigma_groups(sigs, occ_pairs, web_of):
    """Greedy grouping of sigmas into mutually compatible classes; report
    recurrences (compatible, >=MIN_OV agreeing keys) and classify each as
    'shared-occurrence', 'same-web' or 'cross-web'."""
    n = len(sigs)
    links = []
    for i in range(n):
        for j in range(i + 1, n):
            a, c, o = best_compat(sigs[i], sigs[j])
            if c == 0 and a >= MIN_OV:
                oi, oj = occ_pairs[i], occ_pairs[j]
                share = any(overlap(x, y) for x in oi for y in oj)
                wi, wj = web_of.get(i), web_of.get(j)
                kind = ('shared-occurrence' if share else
                        'same-web' if wi is not None and wi == wj else 'cross-web')
                nonid = fixed_frac(sigs[i]) < 0.5 and fixed_frac(sigs[j]) < 0.5
                links.append({'i': i, 'j': j, 'agree': a, 'orient': o, 'kind': kind,
                              'nonid': nonid})
    # greedy clique cover on the compatibility relation (any-overlap, no conflict)
    order = sorted(range(n), key=lambda k: -len(sigs[k]))
    groups = []
    for k in order:
        placed = False
        for g in groups:
            if all(best_compat(sigs[k], sigs[x])[1] == 0 and
                   best_compat(sigs[k], sigs[x])[0] >= MIN_OV for x in g):
                g.append(k); placed = True; break
        if not placed:
            groups.append([k])
    return {'n_sigmas': n, 'n_groups': len(groups),
            'groups_multi': [g for g in groups if len(g) > 1],
            'n_links': len(links),
            'nonid_cross_web': sum(1 for l in links if l['nonid'] and l['kind'] == 'cross-web'),
            'links_by_kind': dict(Counter(l['kind'] for l in links)),
            'links': links,
            'mean_fixed_frac': round(sum(fixed_frac(s) for s in sigs) / max(1, n), 3)}


def overlap(x, y):
    return x[0] == y[0] and x[1] < y[2] and y[1] < x[2]


def null_links(sigs, rng, reps=200):
    """Null: same domain sizes, fresh random bijections per sigma -> how many
    >=MIN_OV compatible links arise by chance."""
    tot = 0
    for _ in range(reps):
        fake = []
        for s in sigs:
            perm = list(range(P)); rng.shuffle(perm)
            fake.append({k: perm[k] for k in s})
        tot += sum(1 for i in range(len(fake)) for j in range(i + 1, len(fake))
                   if (lambda r: r[1] == 0 and r[0] >= MIN_OV)(best_compat(fake[i], fake[j])))
    return tot / reps


def sigma_stat(S, names):
    """Same pipeline as signatures(): raw isomorphs -> witness trim -> sigmas."""
    raw = [p for p in find_isomorph_pairs(S, names=names) if p[2] - p[1] >= MINLEN]
    if not raw:
        return None
    joint = witness_trim(S, raw)
    trimmed = [e['trim'] for e in joint]
    sigs = extract_sigmas(S, trimmed)
    occ_pairs = [((n1, l1, r1), (n2, l2, r2)) for n1, l1, r1, n2, l2, r2 in trimmed]
    occs, pot, webs, gv, edges = build_webs(trimmed)
    web_of = {}
    for a, b, off, k in edges:
        for wi, w in enumerate(webs):
            if a in w or b in w:
                web_of[k] = wi
    return sigma_groups(sigs, occ_pairs, web_of)


def real_sigma_analysis(rng):
    raw = json.load(open(os.path.join(ROOT, 'data', 'sigmas.json')))
    sigs = [{int(k): int(v) for k, v in e['sigma'].items()} for e in raw]
    occ_pairs = []
    for e in raw:
        a, s1, e1, b, s2, e2 = e['pair']
        occ_pairs.append(((a, s1, e1), (b, s2, e2)))
    # web = connected component over overlapping occurrences
    parent = list(range(len(raw)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i in range(len(raw)):
        for j in range(i + 1, len(raw)):
            if any(overlap(x, y) for x in occ_pairs[i] for y in occ_pairs[j]):
                parent[find(i)] = find(j)
    web_of = {i: find(i) for i in range(len(raw))}
    res = sigma_groups(sigs, occ_pairs, web_of)
    res['pairs'] = [e['pair'] for e in raw]
    res['sizes'] = [len(s) for s in sigs]
    res['n_webs'] = len(set(web_of.values()))
    res['null_links_random'] = null_links(sigs, rng)
    fam = lambda m: 'A' if m in FAM_A else 'B' if m in FAM_B else m
    res['cross_family_links'] = sum(
        1 for l in res['links']
        if {fam(occ_pairs[l['i']][0][0]), fam(occ_pairs[l['i']][1][0])} !=
           {fam(occ_pairs[l['j']][0][0]), fam(occ_pairs[l['j']][1][0])})
    # pairwise matrix of best (agree, conflict) for the record
    res['pairwise'] = [[list(best_compat(sigs[i], sigs[j])[:2]) for j in range(len(sigs))]
                       for i in range(len(sigs))]
    return res


# ---------------------------------------------------------- near-repeat depletion

def depletion(S, names):
    """d_k = #(same letter at distance k); ratio to the mean of d5..d10."""
    def dc(d):
        return sum(1 for n in names for i in range(d, len(S[n])) if S[n][i] == S[n][i - d])
    ref = sum(dc(d) for d in range(5, 11)) / 6 or 1
    return {f'd{k}': dc(k) for k in (1, 2, 3, 4)} | {
        'ref_d5_10': round(ref, 2),
        **{f'r{k}': round(dc(k) / ref, 3) for k in (1, 2, 3, 4)}}


# ------------------------------------------------------------------ candidates

def make_ssm(K, avoid, seed):
    """Candidate 1: K-state machine. c = T[s][p]; s' = g[s][p]. avoid: if the
    output would equal the previous cipher letter, step the state s->s+1 (mod K)
    until it doesn't (deterministic in (s, p, prev))."""
    rng = random.Random(seed * 7919 + K * 31 + avoid)
    T = []
    for _ in range(K):
        t = list(range(P)); rng.shuffle(t); T.append(t)
    idx = {c: i for i, c in enumerate(PA)}
    g = [[rng.randrange(K) for _ in range(P)] for _ in range(K)]

    def enc(s):
        st = rng_state0
        out = [T[st][idx[s[0]]]]
        st = g[st][idx[s[0]]]
        for ch in s[1:]:
            p = idx[ch]
            c = T[st][p]
            if avoid:
                k = 0
                while c == out[-1] and k < K:
                    st = (st + 1) % K; c = T[st][p]; k += 1
            out.append(c)
            st = g[st][p]
        return out
    rng_state0 = rng.randrange(K)
    return enc


def make_hybrid(N, M, mode, seed, trig=6):
    """Candidate 2: fast register x shifting slow component.
    fast: context ctx = hash of previous N-1 plaintext letters picks one of R
          bijections F[ctx] (so equal N-grams -> equal fast output; resync after
          N letters of re-agreement).
    slow: state m in 0..M-1, advanced only by `trig` trigger letters
          ('add': m += a[p] mod M;  'perm': m = pi_p(m), pi_p random perms of M
          -> non-abelian walk). Output c = B[m][F[ctx][p]] with B a bank of M
          random bijections (non-commuting sigmas B[m']B[m]^-1).
    doubles: if c == prev, use B[m][F[ctx][p]+1] (deterministic, sigma-covariant)."""
    rng = random.Random(seed * 104729 + N * 1009 + M * 17 + (mode == 'perm'))
    idx = {c: i for i, c in enumerate(PA)}
    R = 64
    F = []
    for _ in range(R):
        f = list(range(P)); rng.shuffle(f); F.append(f)
    B = []
    for _ in range(M):
        b = list(range(P)); rng.shuffle(b); B.append(b)
    hv = [rng.randrange(1 << 30) for _ in range(P)]
    triggers = set(rng.sample(range(P), trig))
    a = {p: rng.randrange(1, M) for p in triggers}
    pi = {}
    for p in triggers:
        q = list(range(M)); rng.shuffle(q); pi[p] = q
    m0 = rng.randrange(M)

    def enc(s):
        ps = [idx[c] for c in s]
        m = m0
        out = [B[m][F[0][ps[0]]]]
        for t in range(1, len(ps)):
            p = ps[t]
            h = 0
            for q in ps[max(0, t - N + 1):t]:
                h = (h * 1000003 + hv[q]) & 0xFFFFFFFF
            f = F[h % R][p]
            c = B[m][f]
            if c == out[-1]:
                c = B[m][(f + 1) % P]
            out.append(c)
            if p in triggers:
                m = (m + a[p]) % M if mode == 'add' else pi[p][m]
        return out
    return enc


# --------------------------------------------------------------------- scoring

def tests(sig, dep, real_dep):
    """Scored tests (d4 RETIRED). Depletion: r2, r3 no higher than 1.5x real
    (real is depleted) and doubles == 0."""
    return {
        'iso L>=25': sig['maxL'] >= 25,
        'walk dim==1': sig['walk_dim'] == 1,
        'noncommute<.6': sig['commute'] is not None and sig['commute'] < 0.6,
        'zero doubles': sig['doubles'] == 0,
        'd2,d3 depleted': dep['r2'] <= max(0.75, 1.5 * real_dep['r2']) and
                          dep['r3'] <= max(0.75, 1.5 * real_dep['r3']),
        'resync .03-.12': 0.03 <= sig['resync'] <= 0.12,
        'IoC .011-.015': 0.011 <= sig['ioc'] <= 0.015,
    }


def run_config(label, maker):
    rows = []
    for sd in SEEDS:
        enc = maker(sd)
        rng = random.Random(1000 + sd)
        S, names = corpus(enc, rng)
        sig = signatures(S, names, rng)
        dep = depletion(S, names)
        ss = sigma_stat(S, names)
        rows.append({'sig': sig, 'dep': dep,
                     'sig_groups': None if ss is None else
                     {k: ss[k] for k in ('n_sigmas', 'n_groups', 'n_links',
                                         'links_by_kind', 'mean_fixed_frac',
                                         'nonid_cross_web')}})
    return rows


def summarize(rows, real_dep):
    sc = [tests(r['sig'], r['dep'], real_dep) for r in rows]
    keys = list(sc[0])
    pass_rate = {k: round(sum(s[k] for s in sc) / len(sc), 2) for k in keys}
    all_pass = round(sum(all(s.values()) for s in sc) / len(sc), 2)
    mean = lambda xs: round(sum(xs) / len(xs), 3) if xs else None
    sg = [r['sig_groups'] for r in rows if r['sig_groups']]
    return {'pass_rate': pass_rate, 'mean_tests_passed': mean([sum(s.values()) for s in sc]),
            'all7_rate': all_pass,
            'mean_maxL': mean([r['sig']['maxL'] for r in rows]),
            'mean_commute': mean([r['sig']['commute'] for r in rows if r['sig']['commute'] is not None]),
            'mean_resync': mean([r['sig']['resync'] for r in rows]),
            'mean_d4_ratio_RETIRED': mean([r['sig']['d4_ratio'] for r in rows]),
            'mean_r2': mean([r['dep']['r2'] for r in rows]),
            'mean_r3': mean([r['dep']['r3'] for r in rows]),
            'mean_fixed_frac': mean([g['mean_fixed_frac'] for g in sg]),
            'mean_n_sigmas': mean([g['n_sigmas'] for g in sg]),
            'mean_n_groups': mean([g['n_groups'] for g in sg]),
            'mean_cross_web_links': mean([g['links_by_kind'].get('cross-web', 0) for g in sg]),
            'mean_nonid_cross_web': mean([g['nonid_cross_web'] for g in sg]),
            'frac_seeds_nonid_recur': mean([1.0 if g['nonid_cross_web'] else 0.0 for g in sg])}


def main():
    rng = random.Random(7)
    report = {'seeds': SEEDS, 'min_overlap': MIN_OV}

    print('== Part 1: distinct-sigma discriminator (real data/sigmas.json) ==')
    rs = real_sigma_analysis(rng)
    report['real_sigma'] = rs
    print(f"  {rs['n_sigmas']} sigmas (sizes {rs['sizes']}), {rs['n_webs']} webs")
    print(f"  greedy compatible groups: {rs['n_groups']}  (multi-member: {rs['groups_multi']})")
    print(f"  compatible links (>= {MIN_OV} agreeing keys, 0 conflicts): {rs['n_links']} "
          f"{rs['links_by_kind']}; cross-family {rs['cross_family_links']}; "
          f"random-null expectation {rs['null_links_random']:.3f}")
    for l in rs['links']:
        print(f"    {rs['pairs'][l['i']]} ~ {rs['pairs'][l['j']]}  agree {l['agree']} "
              f"{l['orient']} {l['kind']}")
    print(f"  mean fixed-point fraction of sigmas: {rs['mean_fixed_frac']}")

    S = load_msgs()
    real_sig = signatures(S, N9, random.Random(1))
    real_dep = depletion(S, N9)
    real_ss = sigma_stat(S, N9)
    report['REAL'] = {'sig': real_sig, 'dep': real_dep,
                      'sigma_stat_pipeline': {k: real_ss[k] for k in
                                              ('n_sigmas', 'n_groups', 'n_links',
                                               'links_by_kind', 'mean_fixed_frac',
                                               'nonid_cross_web')},
                      'score': tests(real_sig, real_dep, real_dep)}
    print(f'\n  REAL: {real_sig}\n  REAL depletion: {real_dep}')
    print(f"  REAL sigma stat via harness pipeline: {report['REAL']['sigma_stat_pipeline']}")

    configs = {}
    configs['gak41'] = lambda sd: make_gak41(sd)
    for K in (8, 12, 15, 20, 30):
        for av in (0, 1):
            configs[f'ssm K={K}{" avoid" if av else ""}'] = \
                (lambda K, av: lambda sd: make_ssm(K, av, sd))(K, av)
    for mode in ('add', 'perm'):
        for N in (2, 3, 4, 5):
            for M in (5, 10, 15, 20, 30):
                configs[f'hybrid {mode} N={N} M={M}'] = \
                    (lambda N, M, mode: lambda sd: make_hybrid(N, M, mode, sd))(N, M, mode)

    print(f'\n== Part 2: tournament ({len(SEEDS)} seeds/config; d4 retired) ==')
    res = {}
    for label, mk in configs.items():
        rows = run_config(label, mk)
        sm = summarize(rows, real_dep)
        res[label] = sm
        pr = sm['pass_rate']
        print(f"  {label:26s} " + ' '.join(f'{v:.2f}' for v in pr.values()) +
              f"  | all7 {sm['all7_rate']:.2f} mean {sm['mean_tests_passed']}"
              f" maxL {sm['mean_maxL']} rs {sm['mean_resync']} r2 {sm['mean_r2']}"
              f" r3 {sm['mean_r3']} fix {sm['mean_fixed_frac']} xweb {sm['mean_cross_web_links']}"
              f" nonid-xweb {sm['mean_nonid_cross_web']} ({sm['frac_seeds_nonid_recur']})")
    print('  columns: ' + ' | '.join(next(iter(res.values()))['pass_rate']))
    report['candidates'] = res
    out = os.path.join(ROOT, 'data', 'tournament3_report.json')
    with open(out, 'w') as f:
        json.dump(report, f, indent=1)
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
