# Tournament round 3 (tools/tournament3.py -> data/tournament3_report.json)

Run: `python3 tools/tournament3.py` (~10 s, deterministic). 12 seeds per config.
d4 is retired: it is reported as `mean_d4_ratio_RETIRED` but not scored. It is
replaced by a d2/d3 depletion test. Real r2 = d2/mean(d5..d10) = 0.40 and
r3 = 0.72; a candidate passes if r2 and r3 are each <= max(0.75, 1.5x real).

## Distinct-sigma discriminator (real data)
- data/sigmas.json holds 18 sigmas in 2 occurrence webs. Two sigmas are called
  compatible when, in either orientation, they share at least 3 agreeing keys
  and have no conflicts. Greedy grouping by that rule leaves 14 groups.
- There are 5 compatible links, and every one shares an occurrence, so each is
  duplication or composition, not independent recurrence. Two of the five come
  from the same pair (e1[40:49]~e1[68:77]), which appears twice in sigmas.json.
- Cross-web and cross-family recurrence: 0. The random-null expectation is also
  0.000, so a single real recurrence would have been significant.
- The mean fixed-point fraction of the real sigmas is 0.0, meaning no isomorph
  is a straight repeat.
- Limit on power: with only 2 webs, a count of 14 groups does not by itself
  violate a K^2 bound for K >= 4. What rules out the K-state machine is the
  mechanism, below.

## Findings
- Small-state machine, K = 8..30, with or without double-avoid: its isomorphs
  are straight repeats (fixed fraction about 1.0), because a random g moves the
  state every letter, and a non-identity sigma needs the state pair to be held
  for 25+ letters. Walk-dim==1 passes 0/12, noncommute passes 0/12, and resync
  runs 0.43-0.57, far too high. Dead.
- Hybrid register, N = 2..5, M = 5..30, additive or permutation slow part: it
  passes iso, doubles and mostly IoC. walk_dim==1 passes about 0/12 (the null
  space is large, 5-45), resync runs 0.1-0.6 (too high), and sigmas are
  26-74% identity. Several configs show non-identity cross-web sigma
  recurrence (real: none). An exploratory (unreported) trigger-rate sweep (2/6/15/30 trigger letters)
  only trades isomorph length against resync; walk_dim never reaches 1.
- gak41 baseline: mean 4.4 of 7 tests. It fails resync (0.0 in every seed) and
  d2 depletion (r2 about 1.7), and gets walk_dim==1 in 7/12 seeds.
- No config passes all 7 tests in any seed.
