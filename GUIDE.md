# Operator's guide — noiteyes toolchain & attack plan

How to run every tool, how to read what it prints, and the step-by-step plan
for cracking (or decisively characterizing) the Eye Messages cipher. Written
so a fresh session can pick up without re-deriving anything.

Everything runs locally with only the Python standard library. No network, no
API calls, no cost per run. The whole point of the design is that you re-run
freely and trust the output *because the controls pass*, not because a model
told you it worked.

---

## 0. The trust model — read this first

Every result-producing tool cracks **synthetic ciphers of known construction
before it is allowed to report on the real corpus.** These are "control
gates." A gate is a hard assertion: the tool encrypts a message with a
mechanism whose answer we already know, runs the identical analysis on it, and
checks it recovers the known answer (and, for negative controls, that it does
*not* hallucinate structure in mechanisms that lack it).

- If every gate prints `PASS`, the real-corpus numbers in the same run are
  trustworthy and the process exits `0`.
- If any gate prints `FAIL`, the tool prints `real verdicts NOT trustworthy`
  and exits `1`. **Do not believe the real-corpus lines from a failing run.**

This is why you don't need to spend effort (or tokens) double-checking a run:
a green gate line is a machine-checked proof that the method works on a case
with a known answer. When you change a tool, run it; if a gate flips to FAIL,
the change (or an assumption) is wrong — fix it before reading real output.

Determinism: every tool seeds its RNG (`SEED = 1`) so runs are reproducible.
Same code → same numbers.

---

## 1. Tool-by-tool operation

Run all commands from the repo root (`/home/voyd/git/noiteyes`).

### 1.1 `tools/transcribe.py` — pixels → data

```
python3 tools/transcribe.py
```

Rebuilds `data/messages.json` and `data/messages.txt` from the source images
by exact template matching. Prints a verification report: 3,108 eyes located
at zero-tolerance, exactly 5 pupil patterns, the derived digit-significance
order, and the proof that exactly one of the 36 reading orders yields a
contiguous 0–82 alphabet. **You normally never need to run this** — the data
is committed. Run it only if you suspect the transcription or change the
images. If any count in its report differs from the README table, stop:
the dataset changed under you.

### 1.2 `tools/isomorphs.py` — the isomorph catalog

```
python3 tools/isomorphs.py
```

Writes `data/isomorphs.json` (18 maximal isomorph pairs) and prints each with
its length `L` and `evidence_positions` (`ev`). **Read `L` as an upper bound
on the shared phrase and `ev` as its real support** — a window can overrun the
true phrase through letters that occur only once (they carry no repeat
evidence). Any tool that consumes windows must account for this (see the
witness-core repair in §1.4).

### 1.3 `tools/pi_solver.py` — additive-walk solver (historical)

```
python3 tools/pi_solver.py [min_len]
```

The original GF(83) additive-walk ("cut deck") test. **Kept for provenance
only.** Its full-window equation system is contaminated by window overrun, so
its "refuted" verdict — while it happens to be correct — was reached unsafely.
`sigma_web.py` supersedes it with a validated method. Use `pi_solver` only for
its still-good primitives (`find_isomorph_pairs`, `gauss_nullspace`,
`load_msgs`), which the newer tools import.

### 1.4 `tools/sigma_web.py` — σ-algebra battery (mechanism family test)

```
python3 tools/sigma_web.py
```

The third-phase workhorse. Extracts the letter-maps σ between isomorph
occurrences and asks what algebraic family they belong to, then tries to crack
the corpus as a bijective-deck walk over a **witness-core** system (jump-safe
segments that survive window overrun). Writes `data/sigmas.json` and
`data/sigma_report.json`.

Control gates (all must PASS): `strict`, `sparse`, `rowkey` walks must be
cracked end-to-end (consistent system, injective pinned solution, plaintext
IoC ≈ 0.07, recovered deck matching truth 1.0); `affine` and `shuffle` must
yield no exact isomorphs.

How to read the real-corpus block:
- `GF(83) CORES: dim 1` → the jump-safe system is inconsistent → **the whole
  piecewise-translation / group-walk family is refuted** (not just the strict
  walk). This is the headline third-phase result.
- `consensus: UNDERDETERMINED/INCONSISTENT` → no translation-consistent subset
  of segments exists.
- `modulus sweep … all dim 1` → no deck size (prime 59–127) rescues the walk.
- `commutation 6/16` with all successes in the family-48 web → σ maps don't
  commute; they're arbitrary-looking fixed substitutions.
- A subsystem that is "consistent" only by collapsing 50–67 letters onto 2
  values is the **mass-merge pathology** — genuine non-translation structure,
  reported as such, never as a success.

### 1.5 `tools/rowphase.py` — schedule & period extractor (fourth phase)

```
python3 tools/rowphase.py
```

Two independent instruments; writes `data/rowphase_report.json`.

- **Instrument 1 — break phase (row-lock test).** Are the constant-σ stretch
  boundaries (isomorph window ends + closure-violation zones) locked to the
  physical eye rows? Key line is the **nearest-boundary `within1`** z-score
  (are breaks piled *at* a boundary), not just mean distance. Also a
  candidate-grid scan (periods 2/4/13/26) and a residue-mod-k scan.
- **Instrument 2 — keystream period.** Distance-`d` equal-letter recurrence
  (d = 1..8) against an **exact within-message shuffle null**. `d=4` elevated
  = the period-4 cycle.

Gates: `rekey-26` (breaks pile at boundaries when they truly do), `no-schedule`
(they don't when there's no schedule), `shuffle-null` (null calibrated to
|z|<3), `planted-d4` (probe detects a known distance-4 excess).

How to read the real-corpus block (current state):
- Instrument 1: `within1 z≈−0.3` at period 26 and `bracket rate ≈ null` →
  **not row-locked.** The printed structural fact (L=33 isomorph constant
  across boundary 52) is the clinching, synthetic-free evidence.
- Instrument 2: `d=4 z≈4.4` and zone `mod 4 z≈4.3` → **period-4 confirmed**
  from two independent observables. `d=1` (zero doubles) and `d=2` (depressed)
  recur. Single-observable mod-5/mod-7 blips = noise, not signal.
  **CORRECTION (fifteenth phase):** the saved report has zone mod 4 z=0.97
  (mod 5 z=3.64); the "z≈4.3" was a transcription of the d4 figure. And the
  d4 excess itself is a duplication artifact — see §1.12. Period-4 is NOT
  confirmed.

### 1.6 `tools/phase4.py` — Step A: period-4 as an additive key index

```
python3 tools/phase4.py
```

Tests whether the confirmed period-4 cycle is an *additive* keystream index
that splits the cipher into interleaved translation streams. Splits isomorph
pairs by δ = (relative occurrence offset) mod 4 and checks three signatures:
S1 δ-distribution (additive-p4 piles pairs at δ≡0), S2 mod-4 repeat alignment
(additive-p4 forces δ≢0 pairs to ~1.0; chance is ~0.25), S3 stratified solve
(does the δ≡0 subset crack as a translation?). Gates: a planted period-4 walk
must concentrate at δ≡0 and let δ≡0 recover the deck; a plain walk must spread
and crack whole; shuffled text cracks nothing. Writes `data/phase4_report.json`.

Result (all gates pass): **refuted on all three signatures** — real pairs
spread 6/5/5/2 over δ, δ≢0 alignment is at chance, and δ≡0 pins only 3 letters.
The period-4 is real but not an additive key index; it must enter through the
non-commuting substitution structure. Routes to Step B/C. *(The premise
has since dissolved: §1.12 shows period-4 is a duplication artifact.)*

### 1.7 `tools/resync.py` — Step C: resync algebra / effective state size

```
python3 tools/resync.py
```

Measures how often two messages, after diverging over a gap of different
content, snap back to letter-identical output (a "realignment," chance
≈1.7e-6 per run≥3). The realignment rate bounds the cipher's effective
*accumulating* state: `S_eff ≈ 1/rate`. Calibrated on four ciphers of known
state on one shared-template plaintext (gates: position-keyed → S_eff≈1.8,
small-state M=5 → ≈3.7, additive S=83 → ≈54, shuffle autokey → ∞). Writes
`data/resync_report.json`.

Result (gates pass): the corpus shows 4 structural realignments where a
large-state autokey predicts ~0 → **the community deck-shuffle / S₈₃
group-autokey (≈82! state) is refuted**; S_eff≈15 (< 83, upper bound at 4
events) → small effective state; position 0 is off-chain (all 9 differ there
yet re-sync at position 1). Read the ROBUST line (large-state refuted) as
solid; the SUGGESTIVE line (size ~15) as an upper bound.

### 1.8 `tools/tournament.py` — Step B: small-state generator tournament

```
python3 tools/tournament.py
```

Consolidates a single seven-signature fingerprint (`maxL, walk_dim, commute,
same_web, doubles, d4-ratio, resync, IoC`) and scores candidate mechanisms
against the real corpus with the same detectors used throughout. Candidates
are fixed-key encryptors; each gets phrases rejection-sampled to be repeat-
dense under its *own* encryption (a fair shot at isomorphs). Gates validate the
detectors (additive → crackable + commuting; single-rotor → no long isomorphs).
Writes `data/tournament_report.json`.

Result round 1 (gates pass): **no simple small-state generator matches all
seven** — best was the additive walk at 3/7, and the rotor-cycle candidate is
rejected (it produces no isomorphs, because occurrences land in different
rotor phases).

Round 2 added the quagmire family and the community-derived gak41 (see
README): `quagC bumped` (notch pool + collision bump) ties the additive walk
at 3/7 with structural zero doubles (2/7 in the merged harness); and
**`gak41 affine autokey` (C83:C41 group autokey) scores 5/7 — the only
candidate ever to reproduce walk_dim==1** (robust over 20 seeds) — missing
resync rate (a parameter choice) and d4. A new construction gate (bump
variants must give zero doubles) also passes. `quagC + p4 plaintext`
originally suggested the d4 excess could come from period-4 PLAINTEXT
structure, but that was a single run under the pre-fix harness; merged and
seed-swept it clears d4>1.5 in only 5/20 (see README's merge note). Use the tool as the scaffold for candidates —
add a `make_*` encryptor and re-run.

**Three more candidates (eleventh phase)**, targeting that exact tension:
`make_deck_rebuild` ("deck rebuilt, not shifted" — a bank of `M=15`
independently-random tables, switched by a plaintext-notch autokey, applied
UNCHANGED across a whole stretch instead of progressively rotated);
`make_period4_position` (candidate 1 plus an outer period-4 layer keyed by
ABSOLUTE character position); `make_period4_content` (candidate 1 plus a
period-4 layer keyed by a second, faster plaintext notch instead of raw
position). A new diagnostic, `delta_mod4_spectrum`, histograms each raw
isomorph pair's relative offset mod 4 — the same statistic Step D ran on the
real corpus (offsets at every residue: 0:6,1:2,2:5,3:5).

Result: `make_deck_rebuild` is the **first tested mechanism to produce
abundant, offset-independent exact isomorphs at all** (maxL 39, spectrum
spanning residues 0/1/2/3) — the rotor-cycle made zero. But it fails zero-
doubles (a bare table-swap has no collision-avoidance, so plaintext doubles
still produce ciphertext doubles) and period-4 (no elevation; it has no
period-4 component at all). Composing a period-4 layer on top exposes a
sharp, *predicted-and-confirmed* split: `make_period4_position` collapses
the isomorph yield from dozens of pairs to 2, both at delta≡0 mod 4 — a
position-locked period-4 layer destroys almost all cross-offset isomorphs,
exactly as Step D's real spectrum (spread across every residue) says a real
mechanism must NOT do. `make_period4_content` (the same period-4 idea but
keyed to plaintext-letter identity instead of raw position) keeps isomorphs
at every residue (0/1/2/3, matching real) — but its d4-ratio comes back
~1.0 (not elevated): a content-triggered state change is an irregular
renewal process, too loose to produce the real corpus's crisp lag-4 excess.
**This sharpens the binding tension from two-way to three-way**: abundant
offset-independent isomorphs, non-additive/non-commuting structure, and a
crisp period-4 recurrence cannot yet be reproduced simultaneously by any
tested construction — a mechanism satisfying any two of the three, in every
candidate tried so far, loses the third. (Caveat: `walk_dim`/`commute` for
the new table-swap candidates read as suspiciously high — likely a small-n
artifact of this harness's few seeded repeat-phrases landing in the same
one or two states, not evidence the tables themselves commute; the
isomorph-count and delta-spectrum results rest on more direct, less
model-dependent counts and are the ones to trust.)

### 1.8b `tools/wak_sweep.py` / `tools/eye_mural_scan.py` — Step E, non-bounded sweep

```
python3 tools/wak_sweep.py
python3 tools/eye_mural_scan.py
```

Two companion instruments that remove the "bounded" caveat from the original
Step E (§3, `tools/wak_unpack.py`). `wak_sweep.py` byte-scans **every one**
of `data.wak`'s 14,745 entries (not just filename matches) for an English +
Finnish cipher/eye keyword list; gate: a planted marker string must be found
and attributed to the right file. `eye_mural_scan.py` pixel-scans **every**
PNG/BMP/PSD in the archive (9,046 raster files total) for the literal
eye-glyph outline template from `transcribe.py` (11×7, ≥18/20 outline pixels
lit), looking for a GRID of matches (not a lone eye sprite, of which the
game has several unrelated ones). Three gates: the known source sheets must
recover a huge regular grid (proves the detector works); seeded random noise
must not (calibrates the raw threshold); and the catalogued lone-eye Easter
eggs must never register as a grid — this last gate initially FAILED on
`eyespot.png` (a solid-filled icon trivially satisfies an outline-only
template at every interior pixel) until the detector was strengthened to
also require the eye's always-background corner pixels stay clear, which a
solid fill can't satisfy. Both write reports to `data/`.

Result (all gates pass): still **no shipped decoder** — the keyword sweep's
only hits are unrelated engine/asset strings (a Twitch-integration
`decode`, a trailer PSD layer named "glyph", `Runestone` wand items, a boss
`pupil` sprite, the ubiquitous XML attribute `orientation`); the pixel scan
finds **zero** PNG/BMP/PSD in the whole archive containing an eye-outline
grid anywhere, at native resolution, other than confirming the two already-
catalogued non-mural Easter eggs correctly score zero. `strings` on both
`noita.exe` and `noita_dev.exe` (never searched before) turns up nothing
beyond generic engine text (`Message_*` event names, a PNG codec error, the
`ThreeEyesAreWatchingYou` perk). This makes Step E's "authored/externally-
keyed, not derivable from shipped data" conclusion exhaustive rather than
bounded.

### 1.9 `tools/rotorfit.py` — non-abelian-σ test (general rotor?)

```
python3 tools/rotorfit.py
```

Tests whether σ is a general-rotor power: reads T off a unit-offset edge in an
isomorph web, then checks the OTHER edges satisfy σ = T^δ. Gates validate it
both ways (a true rotor → independent edges fit 1.0; a per-message
monoalphabetic cipher, offset-independent σ → ~0.17). Writes
`data/rotorfit_report.json`.

Result (gates pass): the general-rotor hypothesis is **not supported** — the
one real web with a unit edge fits 0/4 (vs 1.0 for a rotor), consistent with
offset-INDEPENDENT σ. Caveat: only 4 discriminating comparisons. Read it as:
σ is a per-occurrence substitution difference (Pⱼ∘Pᵢ⁻¹), not a positional
rotor, and 18 isomorphs under-constrain the Pᵢ.

### 1.10 `tools/alphacount.py` — minimum-alphabet factorization

```
python3 tools/alphacount.py
```

How many distinct alphabets do the σ maps force? Occurrences are partitioned
into shared-alphabet classes; a partition is consistent iff the class
relations, closed under inverse+composition, stay functional, injective and
identity-on-diagonal. Exact branch-and-bound gives m_min; every pair is also
classified DISTINCT / SUPPORTED-EQUAL / FREE from the all-singleton closure.
σ entries come from each pair's LARGEST witness-span block only (overrun
repair — full-trim σ fails the controls). Gates: quag4 recovered exactly
(accept truth + m_min + unique partition), fresh and translate controls show
no false merges. Writes `data/alphacount_report.json`.

Result (gates pass): the real web is FULLY constrained — web 0 needs exactly
3 alphabets, web 1 exactly 5 of 6 with ONE provable reuse
(west-1[30:52] ~ east-1[30:49], in the known resync region). No small
reused pool; per-stretch alphabets are near-fresh.

### 1.11 `tools/exesearch.py` — EXE re-derivation + wak sweeps

```
python3 tools/exesearch.py
```

Extracts the eye data embedded in `noita.exe` (mov-immediate u64 chunks,
base-7 digits −1, '5' = row/message separator; see README for the encoding)
and compares with the transcription; also sweeps every data.wak text file for
eye-related keywords and every wak PNG for the exact eye template. Gates: the
documented 2021 anchor u64 decodes to our east-1 prefix; a planted synthetic
stream is recovered exactly; the extracted stream equals the transcription.
Needs the Noita install (path constant at the top). Writes
`data/exesearch_report.json`.

Result (gates pass): byte-for-byte match, all 9 messages, interleaved order
E1 W1 E2 W2 E3 W3 E4 W4 E5; NO other decodable runs (no west-5 / tenth
message); wak sweeps find only `data/particles/eye.png`. Provenance closed.

### 1.12 `tools/diagscan.py` — exhaustive isomorph scan + deduplicated recurrence

`python3 tools/diagscan.py [shuffled_trials=20]`. Walks every message pair
at every relative offset, taking maximal consistent-bijection windows
(evidence = repeats matched on both sides); this finds sparse isomorphs that
`find_isomorph_pairs`' 12-letter seed misses. Nulls: shuffled corpora and a
reversed message both give 0 windows at evidence ≥3, so every real window
(28 at ≥3) is genuine; 16 are NEW, including a family-B tail isomorph
(east-4[68:102]~west-4[71:105], ev 7) and east-3 linked to that tail and to
west-1. Then union-find over (msg,pos) merges shared-prefix copies and
isomorph-aligned positions, and counts each distance-d recurrence class once
against Poisson(slots/83). Result: d4 26 raw → 11 dedup vs 7.8 (p=0.16),
level with d5/d9/d13 — **period-4 is a duplication artifact**; what survives
is near-repeat depletion (d1=0, d2=5, d3=3). Writes
`data/isomorphs_diag.json` and `data/dedup_recurrence.json`; leaves
`data/isomorphs.json` untouched. The tournament's period-4 test is moved to
`RETIRED_TESTS` accordingly.

### 1.13 `tools/tournament3.py` — tournament round 3

`python3 tools/tournament3.py` (~10 s, deterministic, 12 seeds per config).
Part 1 checks the real σ in `data/sigmas.json`: 18 σ fall into 14
incompatible groups, all 5 compatible links share an occurrence, nothing
recurs across webs or families, and no σ has a fixed point. Part 2 scores
gak41, a K-state machine (K = 8–30, ± double-avoid) and a hybrid register
(fast context N = 2–5, slow part with M = 5–30 states that 6 trigger letters
advance, additively or by permutation) on 7
tests. Period-4 is reported but not scored; it is replaced by d2/d3
depletion (real r2 = 0.40, r3 = 0.72 against the d5–d10 mean). No config
passes all 7. The state machine is dead (its isomorphs are plain repeats);
the hybrid fails walk dim, resync and σ-recurrence; gak41 leads at 4.4/7 but
fails resync (0.0) and d2 depletion (r2 ≈ 1.7). Writes
`data/tournament3_report.json`; imports the harness from `tournament.py`.

### 1.14 `tools/tournament4.py` — tournament round 4 (repairing gak41)

`python3 tools/tournament4.py` (deterministic; 53 configs × 12 seeds, then
the leaders re-scored over 48 seeds as `confirm_48`). Adds three
modifications to gak41, with the same key material per seed: **reset**
(a trigger plaintext letter returns the state to (1, 17)), **avoid w** (bump
a until the output misses the last w outputs; 'state' keeps the bump,
'out' does not), and **short steps K** (additive steps are distinct values
in 1..K, so on untwisted stretches no d-step sum is 0 mod 83 while dK < 83,
which gives soft d-repeat suppression). Results over 48 seeds:
- **short K=27** scores 5.42/7, with r2/r3 = 0.20/0.44 against the real
  0.40/0.72, and doubles pass in every seed.
- **short K=27 + reset q z j x** gives a mean resync of 0.061 against the real
  0.065. Only 15% of seeds fall inside the resync band, though, and walk dim
  and non-commutation both drop.
- **Avoidance** zeroes d2/d3 entirely instead of depleting them.
- **No config passes all 7 in more than 8% of seeds.**

Writes `data/tournament4_report.json`. It imports the harness from
`tournament.py` and `tournament3.py`.

### 1.15 `tools/tournament5.py` — tournament round 5 (lossy updates)

`python3 tools/tournament5.py [n_seeds] [n_procs]` (parallel; 23 configs ×
200 seeds; ~20 s on enter). It tests four ways to vary the short-step gak41:
- state-driven twists (`stw`), where b twists when a lies in an m-set;
- plaintext overwrite of b (`bset`);
- a coarse-step autokey (`coarse C`), updated by letter class;
- fewer twist letters (`ptwist tw=t`).

Two side tests run with each config: unaligned exact repeats ≥4 (real: 0),
and resync under 1–4-letter plaintext edits, matching the real mismatch
islands. Findings:
- State-driven loss destroys the isomorphs (walk dim ≤ 0.16), so any loss
  must be plaintext-driven.
- bset kills non-commutation (0.04).
- Resets predict unaligned repeats; resync band with no unaligned repeat
  occurs in only 9% of seeds.
- Under short edits, lossless gak41 with 2–3 twist letters realigns into the
  band in 12–21% of seeds (resets: 5%). Resync therefore does not require a
  lossy update.
- Coarse-step autokeys fail zero doubles.

Writes `data/tournament5_report.json`.

---

## 2. Current state of knowledge (one-screen summary)

Confirmed (reproduced from primary data, controls passing):
- 83-letter contiguous alphabet, unique reading order; zero adjacent doubles;
  flat frequencies; 18 certain isomorph pairs.
- **Refuted:** monoalphabetic substitution; the entire bijective-deck
  group-walk family (any additive/translation "cut-deck," any group order
  59–127); collision-skip walks (both sparse and dense); **row-locked
  re-keying**; **period-4 as an additive key index** (Step A — the cycle is
  real but does not split the cipher into interleaved translation streams);
  **the large-state deck-shuffle / S₈₃ group-autokey** (Step C — the resyncs
  need a small state, ~82! is impossible).
- **Small effective state:** the accumulating state is order ~10-20 (Step C,
  upper bound), the header (position 0) is off-chain. ~~There is a period-4
  component~~ (retracted: duplication artifact, §1.12).
- **Confirmed structure:** non-commuting fixed substitutions constant over
  18–33-letter stretches, occasionally identity (resyncs); ~~a period-4 keystream
  sub-cycle~~ (retracted, §1.12); near-repeat depletion d1–d3.
- Community consensus (polyalphabetic, non-cyclic-group, deck/shuffle autokey,
  period-4, unsolved) agrees with all of the above; their S₈₃ Group-Autokey
  model is an unproven superset of our surviving family. See README
  "Community cross-check."

- **Consolidated fingerprint** (Step B, all present): exact isomorphs to L=33,
  dim-1 additive solve, non-commuting σ (0.38; but same-stretch σ tend to
  commute), zero doubles, ~~period-4 (d4-ratio 2.6)~~ (duplication artifact, §1.12), resync ~0.065, IoC 0.0128.
- **Alphabet count** (alphacount): the σ web forces exactly 3 + 5 distinct
  alphabets over its two components, with ONE provable reuse (e1~w1 resync
  region) — near-fresh alphabets, no small reused pool.
- **Provenance closed** (exesearch): the eye stream extracted from noita.exe
  equals the transcription byte-for-byte; no extra messages exist in the
  binary; the wak contains nothing eye-related beyond one particle sprite.

Open: the exact surviving mechanism. Step B isolated the crux; the non-abelian-σ
test (`rotorfit`) then ruled out the cleanest resolution (a general progressive
rotor: σ = T^δ fits 0/4 on real vs 1.0 for a true rotor). Current best
characterization: σ is a fixed PER-OCCURRENCE substitution difference
(Pⱼ∘Pᵢ⁻¹), constant within a stretch, changing per stretch with NO positional
(rotor/progressive) organization — and 18 isomorphs under-constrain the Pᵢ. The
ciphertext alone appears to under-determine the schedule; further progress
likely needs an external constraint (a crib on the family-48 web, or the game's
generation code) — and the eleventh-phase game-data sweep below found neither.

The eleventh-phase tournament additions turned that qualitative
characterization into a **three-way quantitative tension**: (1) abundant,
OFFSET-INDEPENDENT exact isomorphs (real spans all 4 residues mod 4 —
0:6,1:2,2:5,3:5), (2) non-additive/non-commuting structure, (3) a crisp
period-4 recurrence excess (d4-ratio 2.6). A plaintext-notched autokey that
REBUILDS (not shifts) its substitution table gets (1)+(2) for the first
time in this project (a prior naive rotor-cycle got neither) but has no
period-4 component at all, so it misses (3) outright and also fails
zero-doubles (no collision-avoidance). Composing a period-4 layer on top
gets (3) only by keying it to absolute POSITION, which then destroys (1)
(isomorphs collapse to delta≡0 mod 4 only, contradicting the real spectrum);
keying it to plaintext CONTENT instead preserves (1) but is too irregular
to reproduce (3)'s crisp lag-4 statistic. No tested construction gets all
three at once. See §1.8 for the full readout. **Leg (3) is now retired:**
§1.12 shows the d4 excess is duplicated material counted many times (11
deduplicated vs 7.8 expected, p=0.16). The tension is two-way: large /
group-structured state (long non-identity isomorphs) vs fast resync
(S_eff ≲ 15, from only 4 realignments).

Since then (parallel second-machine session, merged 2026-09-28): `alphacount`
shows the alphabets are near-fresh per stretch (≥8 distinct across 9
occurrences, ONE provable reuse in the e1~w1 resync region), and `exesearch`
closes provenance — the shipped binary contains exactly our ciphertext and
nothing else, so the answer must come from the ciphertext.

**Current leading candidate: the C83:C41 affine group-autokey (gak41).** From
the community's GAK classification (every other group ruled out — theirs or
ours), it scored 5/7 on the fingerprint (**5/6 now that period-4 is retired —
its only miss is resync**) including walk_dim==1, which nothing
else has ever matched (robust: median 5/7 over 20 seeds). It holds legs (1)
and (2) of the tension above; its misses are resync rate (a free parameter)
and period-4 — leg (3) — which is still unexplained: the earlier
"period-4 lives in the plaintext" reading did not survive the merged
harness (quagC+p4 clears d4>1.5 in only 5/20 seeds, and p4 plaintext breaks
gak41's walk_dim==1). Moot now: period-4 was never real (§1.12).

Rounds 3–4 (§1.13–1.14) re-scored gak41 on 7 tests, with period-4 replaced
by d2/d3 depletion. It fails resync (0 by construction, since its state
update is invertible) and depletion (r2 ≈ 1.7 against the real 0.40).
**Short additive steps (1..27) repair depletion** softly, bringing r2/r3 to
0.20/0.44, and lift the mean to 5.4/7. **Rare plaintext resets** (≈0.5% of
letters) reach the real resync rate on average but not reliably per seed,
and they cost walk dim and non-commutation. The best configs pass all 7
tests in at most 8% of seeds. The open question is now narrow: what lossy
state update yields resync ≈0.065 while keeping a dim-1 walk and
non-commuting σ.

Round 5 (§1.15) answers that question: nothing needs to lose information.
State-driven loss breaks the isomorphs, and plaintext-driven loss breaks
non-commutation or predicts unaligned repeats. Meanwhile the *lossless*
short-step gak41 with 2–3 twist letters reaches the real resync rate by
coincidence, as long as the diverging plaintext is short (1–4-letter edits,
like the real islands). **Leading candidate: lossless gak41, short additive
steps 1..27, 2–3 twist letters (5.5/7).** Resync now constrains the plaintext
edits more than the cipher.

---

## 3. The attack plan (phased, with decision points)

Each step is a self-validating script (same gate discipline). Do them in
order; each decision point routes to the next. "STOP-and-decrypt" means: if a
step yields an affine-gauge pi with plaintext-level IoC > 0.05 that is
internally consistent across all nine messages, we have layer 1 — hand off to
readability/plaintext analysis instead of continuing.

**Step A — Period-4 as an additive key index. DONE (`tools/phase4.py`):
REFUTED.** The cycle does not split the cipher into interleaved translation
streams (δ-distribution spread, δ≢0 repeat alignment at chance, δ≡0 subset not
a translation). The period-4 must be non-additive — it enters through the
substitution structure itself. This eliminates the additive branch and makes
Step B/C the live path; any future period-4 model must be non-additive (an
order-4 key STATE feeding the non-commuting substitution, not an order-4
offset on a state walk).

**Step C — Exploit the resync algebra. DONE (`tools/resync.py`).** The
realignments bound the effective accumulating state: the large-state
deck-shuffle / S₈₃ autokey is refuted (it would essentially never realign),
the state is small (S_eff≈15, upper bound), and the header is off-chain. This
reframes Step B — the target is a SMALL-state generator, not a large shuffle.

**Step B — Small-state generator tournament. DONE (`tools/tournament.py`).**
No simple small-state generator reproduces the full fingerprint; the naive
rotor-cycle is rejected (no isomorphs). The tournament isolated the crux (see
below) and remains the scaffold for testing new candidate `make_*` encryptors.
**Eleventh phase: three more candidates tried, tension sharpened from two-way
to three-way (see §1.8/§2)** — a plaintext-notched table-swap autokey
("deck rebuilt, not shifted") is the first mechanism to get abundant
offset-independent isomorphs, but composing a period-4 layer on top gets the
confirmed d4 statistic only by breaking those isomorphs (position-locked) or
keeps the isomorphs but loses the d4 statistic (content-locked). Resolved
in the fifteenth phase: the d4 statistic was a duplication artifact (§1.12).

**Non-abelian-σ step — general rotor? DONE (`tools/rotorfit.py`): NOT
SUPPORTED.** σ ≠ T^δ (fits 0/4 on real vs 1.0 for a true rotor); σ is a
per-occurrence substitution difference with no positional organization, and 18
isomorphs under-constrain it. The ciphertext under-determines the schedule, so
the next leverage is external, not another schedule diagnostic.

**Step D — Family-48 web focus + crib. DONE: route is structurally closed.**
The shared prefix is stronger than "isomorphic": positions 1-2 are the *same
ciphertext* `[66,5]` in all nine messages, family A {east-1, east-2, west-1}
shares an *identical* 24-letter prefix (pos 1-24), and family B {east-4,
east-5, west-4} an identical 20-letter prefix (pos 1-20). But two facts kill
the crib plan: (1) **no isomorph edge touches the prefix** — every catalogued
occurrence starts at position ≥29, so a prefix crib cannot propagate into the
σ-web (the prefix is an isolated island, and the catalog explicitly excludes
identical same-position windows). (2) The prefix has **no internal
self-consistency gate**: its own repeats sit in different keystream phases, so a
guessed plaintext there can't be validated against anything — it would be pure
speculation, which our discipline forbids. Two byproducts worth keeping:
- **Pure period-4 keystream is decisively refuted.** Isomorph offsets span all
  residues mod 4 (0:6, 1:2, 2:5, 3:5) with non-identity σ; a period-4 keystream
  would force every long isomorph to offset ≡0 mod 4 with σ=identity. So
  period-4 is a keystream *component*, not the whole keystream (consistent with
  the small extra accumulating state, S_eff≈15).
- **Families A and B diverge in content at position 3** (mechanism-agnostic:
  they share plaintext pos 1-2, so any prior-plaintext-determined keystream
  gives identical c[3] iff p[3] identical — but c[3] differs). The nine messages
  are genuinely different texts sharing only a short common opening.

**Step E — Game-install / source verification. DONE, and now EXHAUSTIVE: no
shipped decoder.** Noita *is* installed here; `tools/wak_unpack.py`
(pure-Python, format decoded and validated: 14,745 files) reads `data.wak`
without launching the game. The original pass was a bounded filename search;
the eleventh phase removed that bound (`tools/wak_sweep.py` +
`tools/eye_mural_scan.py`, §1.8b): every one of the 14,745 files' raw bytes
keyword-scanned, every one of the 9,046 raster images (PNG/BMP/PSD)
pixel-scanned for the literal eye-outline template, and both `noita.exe` /
`noita_dev.exe` (never searched before) `strings`-scanned. **No decoder,
plaintext table, or eye-message mural anywhere in the shipped game or its
executables.** The eye-adjacent assets remain the same two unrelated Easter
eggs already catalogued: `caves/eye_0*.png` (grayscale "watching eye" cave
decoration) and `eyespot_a..e`/`book_s_a..e` (the *tripping* books, readable
English "Notes on Grand Alchemy" flavor text). This bounds the project for
good: the shipped game does not contain the answer, so the eye-message
cipher is authored/externally-keyed and the ciphertext really is all we
have to work from. (What's still technically unsearched: the executables'
non-string binary logic itself — e.g. a decoder implemented as code with no
giveaway string — which would need disassembly, a much larger undertaking
not attempted here.)

**Step F — affine-π solver. DONE (`tools/affinefit.py`): INCONCLUSIVE.** The
real web fits the affine-QR model only once the east-1/west-1 excursions are
cut, and then with 2,050+ assignments and nothing pinned beyond the gauge —
too sparse to decide (see README). Original plan, kept for reference: Every solve so
far was translation-only: σ modeled as π(B) = π(A) + δ. The gak41 candidate
(5/7) predicts σ_k = π⁻¹(β_k·π + α_k) with β_k in the order-41 subgroup of
Z83* — never tested. Build it with the usual gates (a synthetic gak41 corpus
must be cracked end-to-end; translate and shuffle controls must behave):
unknowns π (83 values) + (β_k, α_k) per pair; eliminate α by differencing
entries within a pair; search β over the 41 subgroup values per web with the
linear π-system as the consistency oracle. STOP-and-decrypt rule applies: an
affine-consistent π must survive the gauge-invariant IoC test and produce
readable plaintext across all nine messages before belief. A negative result
exhausts the community's GAK classification for pure group-autokeys —
publishable either way.

**Parked / opportunistic:**
- **More ciphertext.** West-5 is absent (now proven absent from the shipped
  binary too — `exesearch` found no extra decodable runs); if a mechanism
  narrows to a small key space, the community's "need more messages" caveat
  applies — note it, don't fabricate around it.
- ~~**Period-4 as plaintext structure.**~~ Dropped — period-4 is a
  duplication artifact (§1.12). Original note: A single quagC+p4 run suggested the
  d4 excess could come from list-like plaintext, but it did not hold up over
  seeds (5/20). Still cheap to check: if Step F yields a π, look at the
  decrypted stream for a 4-periodic delimiter before hunting more key
  structure.

Ground rule for all steps (unchanged): **trust nothing we haven't reproduced
from primary data**, and believe no breakthrough until it survives a
gauge-invariant IoC test and produces internally consistent, readable
plaintext across all nine messages.
