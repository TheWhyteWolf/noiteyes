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
non-commuting substitution structure. Routes to Step B/C.

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

Result (gates pass): **no simple small-state generator matches all seven** —
best is the additive walk at 3/7, and the rotor-cycle candidate is rejected
(it produces no isomorphs, because occurrences land in different rotor phases).
The tool prints the binding tension for the next iteration: a mechanism that
makes abundant cross-offset EXACT isomorphs (additive-like) yet is
non-additive, non-commuting, and realigns at ~0.065. Use it as the scaffold
for testing new candidate generators — add a `make_*` encryptor and re-run.

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
  upper bound), the header (position 0) is off-chain, and there is a period-4
  component.
- **Confirmed structure:** non-commuting fixed substitutions constant over
  18–33-letter stretches, occasionally identity (resyncs); a period-4 keystream
  sub-cycle.
- Community consensus (polyalphabetic, non-cyclic-group, deck/shuffle autokey,
  period-4, unsolved) agrees with all of the above; their S₈₃ Group-Autokey
  model is an unproven superset of our surviving family. See README
  "Community cross-check."

- **Consolidated fingerprint** (Step B, all present): exact isomorphs to L=33,
  dim-1 additive solve, non-commuting σ (0.38; but same-stretch σ tend to
  commute), zero doubles, period-4 (d4-ratio 2.6), resync ~0.065, IoC 0.0128.

Open: the exact surviving mechanism. Step B isolated the crux; the non-abelian-σ
test (`rotorfit`) then ruled out the cleanest resolution (a general progressive
rotor: σ = T^δ fits 0/4 on real vs 1.0 for a true rotor). Current best
characterization: σ is a fixed PER-OCCURRENCE substitution difference
(Pⱼ∘Pᵢ⁻¹), constant within a stretch, changing per stretch with NO positional
(rotor/progressive) organization — and 18 isomorphs under-constrain the Pᵢ. The
ciphertext alone appears to under-determine the schedule; further progress
likely needs an external constraint (a crib on the family-48 web, or the game's
generation code).

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

**Step E — Game-install / source verification. DONE: no shipped decoder.**
Noita *is* installed here; `tools/wak_unpack.py` (pure-Python, format decoded
and validated: 14745 files) reads `data.wak` without launching the game. A
thorough bounded search found **no decoder, plaintext table, or eye-message
mural** in the shipped data. The eye-adjacent assets are unrelated Easter eggs:
`caves/eye_0*.png` are grayscale "watching eye" cave-decoration frames, and the
five `eyespot_a..e` + `book_s_a..e` are the *tripping* books — readable English
"Notes on Grand Alchemy" flavor text, not our cipher. This bounds the project:
the shipped game does not contain the answer, so the eye-message cipher is
authored/externally-keyed and the ciphertext really is all we have to work from.
(Caveat: bounded, not an exhaustive read of all 14745 files — a full text
extraction + eye-mural pixel search remains available if wanted.)

**Parked / opportunistic:**
- **More ciphertext.** West-5 is absent; if a mechanism narrows to a small key
  space, the community's "need more messages" caveat applies — note it, don't
  fabricate around it.

Ground rule for all steps (unchanged): **trust nothing we haven't reproduced
from primary data**, and believe no breakthrough until it survives a
gauge-invariant IoC test and produces internally consistent, readable
plaintext across all nine messages.
