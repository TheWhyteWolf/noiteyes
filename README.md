# noiteyes — independent analysis of Noita's Eye Messages

Workspace for analyzing the unsolved Eye Messages cipher from the game Noita
(9 messages of eye glyphs hidden in the East/West parallel worlds, added 2021,
unsolved as of July 2026).

Working principle: **trust nothing we haven't reproduced from primary data.**
The ciphertext here is transcribed from pixel-art images by exact template
matching, and every structural claim made by the research community has been
independently re-derived from that transcription before being relied on.

## Source images (top directory)

- `1x-eyes-east.webp`, `1x-eyes-west.webp` — native-resolution (1 game pixel =
  1 image pixel) sheets of all 9 messages.
- `Eyeorientation.png` — the community's value key (pupil position -> 0..4),
  4x nearest-neighbor upscale, decoded programmatically.
- `2560px-Eye_message_alignments_and_gaps.jpg` — community "Alignments and
  Gaps" chart (Lymm), used ONLY as a cross-check, never as input.
- `Eye Messages - The Noita Wiki.html` — saved wiki page (2026-05-13 revision).

## Pipeline

`tools/transcribe.py` — images -> `data/messages.json` + `data/messages.txt`.
Locates all 3,108 eyes by exact 20-pixel outline match (zero tolerance
failures), classifies pupils (exactly 5 patterns occur, no noise), derives the
value mapping from the key image, groups row-pairs into down/up triangles, and
proves by exhaustion that exactly 1 of the 36 digit-significance orders yields
trigram values covering 0..82 contiguously:

- down-triangle (top,bottom,top): msd=top-left, mid=top-right, lsd=bottom
- up-triangle (bottom,top,bottom): msd=bottom-right, mid=bottom-left, lsd=top

(both clockwise walks). Letter = 25·msd + 5·mid + lsd ∈ [0,82].
Sequence order: row-pairs top→bottom, trigrams left→right.
ASCII rendering (community convention): `chr(32 + letter)`.

## Independently verified facts (all reproduced from raw pixels)

| Fact | Result |
|---|---|
| Corpus size | 9 messages, 3,108 eyes, 1,036 letters (99/118/137/119/114 east, 103/102/124/120 west) |
| Alphabet | exactly 83 distinct letters, contiguous 0–82; unique among all 36 reading orders |
| Row structure | ≤39 eyes/row; bottom two rows shorter, lengths differ ≤1 |
| First letters | all 9 distinct: 50,36,63,27,33 (east), 80,76,34,77 (west) |
| Second letter | identical (66) in all 9 messages |
| Adjacent doubles | **0** observed vs ~12.4 expected — structural, not chance (p≈4e-6) |
| Distance-2 repeats | depressed: 5 vs ~12.3 expected (0.41×) |
| Distance-4 repeats | elevated: 26 vs ~12.0 expected (2.16×) |
| Repeated 3-grams | all 184 repeat-pairs are position-aligned (shared sections); zero elsewhere |
| Shared sections | all 9 share positions 1–2; {e1,e2,w1} share 1–24; {e3,e4,e5,w2,w3,w4} share 1–5; {e4,e5,w4} share 1–20; e3 shares 1–9 with them |
| Re-alignment | e1~w1 diverge at 25–28, re-align 29–32 and 36–48 (equal-length substitutions) |
| Cross-check | matches Lymm's chart character-for-character (east-3 and west-4 verified fully; all 9 at overview) |
| Letter frequency | min 3 / max 26 / mean 12.5, std 4.76 (uniform would be ~3.5); IoC 0.01284 vs 0.01205 uniform |

## Interpretation constraints (for any candidate cipher)

1. Zero adjacent doubles across 1,027 pairs is structural: the mechanism
   cannot emit the same letter twice in a row.
2. Same plaintext at same position must encrypt identically across messages
   (shared sections + re-alignment after equal-length divergent chunks) —
   yet the first letter differs in all 9 without desynchronizing position 1+.
3. Whatever state the cipher keeps, it produces flat-ish frequencies, no
   periodicity, isomorphs, depressed distance-2 and elevated distance-4
   letter recurrence.

## Cipher-analysis findings (second phase)

1. **Isomorph catalog** (`tools/isomorphs.py` → `data/isomorphs.json`):
   18 maximal isomorph pairs. Group A is a mutually-isomorphic web — the same
   phrase appears TWICE in each of east-1, east-2, west-1 (six segments,
   pairwise aligned lengths 18–28). Group B: east-5[33:66] ~ west-4[34:67]
   (L=33), east-4 ~ east-5 and east-4 ~ west-4 overlapping the same region.
   Individual significance 10^-10..10^-15 — these are certain, though maximal
   extents can overrun the true phrase through singleton letters (the JSON
   records `evidence_positions` for this reason).
2. **Strict additive-walk ("cut deck") model REFUTED** (`tools/pi_solver.py`):
   the model c_i = D[x_i], x_i = x_{i-1} + v(p_i) mod 83 with bijective D
   turns each isomorph pair into linear constraints π(B_j) = π(A_j) + δ_pair
   over GF(83), π = D⁻¹. The resulting system is solvable only with a
   non-injective π (~47 forced letter merges), even restricted to trimmed
   pattern-verified cores. Since each pair is individually far beyond chance,
   the pairs are real and the strict model is wrong.
3. **Resync structure**: east-1 ~ west-1 share [1,25), diverge, re-align at
   [29,33) and [37,50); east-4 ~ east-5 share [1,21) then re-align three
   times after divergence gaps of 4, 1, and 6 letters. Any surviving
   mechanism must return to identical output after unequal-content chunks of
   equal length.
4. **Position facts**: position 0 distinct in all 9 messages; positions 1–2
   identical in all 9 ("b%"); the two message families split at position 3 on
   ADJACENT letter values (48 vs 49). Consistent with numbered-list
   plaintexts whose header letter sits outside the chained state.
5. **Leading hypothesis — additive walk + collision-skip**: "if the output
   would repeat the previous letter, step again." Structurally explains the
   zero doubles, predicts near-isomorphs ("slightly differing composition",
   as the community reports), and a few skip-corrupted equations propagating
   through exact Gaussian elimination would produce exactly the mass-merge
   failure in (2). Next-phase experiments (σ-composition diagnostics,
   skip-tolerant refitting, mechanism-enumeration harness) are specified in
   `PLAN.md`. **REFUTED in the third phase — see below.**

A same-day, less rigorous cross-check (`tools/sigma_diagnostics.py`,
predecessor to `sigma_web.py` below, kept for provenance only) reached a
compatible but weaker result without the witness-core repair: composing
the 18 isomorph pairs' letter-maps sigma_AB pairwise agrees in 100+ of
106 tested triples, the 6 disagreements all traced to the maximal-window-
overrun artifact. Composition-holds is a necessarily weaker test than the
third phase's — true for any history-independent state-difference
mechanism, abelian or not — so it does not by itself distinguish the
group-walk family the third phase goes on to refute.

## Sigma-web diagnostics (third phase, `tools/sigma_web.py`)

Model-free study of the letter-mappings σ between isomorph occurrences,
with the instrument validated on synthetic corpora of KNOWN mechanism before
any real-data verdict is read. Five control gates must pass on every run:
strict walk, sparse-jump walk, and row-rekeyed walk must each be cracked
end-to-end (consistent system, injective pinned solution, plaintext-level
difference IoC ≈ 0.07, recovered deck matching ground truth 1.0); affine
walks and shuffled text must yield no exact isomorphs. All five pass.

Methodological repairs the controls forced (both matter for prior work):

- **Overrun contamination**: maximal windows extend through repeat-free
  stretches beyond the true shared phrase, so full-window equation systems
  go inconsistent even when the walk model is TRUE (the strict control
  reproduces `pi_solver.py`'s "refuted" verdict on synthetic data where the
  model holds). Finding 2's original methodology was therefore unsafe.
- **Witness cores**: a letter repeating inside a window at positions p<q
  forces equal state offsets at p and q (its repeat is mirrored on the other
  side), so no net state jump can sit strictly inside a witness span without
  breaking the window. Witness spans, merged only on strict overlap, are
  jump-safe segments for ANY piecewise-translation mechanism at any jump
  density; a consensus prune (drop the most pairwise-conflicted segment)
  additionally removes spans corrupted by exactly-cancelling jump pairs.
  On the controls this pipeline recovers the true deck perfectly.

Real-corpus verdicts (validated instrument, `data/sigma_report.json`):

1. **The whole bijective-deck group-walk family is refuted**, not just the
   strict walk — finding 2's conclusion survives its methodology repair, much
   strengthened. The jump-safe core system is inconsistent (nullspace dim 1);
   consensus pruning finds no usable translation-consistent subset (20 pinned
   letters with 4 forced merges after dropping 5/19 segments); subsystems
   (cross-message, group-A web) are "consistent" only by collapsing 51–67 of
   ~70 letters onto 2 values — the mass-merge pathology of finding 2, now
   shown to be genuine non-translation structure rather than skip corruption.
2. **No cyclic group order rescues it**: the same core system is inconsistent
   over GF(p) for every prime p in 59–127 (deck sizes other than 83 don't
   help within the walk family).
3. **Collision-skip (finding 5) is dead in both regimes**: the sparse-jump
   control proves the pipeline cracks a rarely-skipping walk (so the real
   corpus would have cracked if it were one), and the dense-jump control
   (signed ±8 at 16%/step) produces only 2 short isomorph pairs — dense
   jumping cannot reproduce the real catalog's 18 pairs with an exact L=33
   window. No jump density explains the corpus.
4. **The σ maps do not commute** (6/16 vs 74/74 on every walk control), have
   no fixed points and no closed cycles — they are arbitrary-looking partial
   injections, constant over stretches of 18–33 letters, occasionally exactly
   identity (the resyncs of finding 3). Every commutation success lies inside
   the family-48 web {east-1, east-2, west-1}; family-49 and mixed pairs
   never commute (small n — a lead, not a conclusion).
5. **Breaks correlate with the physical eye rows**: isomorph window ends and
   closure-violation zones sit closer to row-pair boundaries than chance
   (z ≈ −2.9 and −3.0). Note the row-rekey control with row-ALIGNED phrase
   copies shows z ≈ 0 (aligned jumps cancel pairwise), so the real signal
   points at relations breaking where row structure crosses the pair
   MISALIGNED — consistent with per-row keying plus occurrence offsets that
   are not row-multiples.

Surviving mechanism family: per-stretch FIXED substitutions between
occurrences (exact isomorphy over 20–33 letters) that change at row-
correlated boundaries, are not translations of any group of order ≤127,
mostly do not commute, and occasionally return to identity (resyncs).
Natural candidates: a schedule of mixed substitution alphabets keyed by
row (progressive/quasi-periodic, family-dependent), or a deck that is
REBUILT (not shifted) at row-ish boundaries. Next phase: fit alphabet-
schedule models directly; split commutation/consistency by row phase;
test whether σ relations repeat when row-schedule offsets coincide.

## Row-phase & period cross-check (fourth phase, `tools/rowphase.py`)

Two independent instruments, each validated on synthetic corpora of KNOWN
schedule before any real verdict is read (gates: a re-key-every-26 walk must
show its breaks piling at boundaries; a no-schedule walk must not; an exact
within-message shuffle must calibrate to |z|<3; a planted distance-4 excess
must be flagged). All gates pass.

- **Row-locked re-keying is refuted.** The third phase's "breaks correlate
  with rows" (mean distance z≈−2.9) does NOT survive as a *schedule*: real
  break events are not piled at row-pair boundaries (within-±1 z≈−0.3),
  closure zones do not bracket a row crossing above chance (0.36 vs 0.38),
  and — decisively, needing no synthetic — the L=33 isomorph
  east-5[33:66]~west-4[34:67] is exactly constant across row-pair boundary 52
  in *both* occurrences. A full substitution rebuild at row boundaries is
  therefore impossible; the row-distance effect is a weak, diffuse residue,
  not a re-key point. (The natural candidate "deck rebuilt at row-ish
  boundaries" from the third phase is thus downgraded.)
- **[RETRACTED 2026-09-28 — see "Period-4 is a duplication artifact"
  below.** The d4 excess vanishes once duplicated material is counted once,
  and the "mod 4 z≈4.3" zone figure was a transcription error: the saved
  `data/rowphase_report.json` has zone mod 4 z=0.97 (mod 5 z=3.64). Original
  text kept for the record:] A period-4 keystream cycle is confirmed, from two independent
  observables. Distance-4 equal-letter recurrence is elevated (obs 26 vs
  11.6, z≈4.4) against an exact within-message shuffle null, and closure-zone
  break positions concentrate at a fixed residue mod 4 (z≈4.3). Zero doubles
  (distance-1) and depressed distance-2 recur as before. Because the corpus
  also contains L=33 exact isomorphs (long constant-σ runs), this period-4 is
  a keystream *sub-cycle*, not the alphabet-change period. Isolated mod-5/mod-7
  peaks appear in only one observable each and are treated as multiple-
  comparison noise, not signal.

Net redirection: the surviving mechanism is a non-translation fixed-
substitution structure (third phase) carrying a period-4 keystream component
(fourth phase), with no row-locked schedule. Next candidates: a keystream
with an order-4 positional cycle (e.g. a 4-state sub-key or a length-4 rotor
stage) composed with the non-commuting substitution; fit the period-4 phase
directly and test whether σ segments sharing a mod-4 phase are related.

## Step A — period-4 as a key index (fifth phase, `tools/phase4.py`)

Direct test of the plan's Step A: is the confirmed period-4 cycle an ADDITIVE
key index that reduces the cipher to interleaved translation streams? The only
period-4 form compatible with L=33 exact isomorphs is a single deck with a
period-4 additive state offset `k[i mod 4]`; under it, isomorph pairs at
relative offset δ split sharply — δ≡0 mod 4 pairs are pure translations, δ≢0
pairs can only exist with fully mod-4-aligned internal repeats. Three
predictions, each validated on synthetic corpora of known construction (a
planted period-4 walk must concentrate its pairs at δ≡0 and let the δ≡0 subset
recover the deck; a plain walk must spread and crack whole; shuffled text must
crack nothing). All gates pass.

**Refuted, on three independent signatures.** The additive period-4 key index
does not hold:

- **δ-distribution:** real pairs spread over δ mod 4 as 6/5/5/2 (fraction at
  δ≡0 = 0.33), exactly like a plain walk; the planted additive-period-4
  control concentrates 9/9 at δ≡0. Real shows no δ≡0 suppression.
- **Repeat alignment:** real δ≢0 pairs have internal repeats mod-4-aligned at
  chance (0.28/0.26/0.21 for δ = 1/2/3), whereas the additive model *requires*
  ~1.0. The L=33 pair (δ=1) sits at 0.33 — an additive period-4 offset cannot
  produce it.
- **Stratified solve:** the pooled core system is dim-1 (as before) and the
  δ≡0 subset pins only 3 letters (the control's genuine δ≡0 translation web
  pins 74) — δ≡0 is not a hidden consistent translation.

So the period-4 statistic is real (fourth phase) but is **not an additive /
translation keystream index**: it does not decompose the cipher into four
simpler additive streams. Combined with the third phase (whole translation
family refuted) this closes the additive-period-4 branch. The period-4 must
enter through the non-commuting substitution structure itself (a non-additive
order-4 key state), which is Step B/C territory — enumerate small non-additive
period-4 generators and exploit the resync equations, rather than any further
additive fit.

## Step C — resync algebra bounds the state (sixth phase, `tools/resync.py`)

The messages repeatedly RE-SYNCHRONISE: two messages are letter-identical on a
shared prefix, diverge over a gap of different content, then become identical
again for a run whose chance probability is ~(1/83)^run (negligible at run≥3).
Reproduced catalog: the two families share a prefix after a distinct,
off-chain header — {east-1,east-2,west-1} share [1,25), {east-4,east-5,west-4}
share [1,21) — and there are 4 structural mid-realignments (east-1~west-1:
gap4→run13@37 and gap4→run4@29; east-4~east-5: two run-3 realignments).

What a realignment costs bounds the cipher's effective accumulating state: for
a position-keyed or tiny-state cipher realignment is automatic wherever
plaintext re-agrees; for an accumulating autokey of state size S it needs a
state collision (rate ~1/S). The instrument measures the realignment rate and
is calibrated on ciphers of KNOWN state, all on one shared-template plaintext
family (gates pass): position-keyed S_eff≈1.8, small-state M=5 S_eff≈3.7,
additive S=83 S_eff≈54, shuffle autokey S_eff→∞ (never realigns).

- **The large-state autokey is refuted (robust).** The community's leading
  model — a per-letter deck shuffle / S₈₃ group-autokey with ~82! hidden
  states — predicts essentially zero realignments; the corpus shows four
  structural ones. A large accumulating state cannot re-collide to a
  13-letter identical run after four different letters.
- **The effective state is small (suggestive).** Real S_eff≈15 (< 83); the
  corpus realigns more often than an S=83 additive walk would (expected 1.15,
  observed 4, p≈0.03). With only four events, treat 15 as an upper bound, not
  a point estimate — but it is firmly small and consistent with the period-4
  scale.
- **The header sits outside the state chain**, now shown directly: all nine
  differ at position 0 yet re-sync at position 1 (an accumulator would carry
  the header difference forever and never share the [1,25) / [1,21) prefixes).

Net: the cipher is NOT a large-hidden-state autokey; it has a small effective
state (order ~10-20, with a period-4 component) and an off-chain header. This
reframes the mechanism hunt toward small-state generators (Step B) and away
from the S₈₃ deck-shuffle framing.

## Step B — small-state generator tournament (seventh phase, `tools/tournament.py`)

With the state known small (Step C), non-additive and period-4, the remaining
signatures pin the *class*: L=33 EXACT isomorphs force σ constant over 33
positions, which only a progressive key E_i = Tⁱ·B produces (σ = T^δ); the
σ-maps commute within a stretch but not across (verified: same-web 6/14 vs
cross-web 0/2), so T must change between stretches; the dim-1 additive solve
forces T to be a general (non-cyclic) permutation; period-4 suggests a small
cycle of rotors. The derived candidate is a **cycle of 4 general-permutation
rotors, progressive within each stretch, re-keyed on a plaintext notch.**

The tool consolidates a single real fingerprint and scores candidates against
it with the same detectors used throughout (gates validate the detectors: the
additive baseline must be crackable and commuting, a single-rotor Enigma must
make no long isomorphs — both hold). Real fingerprint (all seven present):
`maxL 33, walk_dim 1, commute 0.38, doubles 0, d4-ratio 2.6, resync 0.065,
IoC 0.0128`.

**Result — no simple small-state generator reproduces the full profile, and
the naive rotor-cycle is rejected.** The additive walk (3/7) is the only
candidate that makes abundant cross-offset exact isomorphs, but its σ commute
(1.0) and it is crackable — unlike the real corpus. The single-rotor Enigma
and the 4-rotor progressive make **no** isomorphs at all: two occurrences of a
phrase land in different rotor phases, so their cipher windows are not
isomorphic. This isolates the binding tension for the next iteration: the
mechanism must simultaneously (i) produce many cross-offset EXACT isomorphs
(among tested generators only the additive walk does), (ii) be non-additive
and non-commuting, and (iii) realign at a moderate rate (~0.065). The obvious
rotor/progressive constructions satisfy (ii)–(iii) but fail (i); the additive
walk satisfies (i) but fails (ii). Resolving that tension — a non-additive
mechanism that still yields offset-invariant constant-σ isomorphs — is the
sharply-posed open problem the tournament hands to Step D and beyond.

## Non-abelian-σ test — is it a general rotor? (eighth phase, `tools/rotorfit.py`)

Step B's tension (offset-invariant exact isomorphs, yet non-additive and
non-commuting) has one clean resolution: a progressive key with a GENERAL
(non-cyclic) permutation rotor T, c_i = Tⁱ(B(p_i)), so σ = T^δ between
occurrences at potential offset δ. This is directly testable: read T off a
unit-offset edge in a web, then verify the OTHER edges satisfy σ = T^δ.
Validated both ways (gates pass): a true single rotor fits its independent
edges 42/42 = 1.0; a per-message monoalphabetic cipher (which also makes exact
isomorphs, but with an OFFSET-INDEPENDENT σ) fits only 0.17.

**Result: the general-rotor resolution is not supported.** On the one real web
with a unit-offset edge ({east-5, west-4, east-4}), the independent edges agree
**0/4** with σ = T^δ — nothing like the 1.0 a rotor gives, and consistent with
the offset-INDEPENDENT control. The other web has no unit-offset edge, so it
cannot be read directly. Caveat: this rests on only four discriminating
letter-comparisons (the message pairs share few letters), so it is suggestive,
not iron-clad.

Interpretation: the isomorphs' offset-invariance does NOT come from a
progressive rotor. It is consistent with σ being a fixed PER-OCCURRENCE
substitution difference (σ = Pⱼ∘Pᵢ⁻¹) that is constant within a stretch and
changes per stretch in a POSITION-INDEPENDENT way — the third phase's
"per-stretch fixed substitution", now sharpened: the pieces are not organized
by position/offset, so no rotor or progressive schedule ties them together.
With only 18 isomorph pairs, the per-occurrence substitutions Pᵢ are
under-constrained (each occurrence touches few edges) — a concrete measure of
why the cipher resists: the corpus genuinely under-determines the schedule.
The next leverage must come from constraining the substitution externally
(Step D family-48 focus with a crib, or verifying the game's generation code),
not from more schedule diagnostics on the existing ciphertext.

## Step D — prefix crib route (ninth phase): structurally closed

The shared prefixes are stronger than "isomorphic": positions 1–2 are the *same
ciphertext* `[66,5]` in all nine messages, family A {east-1, east-2, west-1}
shares an *identical* 24-letter prefix (pos 1–24), and family B {east-4,
east-5, west-4} an identical 20-letter prefix (pos 1–20). Two facts kill the
crib plan:

1. **No isomorph edge touches the prefix** — every catalogued occurrence
   starts at position ≥29, so a prefix crib cannot propagate into the σ-web
   (the prefix is an isolated island; the catalog rightly excludes identical
   same-position windows, which carry no σ information).
2. **The prefix has no internal self-consistency gate**: its own repeats sit
   in different keystream phases, so a guessed plaintext there cannot be
   validated against anything — pure speculation, which our discipline forbids.

Two byproducts worth keeping:

- **A pure period-4 keystream is decisively refuted.** Isomorph offsets span
  all residues mod 4 (0:6, 1:2, 2:5, 3:5) with non-identity σ; a period-4
  keystream would force every long isomorph to offset ≡0 mod 4 with
  σ = identity. Period-4 is a keystream *component*, not the whole keystream
  (consistent with the small extra accumulating state, S_eff≈15).
- **Families A and B diverge in content at position 3** (mechanism-agnostic:
  they share plaintext at pos 1–2, so any prior-plaintext-determined keystream
  gives identical c[3] iff p[3] is identical — but c[3] differs). The nine
  messages are genuinely different texts sharing only a short common opening.

## Step E — game-install verification (tenth phase, `tools/wak_unpack.py`): no shipped decoder

Noita *is* installed here; `tools/wak_unpack.py` (pure-Python, format decoded
and validated: 14,745 files, every offset bounds-checked) reads `data.wak`
without launching the game. A thorough bounded search found **no decoder,
plaintext table, or eye-message mural** in the shipped data. The eye-adjacent
assets are unrelated Easter eggs: `caves/eye_0*.png` are grayscale "watching
eye" cave-decoration frames, and the five `eyespot_a..e` + `book_s_a..e` are
the *tripping* books — readable English "Notes on Grand Alchemy" flavor text,
not our cipher. This bounds the project: the shipped data files do not contain
the answer, so the eye-message cipher is authored/externally-keyed and the
ciphertext really is all we have to work from. (Caveat: bounded — the
executable itself and a full text/pixel sweep of all 14,745 files remain
unsearched.) **Superseded by the eleventh phase below, which removes the
bound.**

## Step E, made exhaustive (eleventh phase, `tools/wak_sweep.py` + `tools/eye_mural_scan.py`)

Two instruments that remove the tenth phase's "bounded" caveat, run on a
second machine (the WAK archive lives at a different, machine-specific Steam
path here — `tools/wak_unpack.py`'s `find_wak()` now probes several known
locations instead of one hardcoded path).

1. **Byte-level keyword sweep, all 14,745 files** (`wak_sweep.py`). Gate: a
   synthetic marker string planted at a random byte offset inside a real
   entry's span (SEED=1) must be found and attributed to the right file —
   PASS. Real sweep (English + Finnish cipher/eye terms — Nolla Games is a
   Finnish studio): every hit traces to an unrelated engine or asset string
   — `decode` in a Twitch-integration script, `glyph` in a *trailer* PSD's
   layer name, `rune` in `Runestone` wand-item files, `pupil` in a boss
   enemy's eye sprite, `orientation` as the ubiquitous XML sprite attribute
   (859 hits, 432 files — a generic engine term), `translat` in unrelated
   scripts (essence pickups, perks, the debug menu). No hit for `cipher`,
   `decrypt`, `encrypt`, `alphabet`, `trigram`, `silmä`, `viesti`,
   `salakirjoitus`, or any other cipher-specific term, in English or Finnish.

2. **Pixel-level eye-outline GRID scan, all 9,046 raster files** (9,030 PNG
   + 6 BMP + 6 PSD, `eye_mural_scan.py`). Reuses `transcribe.py`'s own exact
   template (11×7 outline, ≥18/20 outline pixels lit at native resolution)
   and looks for many matches mutually paired at the 12×7 pitch — a GRID,
   not a lone eye sprite (the game has several of those, unrelated to our
   cipher). Three gates: (a) the known source sheets must recover a huge,
   >99%-paired grid (1,347–1,761 matches) — proves the detector works on
   the real thing; (b) seeded random noise at several scales must register
   no grid — calibrates that the raw 18/20 threshold isn't trivially
   satisfied by chance; (c) the catalogued lone-eye Easter eggs
   (`caves/eye_0*.png`, `eyespot.png`, boss eye sprites) must never
   register as a grid. Gate (c) **initially failed**: `eyespot.png` (a
   flat, solid-filled navy circle+triangle map icon) trivially satisfies an
   outline-only template at every interior pixel, since a uniform fill
   lights up any subset of points — and the mechanical de-overlap step
   imposed an artificial 300+-match "grid" on top of that. Fix: the
   detector now also requires the eye shape's always-background corner
   pixels (28 positions per box, computed from `transcribe.py`'s own
   OUTLINE/INSIDE convention) to stay clear, which a solid fill can't
   satisfy but a real almond-shaped eye always does; all three gates pass
   after the fix. Real scan: **zero** of the 9,046 images anywhere in the
   archive register as an eye-outline grid, at native resolution, other
   than the two already-catalogued Easter eggs (which now correctly score
   0 raw matches).

3. **Executable strings** (`strings` on both `noita.exe` and
   `noita_dev.exe`, never searched before): no cipher/eye-message hits
   beyond generic engine text — the `Message_*` entity-event system, a PNG
   decoder's "invalid decoded scanline length", a `base16_decode` debug
   utility, and the `ThreeEyesAreWatchingYou` perk name (an unrelated wand
   perk).

**Net: the exhaustive sweep changes nothing about the tenth phase's
conclusion, but removes its caveat.** The eye-message cipher is not
recoverable from any file in the shipped game or its executables' visible
strings — it is authored/externally-keyed, and the ciphertext (plus, now,
this fully-searched absence) is what the project has to work with.
(Residual caveat: the executables' binary *logic* itself — code implementing
a decoder with no giveaway string — is unsearched; that would need
disassembly, out of scope here.)

## Step B extended — three more candidates (twelfth phase, `tools/tournament.py`)

With the shipped game ruled out as an external source (eleventh phase),
effort returned to the substitution-schedule problem Step B (seventh phase)
left open: no tested small-state generator simultaneously (i) makes
abundant cross-offset EXACT isomorphs, (ii) is non-additive/non-commuting,
and (iii) realigns at the observed rate. Three new `make_*` candidates,
plus a new diagnostic (`delta_mod4_spectrum`, histogramming each raw
isomorph pair's relative offset mod 4 — the same statistic Step D ran on
the real corpus). Gates (existing detector-validation gates) still pass
after adding the candidates and after a fix to the synthetic corpus itself
(below).

**Methodology fix first.** The tournament's shared synthetic-corpus builder
(`corpus()`) inserts each candidate-seeded repeated phrase at hardcoded
positions. Checking those positions: every pairwise gap between insertions
of the *same* phrase (32, 4, 28, 4, 16, 16) was **coincidentally a multiple
of 4** — meaning the harness could never tell a period-4 layer keyed to
absolute position apart from one that isn't, since it never sampled a
cross-message offset outside residue 0. Fixed by shifting one insertion per
phrase group by 1-2 positions so the pairwise gaps span multiple residues
mod 4; re-ran the existing candidates to confirm the detector-validation
gates still pass (they do) before trusting any new candidate's result.

**Candidate 1 — `make_deck_rebuild`: "deck rebuilt, not shifted."** A bank
of 15 (~ Step C's S_eff upper bound) independently-random permutation
tables; the encryptor applies ONE table unchanged across an entire
plaintext-notch-delimited stretch (no progressive exponent, unlike
`make_rotorcycle`'s T^m), and switches to the next table on a notch —
tests the third-phase idea, previously downgraded only for row-locked
timing, with a plaintext-driven trigger instead.

Result: **first tested mechanism to produce abundant, offset-independent
exact isomorphs** (maxL 39; raw-pair delta-mod-4 spectrum {0:34, 1:1, 2:3,
3:4} — spread across residues, like real, because composing two
INDEPENDENT fixed substitution tables is a single well-defined map
regardless of position). The prior naive rotor-cycle produced zero. But it
fails zero-doubles (70 ciphertext doubles — a bare table-swap has no
collision-avoidance mechanism) and period-4 (d4-ratio 1.05 — it has no
period-4 component at all). (Caveat: `walk_dim` and `commute` for this
candidate read high, like the additive walk's — plausibly a small-n
artifact of this harness's few seeded phrases mostly landing in the same
one or two states, since the isomorph-count and delta-spectrum evidence,
built on more direct counts, already shows the tables are NOT globally
translation-equivalent the way the additive walk's deck is.)

**Candidate 2 — `make_period4_position`.** Candidate 1 plus an outer
period-4 layer keyed by ABSOLUTE character position (`Q[i mod 4]` composed
outside the table lookup) — the literal reading of "a period-4 sub-key that
selects among 4 fixed substitutions, composed with a small autokey state."

Result: **isomorph yield collapses from dozens of pairs to 2**, and both of
those 2 sit at delta≡0 mod 4 (spectrum `{0: 2}`) — exactly the predicted
failure mode. Because the period-4 phase depends on absolute position, two
occurrences of the same phrase at a relative offset not divisible by 4 pick
up a DIFFERENT composed map at each position within the window instead of
one constant map, so the window stops being an exact isomorph almost
always. Real isomorph offsets span every residue mod 4 (0:6,1:2,2:5,3:5,
Step D) with non-identity σ — a position-locked period-4 composition
predicts the opposite of what's observed, now demonstrated concretely
rather than just argued.

**Candidate 3 — `make_period4_content`.** Candidate 1 plus a period-4 layer
keyed by a SECOND, faster plaintext notch (mean spacing ~4 letters, tied to
letter identity like the slow notch, not to position) instead of raw
position — operationalizing the open idea from Step A's own writeup: "an
order-4 key STATE feeding the non-commuting substitution, not an order-4
offset on a state walk."

Result: **keeps isomorphs at every residue mod 4** (spectrum `{0:12, 1:1,
2:1, 3:3}` — matches the real pattern qualitatively, unlike candidate 2),
because the fast state advances identically in both occurrences of an
identical phrase (it depends on the phrase's own letters, not on where the
phrase sits), so two occurrences that merely start in the same fast-phase
(~1-in-4 chance) stay synchronized regardless of absolute offset. But its
d4-ratio comes back **1.17 — not elevated** (target >1.5): a
content-triggered state change is an irregular renewal process (mean
interval ~4, not a strict clock), too loose to produce the real corpus's
crisp, single-lag-4 recurrence excess.

**Net — the two-way tension is now three-way, and demonstrated rather than
argued.** Abundant offset-independent isomorphs (i), non-additive structure
(ii), and a crisp period-4 recurrence (iii) have never been produced
together by any construction tried in this project. Every candidate that
gets two of the three loses the third: the rotor-cycle got (ii)+(iii) and
lost (i); the additive walk gets (i) but is additive, failing (ii);
deck-rebuild gets (i)+(ii) but has no (iii); position-locked period-4 adds
(iii) but destroys (i); content-locked period-4 keeps (i) but the (iii) it
adds doesn't reach the observed sharpness. Zero-doubles remains unexplained
by every substitution-swap design tried (none has a collision-avoidance
rule); whether one can be added without re-breaking the isomorph catalog —
the reason collision-skip was refuted for the WALK family — is untested for
the deck-rebuild family specifically and is the most concrete next step.

## σ-web minimum-alphabet factorization (thirteenth phase, `tools/alphacount.py`)

The four sections from here to the community cross-check were done in a
parallel session on a second machine (2026-07-16, branched from the ninth
phase) and merged on 2026-09-28; they post-date nothing above in logic, only
in merge order.

The one question the ciphertext could still answer alone after rotorfit: how
many distinct substitution alphabets do the 18 σ maps force? Each occurrence
carries an unknown bijection P_i; σ = Pⱼ∘Pᵢ⁻¹ on its support; candidate
partitions of occurrences into shared-alphabet classes are checked by closing
the class relations under inverse+composition (functional, injective,
identity-on-diagonal). Exact branch-and-bound gives the minimum class count.
Methodology repair the controls forced: σ entries are restricted to each
pair's LARGEST strictly-overlapping witness-span block — mirrored carrier
repeats in the overrun form disconnected spans that poison full-trim σ.
Gates (all pass): a pool-of-4 quagmire must be recovered exactly (accepts
truth, m_min == true count, unique minimal partition); fresh-alphabet and
deck-translate controls must show no false merges.

Real corpus — the web is FULLY constrained (no FREE pairs at all), and the
per-occurrence-alphabet model holds exactly (zero base violations):

- **Web 0** (family-49: east-5, west-4, east-4): all 3 pairwise DISTINCT →
  exactly **3 alphabets**.
- **Web 1** (family-48, 6 occurrences): 14 pairs DISTINCT, exactly one
  SUPPORTED-EQUAL — **west-1[30:52] ~ east-1[30:49] provably share an
  alphabet** (4 identity letters, derived through the web) → exactly
  **5 alphabets**, unique minimal partition. The one confirmed reuse sits
  precisely in the known east-1~west-1 resync region [29,50).

So 9 occurrences force ≥8 distinct alphabets with exactly one reuse: any
reused-pool schedule needs a pool well above 4 (birthday logic pushes it to
the S_eff≈15–41 scale), and per-stretch alphabets are near-fresh.

## Tournament round 2 — quagmire family and the C83:C41 affine autokey

Six new candidates added to `tools/tournament.py` (gates still pass, plus a
new construction gate: the collision-bump variants must yield zero doubles):

- **quagC (notch-pool quagmire + collision bump)** — advance the pool index
  when the output would repeat — gives structural zero doubles in a
  non-additive mechanism with realistic resync (3/7, tying the additive walk).
- **quagC + period-4 plaintext** (a soft list-delimiter every 4th unit)
  reproduced the d4 excess in the original run (3.0 vs real 2.6), which
  suggested **the period-4 signature could live in the PLAINTEXT layer**.
  *Not supported after the merge re-run* (see the merge note below).
- **gak41 — Group Autokey over C83:C41** (from the community intake: the ONE
  group in their GAK classification nobody ruled out; 41 hidden states, right
  in our small-state window): affine maps x→bx+a mod 83, b in the order-41
  subgroup, state composed per letter, c = D[a]. a≠0 gives structural zero
  doubles; shared running products give offset-invariant exact isomorphs with
  constant AFFINE (non-translation, non-commuting) σ. **Scores 5/7 — the
  first candidate ever to reproduce walk_dim==1** (every prior candidate
  missed it) together with non-commuting σ, exact isomorphs, zero doubles and
  the IoC band. Its two misses are resync rate (parameterization-sensitive —
  the twist elements are a free choice) and period-4 (d4), which is still
  unexplained. *(Update: period-4 is a duplication artifact and is no longer
  scored; gak41 is 5/6, missing only resync.)*

**Merge note (2026-09-28).** This round was run on a second machine
(`enter`, 2026-07-16) on top of the ninth-phase harness, in parallel with the
eleventh/twelfth phases here, and merged afterwards. Its original numbers
were produced BEFORE the twelfth phase's mixed-residue `corpus()` fix, and
the tournament shares one RNG stream across candidates, so every round-2
score was re-run in the merged harness. Results: gak41 is unchanged (5/7,
walk_dim==1; robust over a 20-seed sweep — median 5/7, walk_dim==1 in
14/20). The period-4-in-plaintext claim did NOT survive: quagC + p4 drops
to d4 1.44 (1/7), and over 20 seeds clears d4 > 1.5 in only 5/20 (median
1.14); adding p4 plaintext to gak41 also breaks its walk_dim==1 (2/20).
Period-4 therefore remains the open third leg of the twelfth phase's
three-way tension, now with gak41 holding the other two. *(Superseded:
the fifteenth phase shows period-4 is a duplication artifact, so the
tension is two-way and gak41 misses only resync.)*

**Killer next test:** every solve so far was translation-only (σ = π+δ). The
gak41 model predicts σ_k = π⁻¹(β_k·π + α_k) with β_k in the order-41
subgroup — an affine-conjugate fit no phase has attempted. If the real σ web
is affine-consistent, that is the break; if not, the community's GAK
classification is exhausted for pure group-autokeys.

## Community deep-doc intake (`community-docs/`, adversarial)

Fetched the deep documents linked from the wiki (Lymm's alignments, Toboter's
Progress 2025-12-28, codewarrior0's Analytical Overview + Isomorphism doc,
Dykoine's cipher-model doc, kaliuresis's decompilation guide, the main
reverse-engineering doc). Reconciliation:

- **New to us, adopted after verification:** the eye data is embedded in
  `noita.exe` (kaliuresis documents the 2021-beta encoding and the first
  east-1 u64, `0xacf686745634505c`) — verified and extended by
  `tools/exesearch.py` below. The GAK classification over transitive groups
  {C83, D166, C83:C41, C83:C82, A83, S83} with all but **C83:C41** ruled out
  (their A83/S83 survivor is refuted by our Step C resync bound) — directly
  produced the gak41 tournament candidate above. Toboter's independent state
  lower bound (≥21 states, 10% chance >26) complements our S_eff≈15 scale
  from the other side; 41 sits in the joint window.
- **Superseded by our phases:** their "P^pos[S[char]] → Maybe?" is the
  progressive rotor rotorfit refuted (0/4); their Alberti/Vigenère/affine-
  positional failures are all inside our refuted families; Dykoine's
  brute-force candidates (polynomial, N-time-pad, Alberti) are all marked
  Failure by the community itself.
- **Colour:** Arvi (the developer) confirmed in an interview that "the eye
  decorations do contain a message"; the 0xacf68674 dword happens to be the
  CRC-32b of "lumikki" (Snow White) — a curiosity, unverified as intentional.

## EXE re-derivation — provenance closed (`tools/exesearch.py`)

The eye data was located and fully decoded in the CURRENT build (Jan 2025) of
`noita.exe`: the spawn function materializes u64 constants as mov-immediates
(lo dword at ebp-0x4c, hi at ebp-0x48, hi omitted when zero); each u64 packs
≤22 stream characters as base-7 digits (char+1, digits 1..6, '5' = row/message
separator) with a final ×7 shift. Gates (all pass): the documented 2021 u64
must appear and decode to our east-1 prefix (two fully independent
derivations agreeing); a planted synthetic stream must be recovered exactly;
and the headline —

- **The extracted stream (3,194 chars, one contiguous run in interleaved
  order E1 W1 E2 W2 E3 W3 E4 W4 E5) is byte-for-byte identical to our
  pixel transcription.** The provenance chain pixels → transcription →
  shipped binary is closed at zero mismatches.
- **No other decodable runs exist in the binary** — no hidden tenth message,
  no west-5, in this encoding.
- Completing Step E's bounded caveat: a keyword sweep over all data.wak text
  files finds nothing eye-related, and an exact eye-template pixel sweep over
  all 9,030 wak PNGs finds exactly one match — `data/particles/eye.png`, a
  mundane single-eye particle sprite (which doubles as a detector sanity
  check). The shipped game contains the ciphertext and nothing else.
  (This independently agrees with the eleventh phase's `wak_sweep` /
  `eye_mural_scan` sweep above, and goes past its residual caveat: the
  eye data inside the executable's code is decoded here, not just its
  strings searched.)

## Step F — affine-π solver for gak41 (fourteenth phase, `tools/affinefit.py`)

gak41 predicts every isomorph σ is affine in the hidden order π:
π(b) = β·π(a) + α with β in the order-41 subgroup of Z83* (= the nonzero
quadratic residues). Every earlier solve was translation-only (β = 1). The
solver branches over β per evidence block (41 choices), keeps an incremental
GF(83) echelon system in (π, α), and prunes on inconsistency or any forced
letter merge (identical null-basis rows). Split over 41 processes; run on the
16-core machine.

Two method repairs the controls forced:

- **Witness depth 2.** The isomorph finder extends windows until a repeat
  MISmatches, so windows reach into the region where the plaintexts have not
  yet converged whenever one repeat there matches by chance. On synthetic gak41
  with known π, depth-1 witness spans leave 9/377 blocks non-affine under the
  TRUE π; requiring two distinct witness spans leaves 0/420.
- **Excursions.** Where two messages are identical on both sides of a short
  differing patch (east-1/west-1 at 25–28 and 33–36), the state left and came
  back; no witness can see that, so those positions are cut from every block
  (detector false-positive rate: 1 in 30 synthetic corpora).

Gates (all pass): synthetic gak41 (274 aligned pairs) → exactly one
assignment, the true one, π recovered exactly up to affine gauge; the additive
walk's all-β=1 assignment is accepted; fresh random alphabets (319 pairs) give
no fit within budget (an exhaustive proof takes hours, so a REFUTED real
verdict would additionally require it).

**Real corpus: INCONCLUSIVE.** Uncut, no assignment fits — but leave-one-out
shows the conflict comes only from the two blocks crossing the east-1/west-1
excursions (an affine map fixing 19, 0, 14, 21 cannot send 40→25). With the
excursions cut, the web fits with 2,050+ assignments, null-space dimension 8,
and only the 2 gauge letters pinned (122 aligned pairs over 60 letters vs the
274 that pinned the synthetic truth). gak41 is neither refuted nor recoverable
from the isomorphs alone. Its standing tension is resync: a 3,403-element
group re-collides after differing snippets only on equal group products, yet
the corpus resyncs four times (Step C: S_eff ≲ 15).

A written status report of the whole project to date is in
`reports/noiteyes-status-report-2026-09-28.pdf`.

## Period-4 is a duplication artifact (fifteenth phase, `tools/diagscan.py`)

The d4 "period-4" signature (distance-4 equal-letter recurrence 26 vs ~11.6
expected; d4-ratio 2.6) drove Step A, the period-4 tournament candidates and
one leg of the three-way tension. It does not survive counting each plaintext
event once.

**Exhaustive isomorph scan.** `find_isomorph_pairs` seeds only on 12-letter
windows with ≥3 repeats, so sparse isomorphs were never seeded. `diagscan`
walks every message pair at every relative offset and takes maximal
consistent-bijection windows; evidence = repeats matched on both sides.
Nulls: 20 shuffled corpora give **0** windows even at evidence ≥3; reversing
any one message gives **0**. Real: 28 windows at ≥3, 9 at ≥5, 2 at ≥7 — all
genuine. 16 of the 28 are NEW (not in `data/isomorphs.json`), notably:

| pair | len | evidence |
|---|---|---|
| east-4[68:102] ~ west-4[71:105] | 34 | 7 |
| east-4[68:98] ~ east-5[69:99] | 30 | 5 |
| east-5[66:99] ~ west-4[68:101] | 33 | 5 |
| east-3[63:91] ~ east-4/east-5/west-4 (72/73/75…) | 28 | 3 each |
| east-3[86:109] ~ west-1[55:78] | 23 | 3 |

So family B's tail (after the shared prefix) is one more isomorphic stretch
across all three messages, and **east-3 — previously unlinked — is tied to
that tail and to west-1**. Output: `data/isomorphs_diag.json`
(`data/isomorphs.json` is left untouched for the other tools).

**Deduplicated recurrence.** Union-find over (message, position) merges
same-position identical letters (the family prefixes) and every aligned
position of every evidence ≥3 window; each (class(i−d), class(i)) recurrence
is counted once, against Poisson(slots/83):

| d | 1 | 2 | 3 | **4** | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| raw | 0 | 5 | 9 | **26** | 11 | 12 | 12 | 11 | 18 | 11 | 10 | 7 | 19 | 12 | 4 |
| dedup | 0 | 5 | 3 | **11** | 11 | 10 | 4 | 7 | 12 | 7 | 7 | 7 | 12 | 10 | 4 |
| expected | 7.5 | 7.6 | 7.7 | **7.8** | 7.8 | 7.9 | 7.9 | 8.0 | 8.0 | 8.1 | 8.1 | 8.1 | 8.2 | 8.2 | 8.2 |
| p(≥) | 1 | .87 | .98 | **.16** | .17 | .27 | .96 | .68 | .11 | .69 | .69 | .70 | .12 | .31 | .96 |

d4 is level with d5, d9, d13 — no period. 15 of the 26 raw d4 hits were
copies of the same event across family members and isomorph partners. What
DOES survive is **near-repeat depletion**: d1=0 (p≈5e-4) and d2–d3 low — the
zero-doubles signature, perhaps extending two steps further.

**Transcription error found.** The fourth-phase bullet cited closure-zone
breaks at "mod 4 z≈4.3" as a second, independent period-4 observable. The
saved `data/rowphase_report.json` has `real_mod_zone` mod 4 z=0.97 and mod 5
z=3.64; 4.3–4.4 is instrument 2's d4 count z. There was only ever one
period-4 observable, and it is the duplication artifact above.

**Consequences.** `tools/tournament.py` no longer scores period-4 (moved to
`RETIRED_TESTS`, d4_ratio still reported); **gak41 is now 5/6, missing only
resync**. The three-way tension collapses to two: long cross-offset
isomorphs with non-identity σ (a large or group-structured state) vs fast
resync (S_eff ≲ 15, resting on four structural realignments). Step A's
refutation stands but its motivation is gone.

**Pupil-geometry test (null).** Community avenue: σ maps preserve eye-digit
structure (letter = 25·msd+5·mid+lsd). Over the 18 catalog isomorph pairs with ≥5
non-identity mappings, the fraction of mappings preserving digit k
(k = msd/mid/lsd) is .268/.174/.192 vs a shuffled null of .288/.202/.204
(±.025), and same-digit coherence is .282/.204/.186 vs .282/.200/.199
(±.016). σ carries no pupil geometry.

## Tournament round 3 — small-state machine and hybrid register (sixteenth phase, `tools/tournament3.py`)

Round 3 tests the two mechanisms that could reconcile the fast resync with the
long non-identity isomorphs: a **K-state machine** c_t = T[s_t][p_t]
(K = 8/12/15/20/30, with and without double-avoidance) and a **hybrid
register**: a fast part keyed by the previous N−1 plaintext letters
(N = 2–5), times a slow part with M = 5–30 states that only 6 trigger letters
advance, either additively or by permutation. gak41 is the baseline.
Scoring uses 7 tests over 12 seeds per config. Period-4 is reported but not
scored. It is replaced by **near-repeat depletion**: real raw d2 = 5 and
d3 = 9 against a d5–d10 mean of 12.5 (r2 = 0.40, r3 = 0.72); a candidate
passes when each ratio is ≤ max(0.75, 1.5× real). Deterministic, ~10 s.

**Distinct-σ check on real data (before any simulation).** The 18 σ in
`data/sigmas.json` (2 occurrence webs) fall into **14 mutually incompatible
groups** (compatible = ≥3 agreeing keys, no conflicts, either orientation).
All 5 compatible links share an occurrence, so they are duplication or
composition, not recurrence; one pair (e1[40:49]~e1[68:77]) is listed twice.
**No σ recurs across webs or families** (null expectation 0.000), and **no
real σ has a fixed point**, so no isomorph is a plain repeat. With only 2
webs, the group count alone does not break a K² bound for K ≥ 4.

| candidate | iso L≥25 | walk dim 1 | noncomm | 0 doubles | d2/d3 | resync | IoC | mean /7 |
|---|---|---|---|---|---|---|---|---|
| gak41 | 1.00 | 0.58 | 0.75 | 0.83 | 0.25 | 0.00 | 1.00 | 4.4 |
| state machine K=8–30, no avoid | 0.42–1.00 | 0 | 0 | 0 | ≤0.17 | 0 | 0–0.17 | 0.5–1.1 |
| state machine K=8–30, avoid doubles | 0.50–1.00 | 0 | 0 | 1.00 | ≤0.17 | ≤0.08 | 0–0.50 | 1.7–2.2 |
| hybrid, additive slow part | 0.67–1.00 | 0–0.08 | 0–0.33 | 1.00 | 0–0.42 | 0–0.25 | 0.42–1.00 | 2.5–3.2 |
| hybrid, permutation slow part | 0.67–1.00 | 0–0.08 | 0–0.17 | 1.00 | 0–0.42 | 0–0.08 | 0.17–1.00 | 2.5–3.0 |

Pass rates; no config passes all 7 in any seed. Per-config rows are in
`data/tournament3_report.json`.

**Verdict.**
- **Small-state machine: dead.** Its isomorphs are plain repeats (fixed-point
  fraction ≈1.0, real 0.0). A random transition table moves the state on
  every letter, and a non-identity σ needs the same state *pair* held for 25+
  letters, which never happens. It also fails walk dim, non-commutation and
  resync (0.43–0.57 against the real 0.065) at every K.
- **Hybrid register: fails.** It gets isomorphs, zero doubles and IoC, but
  walk dim is 5–45 instead of 1, resync runs 0.1–0.6, 26–74% of its σ are
  plain repeats, and several configs recur σ across webs (real: never).
  An off-report sweep of the trigger count (2/6/15/30) only traded
  isomorph length against resync.
- **gak41: still best, but fails two tests.** Resync is 0.0 in every seed,
  and it lacks the d2 depletion (r2 ≈ 1.7 against 0.40). Its walk-dim pass is
  7/12, so that test is noisy under this harness. Its round-2 "5/6" did not
  include depletion.
- **Refined constraint set:** rigid walk (dim 1); non-commuting, non-identity σ
  that never recur across webs; a low but nonzero resync rate; d1–d3
  depletion. None of the three families covers it. The core mismatch is
  that σ look like near-fresh alphabets per occurrence (large state), while
  resyncs need a key state that is occasionally held.

## Tournament round 4 — repairing gak41: resets, avoidance, short steps (seventeenth phase, `tools/tournament4.py`)

Round 3 left gak41 in the lead but failing two tests for structural reasons.
**Resync is 0 by construction**: a group autokey's state update is
invertible, so two diverged messages can never re-agree. The real update must
therefore lose information somewhere. **d2 is not depleted** either (r2 ≈ 1.7
against the real 0.40). Round 4 adds three minimal modifications to gak41,
which keeps the same key material per seed, and scores them with the same 7
tests and 12 seeds as round 3 (53 configs):
- **reset**: a trigger plaintext letter returns the state to the identity
  (b, a) = (1, 17). Trigger sets span 0.17% (q z) to 15% (space) of
  plaintext mass.
- **avoid w**: if the output repeats one of the last w outputs, a is bumped
  by a fixed step until it doesn't. In 'state' mode the bump is kept; in
  'out' mode only the emitted letter changes.
- **short steps K**: the additive steps are distinct values in 1..K. On
  untwisted stretches a d-step sum lies in [d, dK], so it cannot be 0 mod 83
  while dK < 83. d-repeats are then structurally suppressed except through
  the 5 twist letters, which gives *soft* depletion rather than zero.

The leading configs are then re-scored over 48 seeds (`confirm_48` in the
report), because 2 all-7 passes in 12 seeds is within noise. The run is
deterministic and takes about 15 s.

| candidate (48 seeds) | iso | walk | noncomm | doubles | d2/d3 | resync | IoC | mean /7 | all 7 | resync mean | r2 / r3 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gak41 | 1.00 | 0.62 | 0.92 | 0.85 | 0.27 | 0.00 | 1.00 | 4.67 | 0 | 0.001 | 1.30 / 1.11 |
| short K=27 | 1.00 | 0.62 | 0.92 | 1.00 | 0.88 | 0.00 | 1.00 | 5.42 | 0 | 0.001 | 0.20 / 0.44 |
| short K=27 + reset q z | 0.98 | 0.60 | 0.83 | 0.96 | 0.88 | 0.02 | 1.00 | 5.27 | 0.02 | 0.033 | 0.21 / 0.45 |
| short K=27 + reset q z j x | 0.98 | 0.50 | 0.71 | 0.98 | 0.88 | 0.15 | 1.00 | 5.19 | 0.08 | 0.061 | 0.20 / 0.46 |
| reset q z j x + avoid w=3 out | 1.00 | 0.38 | 0.73 | 1.00 | 1.00 | 0.15 | 1.00 | 5.25 | 0.08 | 0.062 | 0 / 0 |
| reset m + avoid w=3 state | 0.98 | 0.15 | 0.44 | 1.00 | 1.00 | 0.33 | 0.88 | 4.77 | 0.06 | 0.182 | 0 / 0 |
| *real* | | | | | | | | 7 | | 0.065 | 0.40 / 0.72 |

The 12-seed sweep behind the choice of leaders:
- **Resets need rare triggers.** q z (mass 0.17%) gives resync 0.053 at no
  cost (4.33 vs 4.42). Common triggers destroy the isomorphs: e (10.7%)
  scores 0.83, and space scores 1.17.
- **Avoidance zeroes d2/d3.** avoid w=3 out scores 5.5, but gives r2 = r3 = 0,
  whereas real depletion is partial.
- **Short steps work best at K = 27.** Mean scores are 5.5 at K = 27, 5.0 at
  K = 41 and 5.0 at K = 55. At K = 41 r3 is already 1.11, because 3K > 83.

Per-config rows are in `data/tournament4_report.json`.

**Verdict.**
- **Short additive steps are the best single repair.** They give soft d2/d3
  depletion near the real ratios (0.20 / 0.44 against 0.40 / 0.72), fix
  doubles (1.00), and cost nothing on the other tests. That fits the real
  profile better than avoidance, which removes near-repeats entirely.
  Mechanistic reading: the per-letter key increments are *small* relative
  to the alphabet.
- **Rare resets reach real-like resync on average but not per seed.** With
  4 trigger letters (≈0.5% of plaintext) mean resync is 0.061 against the
  real 0.065, yet only 15% of seeds land in the 0.03–0.12 band. A reset is
  also costly: walk dim falls from 0.62 to 0.50 and non-commutation from
  0.92 to 0.71.
- **Passing all 7 stays rare: at most 8% of seeds for any config.** No
  modification of gak41 reproduces the full fingerprint robustly. Resync is
  the bottleneck, followed by walk-dim noise. The remaining question is what
  lossy state update gives resync ≈0.065 while keeping walk dim 1 and
  non-commuting σ. A plaintext reset is too blunt.

## Tournament round 5 — lossy state updates, and what the resync test measures (eighteenth phase, `tools/tournament5.py`)

Round 4 left one gap: gak41's update is invertible, so it realigns only by
coincidence, and plaintext resets were too blunt. Round 5 tries updates that
lose information in other ways, all on the round-4 short-step gak41 (K = 27),
with 200 seeds per config, run on enter:
- **stw (state-driven twist)**: b is twisted when the current a falls in a set
  of m values (m = 2–16). It multiplies b (`stw-mul`), with or without the
  plaintext twists, or overwrites it (`stw-set`).
- **bset**: plaintext twist letters overwrite b instead of multiplying it.
- **coarse C**: a coarse-step autokey. The state is applied to a fixed letter
  code, c = D[b·e(p) + a], and then updated by the letter's *class* (C
  classes). A same-class substitution changes one ciphertext letter and
  nothing after it.
- **ptwist tw = 1–3**: the lossless baseline with fewer twist letters.

Two side tests come with each config:
- **Unaligned repeats**: exact ciphertext runs ≥4 at different positions. The
  real corpus has none. A reset predicts them, because the text after a
  trigger becomes a fixed function of the plaintext.
- **Resync under a short-edit plaintext model**: 1–4-letter edits instead of
  the harness's 7-letter blocks. The real mismatch islands are 1–4 letters:
  east-1~west-1 at [25,29) and [33,37); east-4~east-5 at [21,25), 28,
  [32,35) and [36,38).

| candidate (200 seeds) | iso | walk | noncomm | doubles | mean /7 | no unaligned repeat | resync band, short edits |
|---|---|---|---|---|---|---|---|
| short K=27 (5 twists) | 1.00 | 0.64 | 0.92 | 1.00 | 5.43 | 0.91 | 0.05 |
| short K=27, 3 twists | 1.00 | 0.69 | 0.84 | 1.00 | 5.53 | 0.85 | 0.12 |
| short K=27, 2 twists | 1.00 | 0.66 | 0.84 | 1.00 | 5.51 | 0.79 | 0.21 |
| short K=27, 1 twist | 1.00 | 0.55 | 0.67 | 1.00 | 5.32 | 0.59 | 0.39 |
| + reset q z j x | 0.99 | 0.52 | 0.76 | 0.96 | 5.25 | 0.55 | 0.05 |
| bset | 1.00 | 0.10 | 0.04 | 1.00 | 4.31 | 0.51 | 0.54 |
| stw-mul m=2 | 0.96 | 0.13 | 0.47 | 1.00 | 4.65 | 0.69 | 0.21 |
| stw (all 12 variants, m=2–16) | 0.02–0.97 | ≤0.16 | ≤0.62 | 1.00 | 2.7–4.6 | 0.69–0.94 | 0.01–0.21 |
| coarse C=3–27 | 1.00 | 0.54–0.58 | 0.81–0.92 | 0.03–0.06 | 3.7 | 0.94–0.99 | 0.08–0.24 |
| *real* | ✓ | ✓ | ✓ | ✓ | 7 | ✓ (0 runs) | 0.0645 |

**Verdict.**
- **State-driven loss is incompatible with exact isomorphs.** Two occurrences
  of a phrase start in different states, so they twist at different letters
  and σ stops being constant. iso L ≥ 25 falls from 0.97 to 0.02 as m grows,
  and walk dim is at most 0.16. Any loss must be triggered by the plaintext
  (or the position).
- **Plaintext-driven loss fails in other ways.** Overwriting b (bset) makes σ
  nearly commute (non-commutation passes in 0.04 of seeds). Resets stay
  bimodal (56% of seeds show zero resync) and predict unaligned repeats that
  the real corpus lacks. Only 9% of seeds are in the resync band with no
  unaligned repeats.
- **The resync test depends on the plaintext model, and round 4's premise was
  too strong.** With short edits like the real islands, the *lossless*
  short-step gak41 realigns by coincidence often enough. With 2–3 twist
  letters it lands in the resync band in 12–21% of seeds (1 twist: 39%, at a
  cost in non-commutation), against 5% for rare resets. So four realignments
  do not show that the state update is lossy. They fit an invertible autokey
  whose diverging plaintext is short and mostly untwisted. Step C's p ≈ 0.03
  against an additive S = 83 walk is the same observation.
- **Coarse-step autokeys keep isomorphs, walk dim and non-commutation but fail
  zero doubles** (≤6%). The output is no longer a pure step of the state, so
  nothing forbids c_t = c_(t+1).
- **Current best: lossless gak41 with short additive steps (1..27) and 2–3
  twist letters**, at 5.5/7 under the harness. Its remaining misses are
  walk-dim noise (≈0.67) and a resync rate that depends on how the real
  plaintexts differ, which is unknown.

## Community cross-check (adversarial, re-derived before use)

Fetched the public write-ups (Noita wiki, the "Unsolved Puzzles" page, the
`ngraham20/NoitaCryptographyResearch` repo) and reconciled every mechanism
claim against our independent results; nothing was adopted un-reproduced.

- **Full agreement, independently reproduced:** 83-letter alphabet with a
  unique valid reading order; zero adjacent doubles; isomorphs ("Lymm's gap
  patterns"); flat frequencies ruling out monoalphabetic substitution; the
  period-4 recurrence (community z=4.04 ≈ our z≈4.4); and — notably — the
  community's ruling-out of cyclic-group / Alberti ciphers via "chaining
  conflicts (commutativity/order)" *converges with* our independent
  refutation of the entire bijective-deck group-walk family (dim-1 core,
  non-commuting σ, no group order 59–127).
- **Their leading model — Group Autokey over the full symmetric group S₈₃
  (a per-letter deck shuffle, 82! hidden states)** — is an unproven
  hypothesis (no decryption exists). It is a *superset* consistent with our
  surviving non-commuting fixed-substitution family, but it leaves the
  shuffle *schedule* unconstrained, which is exactly the open problem; our
  fourth-phase result narrows it (not row-locked; carries a period-4 cycle).
- **Superseded community leads:** the `ngraham20` repo's "incrementing wheel"
  (rotate the ring one step per character) is the additive walk we refuted;
  the "period 4 = structural cipher cycle" framing is too strong given the
  L=33 constant-σ runs — it is a keystream sub-cycle. Toboter's "~86,000
  reading orders" and the "5.8×10⁻¹⁸⁵" improbability are community figures we
  did not reproduce and do not rely on (our exhaustion over the 36
  digit-significance orders already gives the unique 0–82 order).

Net: no community claim contradicts our findings or requires a course change;
our third/fourth-phase work extends the community consensus with validated
refutations it did not have.

**Addendum (2026-09-14):** pulled in the four primary documents this section's
sources pointed to but this project had previously only seen via paraphrase —
Lymm's alignments doc, Toboter's progress log, and codewarrior0's two
cipher-analysis write-ups — now saved verbatim under `community-docs/`. They
add detail, not new mechanism candidates: the specific elimination chain
behind "their leading model is Group Autokey over S₈₃" is a group-theoretic
argument by Lymm/Simplesmiler ruling out the cyclic group C83 (equivalent to
our additive walk) via alphabet-chaining contradictions, then the dihedral
group D166, then (tentatively) the affine group C83:C82 — leaving A83/S83 as
what they had not yet ruled out. Step C above already goes further and rules
that out too, from our own data. Lymm's doc also explicitly retracts its own
"same position implies same plaintext character" assumption as a dead end —
worth remembering given the shared-prefix reasoning in Step D. See
`community-docs/README.md` for the full per-document index.

## Files

- `tools/transcribe.py` — reproducible pipeline + verification report
- `tools/pi_solver.py` — isomorph finder + GF(83) walk-model solver
- `tools/isomorphs.py` — dumps the isomorph catalog to data/
- `tools/sigma_diagnostics.py` — same-day predecessor to `sigma_web.py` below;
  kept for provenance only, see the note under finding 5
- `community-docs/` — the fetched community documents (secondary,
  unverified — see `community-docs/README.md`)
- `tools/sigma_web.py` — σ-diagnostics battery (jump-safe cores, consensus
  prune, modulus sweep) with five synthetic control gates that must pass
  before any real-corpus verdict is read
- `data/messages.json` — full dataset (eye grids + letter sequences + conventions)
- `data/messages.txt` — the 9 letter sequences in ASCII+32
- `data/isomorphs.json` — maximal isomorph pairs with evidence counts
- `data/sigmas.json` — the 18 σ maps with trimmed pairs and witness spans
- `data/sigma_report.json` — full third-phase battery report (controls + real)
- `tools/rowphase.py` — row-phase / keystream-period extractor (row-lock test
  + distance-d autocorrelation) with synthetic schedule gates
- `data/rowphase_report.json` — fourth-phase report (controls + real)
- `tools/phase4.py` — Step A: δ-stratified test of period-4 as an additive
  key index (δ-distribution, repeat alignment, stratified solve) with gates
- `data/phase4_report.json` — Step A report (controls + real)
- `tools/resync.py` — Step C: resync-rate estimator of the effective
  accumulating state, calibrated on ciphers of known state size
- `data/resync_report.json` — Step C report (controls + real)
- `tools/tournament.py` — Step B: small-state generator tournament scoring
  candidate mechanisms against the consolidated real signature fingerprint;
  twelfth phase added `make_deck_rebuild`, `make_period4_position`,
  `make_period4_content`, and the `delta_mod4_spectrum` diagnostic;
  round 2 added the quagmire family (`make_quagA/B/C`), `make_gak41` and
  the `p4ify` plaintext option
- `tools/diagscan.py` — fifteenth phase: exhaustive diagonal isomorph scan
  (shuffled + reversed nulls) and union-find deduplicated recurrence by
  distance; shows period-4 is a duplication artifact
- `data/isomorphs_diag.json` — all evidence ≥3 isomorph windows, NEW-flagged
  against `data/isomorphs.json`
- `data/dedup_recurrence.json` — raw vs deduplicated recurrence, d=1..15
- `tools/tournament3.py` — sixteenth phase: tournament round 3 (K-state
  machine, hybrid register, gak41 baseline) with the distinct-σ check on real
  data and a d2/d3 depletion test replacing period-4
- `tools/tournament3_NOTES.md` — short notes on round 3
- `data/tournament3_report.json` — round-3 report (real σ stats + per-config
  pass rates over 12 seeds)
- `tools/tournament4.py` — seventeenth phase: tournament round 4 (gak41 +
  plaintext resets, output avoidance, short additive steps; 48-seed
  confirmation of the leaders)
- `data/tournament4_report.json` — round-4 report (53 configs × 12 seeds +
  `confirm_48`)
- `tools/tournament5.py` — eighteenth phase: tournament round 5 (state-driven
  twists, b-overwrite, coarse-step autokey, twist count; unaligned-repeat and
  short-edit resync side tests; parallel, 200 seeds, run on enter)
- `data/tournament5_report.json` — round-5 report (23 configs × 200 seeds +
  real realign/unaligned-repeat stats)
- `data/tournament_report.json` — Step B report (candidates + real,
  including the twelfth-phase additions)
- `tools/rotorfit.py` — non-abelian-σ test: fits σ = T^δ per web to test the
  general-rotor hypothesis, validated on rotor vs offset-independent controls
- `data/rotorfit_report.json` — non-abelian-σ report (controls + real)
- `tools/wak_unpack.py` — Step E: pure-Python `data.wak` archive reader
  (list/cat/extract), format self-checked on every read; `find_wak()`
  now probes several known Steam-library paths instead of one hardcoded one
- `tools/wak_sweep.py` — eleventh phase: byte-level keyword sweep of every
  one of `data.wak`'s 14,745 files, gated on a planted-marker control
- `data/wak_sweep_report.json` — eleventh-phase keyword-sweep report
- `tools/eye_mural_scan.py` — eleventh phase: pixel-level eye-outline GRID
  scan of all 9,046 raster files (PNG/BMP/PSD) in `data.wak`, gated on the
  known source sheets (positive), random noise (negative), and the
  catalogued lone-eye Easter eggs (negative)
- `data/eye_mural_scan_report.json` — eleventh-phase pixel-scan report
- `data/sigma_diagnostics.json` — output of the superseded `sigma_diagnostics.py`
- `tools/alphacount.py` — minimum-alphabet factorization of the σ web
  (relation-closure consistency + exact branch-and-bound), gate-validated
- `data/alphacount_report.json` — factorization report (controls + real)
- `tools/exesearch.py` — EXE eye-data extractor + wak keyword/pixel sweeps;
  re-derives all 9 messages from noita.exe byte-for-byte
- `data/exesearch_report.json` — EXE re-derivation report
- `tools/affinefit.py` — Step F: affine-π (β in the order-41 subgroup)
  solver with depth-2 witness blocks, excursion cut, and gak41/translation/
  fresh gates; parallel search
- `data/affinefit_report.json` — Step F report (controls + real, cut/uncut)
- `reports/` — written status reports (PDF)
- `GUIDE.md` — operator's userguide for every tool + the phased attack plan
- `PLAN.md` — session handoff: full findings + prioritized next-phase plan
