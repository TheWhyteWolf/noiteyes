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
- **A period-4 keystream cycle is confirmed, from two independent
  observables.** Distance-4 equal-letter recurrence is elevated (obs 26 vs
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
unsearched.)

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
- `community-docs/` — the four fetched community documents (secondary,
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
  candidate mechanisms against the consolidated real signature fingerprint
- `data/tournament_report.json` — Step B report (candidates + real)
- `tools/rotorfit.py` — non-abelian-σ test: fits σ = T^δ per web to test the
  general-rotor hypothesis, validated on rotor vs offset-independent controls
- `data/rotorfit_report.json` — non-abelian-σ report (controls + real)
- `tools/wak_unpack.py` — Step E: pure-Python `data.wak` archive reader
  (list/cat/extract), format self-checked on every read
- `data/sigma_diagnostics.json` — output of the superseded `sigma_diagnostics.py`
- `GUIDE.md` — operator's userguide for every tool + the phased attack plan
- `PLAN.md` — session handoff: full findings + prioritized next-phase plan
