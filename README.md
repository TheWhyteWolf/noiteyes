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
5. **Sigma-composition diagnostics** (`tools/sigma_diagnostics.py`):
   collapsed the 18 isomorph pairs into 34 directional (msgA,msgB,delta)
   letter-mappings, covering a complete graph on Group A's six occurrences
   and a full triangle on Group B's three. Composing every reachable
   A→B→C triple against the direct A→C mapping agrees exactly in 100+
   cases; the only 6 disagreements trace to one already-flagged
   non-well-defined merge (isomorphs.py's documented maximal-extension-
   overrun artifact), not a new inconsistency. This is necessary but not
   sufficient evidence for a group-action mechanism, and does not
   distinguish an abelian model from a non-abelian one.
6. **Community context** (`community-docs/`, fetched 2026-09-14): the
   Noita Discord community independently reached the same refutation as
   (2), by a different method. Their alphabet-chaining/group-theory
   argument (Lymm, Simplesmiler, in `toboter-progress.txt`) shows the
   additive/cyclic model (C83, equivalent to our offset-1 ciphertext-
   autokey walk) cannot produce the observed chaining contradictions, and
   narrows the surviving candidates to the non-abelian permutation groups
   A83/S83 — a "shuffled 83-card deck" Group-Autokey (GAK) cipher. Not yet
   independently re-derived from our own catalog; see `PLAN.md` for the
   planned commutativity test. `community-docs/README.md` indexes what
   each document claims and what remains unverified.

## Files

- `tools/transcribe.py` — reproducible pipeline + verification report
- `tools/pi_solver.py` — isomorph finder + GF(83) walk-model solver
- `tools/isomorphs.py` — dumps the isomorph catalog to data/
- `tools/sigma_diagnostics.py` — sigma-composition + cycle-structure diagnostics
- `data/messages.json` — full dataset (eye grids + letter sequences + conventions)
- `data/messages.txt` — the 9 letter sequences in ASCII+32
- `data/isomorphs.json` — maximal isomorph pairs with evidence counts
- `data/sigma_diagnostics.json` — pairwise letter-mappings + composition test results
- `community-docs/` — fetched Discord-community documents (secondary, unverified)
- `PLAN.md` — session handoff: full findings + prioritized next-phase plan
