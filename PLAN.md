# Noita Eye Messages — session handoff & next-phase plan (2026-07-14)

NOTE: written as a handoff document — the session is being transferred to
another machine. The portable project state lives in /home/voyd/git/noiteyes;
this file should be copied into that directory (e.g. as PLAN.md) so it
travels with the repo. The project README.md covers everything up to and
including the verified transcription; the findings below from the final
analysis phase are recorded ONLY here so far (plan mode prevented updating
the README).

## Context

Goal: crack (or decisively characterize) the unsolved Eye Messages cipher
from Noita. Ground rule from the user: verify everything from primary data;
community material is "secondary truth" to be re-derived before use.

State: we independently transcribed all 9 messages from 1x pixel-art sheets
(3,108 eyes -> 1,036 letters over a contiguous 0-82 alphabet), proved the
unique 1-of-36 reading order, reproduced every community statistical claim,
and verified character-for-character against the community's chart.
Canonical dataset: data/messages.json (+ messages.txt). Tools:
tools/transcribe.py (full pipeline + PASS/FAIL battery), tools/pi_solver.py.

## Findings from the final analysis phase (post-README, this session)

1. Isomorph catalog (real, precisely located):
   - Six-segment mutually-isomorphic web in group A — the same phrase occurs
     TWICE in each of east-1 (~[30:49] & [60:78]), east-2 ([35:57] & [74:98]),
     west-1 ([30:52] & [64:82]); pairwise lengths 18-28.
   - Group B: east-5[35:66] ~ west-4[36:67] (L=31!), east-4[51:66] ~
     east-5[52:67] (L=15), east-4[46:65] ~ west-4[48:67] (L=19).
   - Individual significance 10^-10..10^-15 — these are certain.
2. Strict additive-walk model ("cut deck": c_i = D[x_{i-1} + v(p_i)],
   D bijective) is REFUTED: encoding the isomorph web as walk-translates
   (pi(B_j) = pi(A_j) + delta_p, pi = D^-1) yields a GF(83) linear system
   whose solutions force ~47 letter merges (pi non-injective). Not caused by
   my earlier window-overrun bug (result persists with pattern-verified
   trimmed cores).
3. Resync structure (verified): e1~w1 share [1,25), diverge, re-align
   [29,33) and [37,50); e4~e5 share [1,21) then re-align three times with
   divergence gaps of 4, 1, and 6 letters. Chance for one 3-run ~ (1/83)^3;
   corpus-wide expectation ~0.006 — all real. Any surviving model must
   deliver positional resync after unequal chunks.
4. Position facts: position 0 distinct in all 9; positions 1-2 identical in
   all 9 (values 66,5 = "b%"); family split at position 3 on ADJACENT values
   (48 vs 49). Consistent with numbered-list plaintexts, header outside the
   chain, families = two boilerplate templates.
5. Leading new hypothesis: additive walk + COLLISION-SKIP rule ("if the
   output letter would repeat the previous one, step again"). Structurally
   explains the zero doubles (12.4 expected, 0 observed), predicts the
   community-reported "isomorphs with slightly differing composition"
   (near-isomorphs), and skip-corrupted equations propagating through exact
   Gaussian elimination would produce precisely the mass-merge failure seen
   in (2). Alternatives worth testing: draw-and-relocate deck mechanics
   (structural no-doubles + ciphertext-feedback state); general group walks;
   pure positional keystreams are disfavored (no structural no-doubles, no
   isomorph formation).

## Next-phase plan (for the main PC)

Priority-ordered experiments:

1. Model-free sigma-algebra diagnostics (cheap, decisive):
   extract the pairwise letter-mappings sigma_AB between all isomorph
   occurrences; test composition consistency (sigma_AC = sigma_BC∘sigma_AB);
   characterize the sigmas' cycle structure. Output: what algebraic family
   the state translations belong to (cyclic? conjugate? arbitrary?).
2. Skip-tolerant walk fitting: refit the additive-walk system allowing a
   small number of skip events — RANSAC at the EQUATION level, or DP
   alignment permitting +1 step insertions inside segments. Success metric
   (gauge-invariant): IoC of the recovered difference stream over the whole
   corpus; ~0.012 = fail, >0.05 = layer 1 broken (plaintext-frequency signal).
3. Mechanism-enumeration harness (tools/harness/): plugin interface for
   candidate mechanisms; template-plaintext generator (numbered lists,
   shared boilerplate, per-entry fill-ins incl. transposed variants; EN+FI);
   constraint battery with effect-size tolerances (0 doubles; dist-2 ~0.4x;
   dist-4 ~2.2x; isomorph & near-isomorph rates; resync behavior; flat IoC;
   no off-alignment 3-gram repeats). Keep/kill each mechanism; for survivors,
   attempt key recovery on the real ciphertext (hill-climb with the
   numbered-list cribs). Decisive success = self-confirming readable
   plaintext (EN/FI) consistent across all 9 messages.
4. Community-doc intake (USER DECISION PENDING): the saved wiki HTML contains
   links to Lymm's alignments doc, Toboter's progress doc, and CodeWarrior0's
   "Analytical Overview" / "Isomorphism in Classical Ciphers". Options:
   fetch now (avoid re-treading years of failed attempts; verify-on-adopt) /
   after experiment round 1 / not at all. Also pending: exact statement of
   the community's deck-cipher model (the video transcript is a paraphrase).
5. Game-install verification (Noita was being installed): locate the eye
   renderer/data in the install (transcript claims the engine draws the
   glyphs; images absent from data files); re-derive the 9 sequences from
   the game itself -> hardens provenance chain to primary.
6. Housekeeping: git init the project; copy this plan into the repo; append
   findings (1)-(5) above to README.md; record the isomorph catalog as
   data/isomorphs.json; decide sequential-vs-orchestrated scale for the
   harness (second USER DECISION).

## Verification approach

Every experiment ships as a standalone script in tools/ printing a PASS/FAIL
or quantitative verdict reproducible from data/messages.json alone; claimed
breakthroughs must survive the gauge-invariant IoC test and produce readable,
internally consistent plaintext across all nine messages before being
believed.
