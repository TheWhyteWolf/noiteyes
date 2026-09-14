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

## Session update (2026-09-14): sigma-diagnostics + community docs

Ran priority-1 experiment (`tools/sigma_diagnostics.py`): collapsed the 18
maximal isomorph pairs into 34 directional (msgA,msgB,delta) relations (17
unordered), covering a complete K6 graph on Group A's six occurrences and
a full triangle on Group B's three. Composed every reachable A->B->C triple
against the direct A->C relation: **100+ triples agree exactly**; the only
6 disagreements all trace to a single already-flagged non-well-defined
relation (east-1's own delta=+/-28 self-map disagrees across its two
merged windows by 1 position each direction) — the known
maximal-extension-overrun artifact (isomorphs.py docstring), not a new
inconsistency. Composition of the raw sigma_AB value-maps holds cleanly
wherever a real structural connection exists.

However: composition-holds is a **weak** test — true for any
history-independent state-difference mechanism, abelian or not — so it
does not by itself validate the additive/Z83 model pi_solver.py refutes.

Per user decision, pulled in the four community docs named in the prior
plan (now saved under `community-docs/`, see its README for details and
caveats). Headline finding: **the community independently reached the
same refutation we did, by a different method.** codewarrior0's
Analytical Overview initially favored offset-1 ciphertext-autokey (CTAK)
to explain the zero-doubles / elevated-distance-4 structure; a later
group-theory argument in Toboter's Progress doc (Lymm, Simplesmiler) rules
CTAK out — CTAK is equivalent to the cyclic group C83, and C83 cannot
produce the alphabet-chaining contradictions actually observed when
chaining the isomorphs. They also rule out the dihedral group D166 and
tentatively the affine group C83:C82, leaving **A83 or S83** (arbitrary,
generally non-commuting permutations of an 83-element "deck") as the
surviving candidate family — a Group Autokey (GAK) cipher: each plaintext
character selects a shuffle/permutation, composed onto a running state;
ciphertext is a fixed function of that state (e.g. "top card").

This exactly matches our own result: our additive-walk model IS the C83
special case, and IS what pi_solver.py's GF(83) solve refutes (47 forced
merges). Two independent methods, ours from the raw transcription, theirs
from manual alphabet-chaining, land on the same rejection. Neither their
D166/C83:C82 eliminations nor the A83/S83 candidate itself have been
independently re-derived here — that re-derivation, not the skip-corrupted
walk idea from the last session (superseded: skip-corruption is one
special case of "non-additive state update," and the community's argument
is more general and already reached by clean elimination), is the load-
bearing next step.

## Next-phase plan (updated 2026-09-14)

Priority-ordered experiments:

1. DONE: model-free sigma-algebra diagnostics (`tools/sigma_diagnostics.py`)
   — see session update above. Composition holds cleanly but doesn't
   distinguish abelian from non-abelian; superseded as the decisive test
   by (2).
2. DONE: community-doc intake — see session update above and
   `community-docs/`.
3. GAK / non-abelian group-action test (new top priority, motivated by (2)):
   independently re-derive the community's alphabet-chaining argument from
   our own catalog rather than trusting it secondhand. Concretely: pick two
   isomorph relations sharing a common message and domain overlap but
   arising from *different* underlying plaintext segments (not just
   different windows of the same delta), and test whether the induced
   value-permutations commute — non-commutation is direct, independent
   confirmation of a non-abelian (A83/S83-style) mechanism and would settle
   what the community reached only by manual chaining. If it holds, the
   real work starts: model state as a permutation of Z83, plaintext
   characters as (possibly a small alphabet of) permutations composed onto
   it, ciphertext as a fixed readout function; attempt alphabet-chaining
   ourselves on the full corpus (not just messages 1-3 as the community
   did) to see if it resolves further under the larger K6/triangle
   structure our sigma_diagnostics.py already extracted.
4. Skip-tolerant walk fitting: DEPROIRITIZED — skip-corruption is a special
   case of "non-additive state update" and the community's group-theory
   argument already supersedes it with a cleaner elimination. Revisit only
   if (3) fails to find a working non-abelian model.
5. Mechanism-enumeration harness (tools/harness/): plugin interface for
   candidate mechanisms; template-plaintext generator (numbered lists,
   shared boilerplate, per-entry fill-ins incl. transposed variants; EN+FI);
   constraint battery with effect-size tolerances (0 doubles; dist-2 ~0.4x;
   dist-4 ~2.2x; isomorph & near-isomorph rates; resync behavior; flat IoC;
   no off-alignment 3-gram repeats). Keep/kill each mechanism; for survivors,
   attempt key recovery on the real ciphertext (hill-climb with the
   numbered-list cribs). Decisive success = self-confirming readable
   plaintext (EN/FI) consistent across all 9 messages.
6. Game-install verification (Noita was being installed): locate the eye
   renderer/data in the install (transcript claims the engine draws the
   glyphs; images absent from data files); re-derive the 9 sequences from
   the game itself -> hardens provenance chain to primary.
7. Housekeeping: decide sequential-vs-orchestrated scale for the harness
   (USER DECISION, still pending).

## Verification approach

Every experiment ships as a standalone script in tools/ printing a PASS/FAIL
or quantitative verdict reproducible from data/messages.json alone; claimed
breakthroughs must survive the gauge-invariant IoC test and produce readable,
internally consistent plaintext across all nine messages before being
believed.
