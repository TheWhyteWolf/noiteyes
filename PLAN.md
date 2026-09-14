# Noita Eye Messages — session handoff & next-phase plan (2026-07-14)

> **SUPERSEDED (2026-07-16): see `GUIDE.md`.** Everything below was executed:
> items 1–3 became `tools/sigma_web.py` / `tools/rowphase.py` /
> `tools/phase4.py` / `tools/resync.py` / `tools/tournament.py` /
> `tools/rotorfit.py`, item 5 became `tools/wak_unpack.py` (Step E), and
> item 6 is done. Item 4 (community-doc intake) was executed 2026-07-16.
> Kept for provenance only.

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

## Session note (2026-09-14): reconciled with the 2026-07-16 work

This session started from a stale local clone (last known commit: the
2026-07-14 findings above) and, not yet aware of the 2026-07-16 work now
merged into this file and `README.md`, redundantly re-ran a same-day
predecessor of `tools/sigma_web.py` (kept as `tools/sigma_diagnostics.py`
for provenance) and separately fetched the four community documents named
in item 4 below — Lymm's alignments doc, Toboter's progress log, and
codewarrior0's two write-ups, now saved verbatim under `community-docs/`
(see its README). That fetch is still net-new: this project had previously
only seen the community's group-theory argument via video-transcript
paraphrase (see README's "Community cross-check" addendum). It adds detail
— the specific elimination chain (C83 ruled out by chaining contradictions,
then D166, then tentatively C83:C82, leaving A83/S83) — but no mechanism
this project hadn't already gone on to independently refute in the
2026-07-16 phases (`GUIDE.md` Step C refutes the large-state S83
group-autokey directly from our own data). **`GUIDE.md` is the current
authoritative plan; treat everything below this note as historical.**

## Verification approach

Every experiment ships as a standalone script in tools/ printing a PASS/FAIL
or quantitative verdict reproducible from data/messages.json alone; claimed
breakthroughs must survive the gauge-invariant IoC test and produce readable,
internally consistent plaintext across all nine messages before being
believed.
