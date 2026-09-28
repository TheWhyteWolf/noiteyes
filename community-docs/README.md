# Community documents (secondary — verify before trusting)

Fetched 2026-09-14 (plain-text export of the public Google Docs linked from
the saved wiki page) at the user's request, resolving the "community-doc
intake" decision deferred in PLAN.md. These are **not** primary data —
per the project's working principle, treat every claim here as a lead to
independently re-derive from `data/messages.json`, not as a fact.

- `lymm-alignments.txt` — Lymm, "Eye Message Alignments and Gap Patterns".
  Manual alignment/gap-size chart across all 9 messages (the same chart as
  `2560px-Eye_message_alignments_and_gaps.jpg`, in prose+ASCII form). Notes
  the shared `b%` header, zero adjacent-repeat trigrams, elevated gap-3
  count, and the "aabbb"-style triple-repeat pattern in msgs 1-2-3 and
  7-8-9 that suggests omitted spaces. Explicitly retracts its own
  "same position -> same plaintext letter" assumption ("Note: This
  assumption is incorrect. Disregard the following conclusions.").
- `codewarrior-isomorphism.txt` — codewarrior0, "Isomorphism in Classical
  Ciphers". General reference on which classical-cipher families produce
  causal isomorphs (Vigenere-family / progressive-alphabet / autokey: yes;
  repeating-key Vigenere, monoalphabetic, transposition: no) and the
  alphabet-chaining recovery technique, with worked examples.
- `codewarrior-analytical-overview.txt` — codewarrior0, "Analytical
  Overview". Applies standard cryptanalytic tests (Kasiski, kappa,
  periodic/positional/autocorrelation IoC) to the corpus and rules out
  simple mono/polyalphabetic and periodic models; lands on a
  ciphertext-autokey-family hypothesis (offset-1 CTAK explains the zero
  doubles and elevated distance-4 recurrence) as the leading "hypothetical
  novel cipher" candidate at the time of writing.
- `toboter-progress.txt` — Toboter, "Progress" (Discord community log,
  2025-12-28). Running collection of the Discord server's findings,
  including a **group-theory elimination** (by Lymm/Simplesmiler) that
  supersedes the CTAK hypothesis above: CTAK is equivalent to the cyclic
  group C83 and is ruled out because it can't reproduce the alphabet-
  chaining conflicts actually observed; the dihedral group D166 and
  (tentatively) the affine group C83:C82 are also ruled out; the group
  claimed to survive is A83 or S83 — i.e. a "shuffled deck of 83 cards"
  where each plaintext character selects a (possibly arbitrary,
  non-commuting) permutation, and the ciphertext is a fixed function of
  the resulting state ("Group Autokey", GAK). `toboter-progress.html` is
  the raw export the `.txt` was extracted from.

Fetched 2026-07-16 on a second machine (merged in later; the four documents
above were fetched there too, byte-identical, and deduplicated):

- `kaliuresis-decompilation-guide.txt` — kaliuresis, "Noita Decompilation
  Guide" (Ghidra walkthrough, March 2021 beta). Documents that the eye data
  is embedded in `noita.exe` and gives the first east-1 u64
  (`0xacf686745634505c`) — independently re-derived and extended by
  `tools/exesearch.py`.
- `eye-messages-main.txt` — the main community reverse-engineering document
  (last substantial update 2021-03-25; Ninji's RE work): message placement
  (internally alternating East/West) and the in-game display mechanics.
- `dykoine-cipher-model.txt` — Dykoine, "Cipher model definition and Brute
  force guidelines". Brute-force candidates (polynomial-with-modulo,
  N-time-pad, Alberti), all marked Failure by the author.
- `noita-documents-directory.csv` — the community's index of Noita
  documents (name, description, author, last-updated), used to find the
  others.

## How this reframes our own work

Our own strict additive-walk model in `tools/pi_solver.py` (c_i =
D[x_i], x_i = x_{i-1} + v(p_i) mod 83, i.e. offset-1 CTAK over the
**cyclic** group Z83) is exactly the C83/CTAK model the community's
group-theory argument independently rules out, and for the same
underlying reason: forcing a single global bijection D to explain every
isomorph pair produces contradictions (our ~47 forced letter merges vs.
their "alphabet chaining ... only leads to contradictions"). Two
independent methods — our brute-force GF(83) linear solve on our own
independently-verified transcription, and their manual chaining/group
argument — converge on the same rejection. Neither the group-theory
argument's D166/C83:C82 eliminations nor the A83/S83 candidate itself
have been independently re-derived here yet; the isomorph catalog and
`tools/sigma_diagnostics.py` results are the connecting piece (see
PLAN.md) — composition of the pairwise letter-mappings sigma_AB holds
cleanly wherever tested, which is necessary but not sufficient evidence
for a group-action model of either kind, and does not yet distinguish
abelian from non-abelian.
