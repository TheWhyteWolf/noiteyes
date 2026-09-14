#!/usr/bin/env python3
"""Step E, non-bounded: full keyword + structural sweep of ALL data.wak entries.

Extends `wak_unpack.py` (index reader) rather than replacing it. The prior
session's Step E was explicitly bounded ("a thorough bounded search" per
GUIDE.md); this script removes the bound by scanning every one of the
14,745 entries' raw bytes for cipher/eye-related keywords (English + Finnish,
since Nolla Games is Finnish). See `tools/eye_mural_scan.py` for the
companion pixel-level sweep (does a real eye-outline GRID -- not just
these keywords -- appear in any image in the archive).

Self-check (the gate, in the same spirit as every other tool here even
though this is a search/audit rather than a statistical test): before
scanning the real archive, we plant a synthetic marker string at a byte
offset chosen uniformly inside a real entry's span (RNG seeded, SEED=1) in
an IN-MEMORY COPY of the blob, run the identical scanner, and assert it is
found and correctly attributed to that entry. If this control fails, the
scan pipeline itself is broken and the real results must not be trusted.

Usage: python3 tools/wak_sweep.py
"""
import bisect
import random
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import wak_unpack as W

SEED = 1

KEYWORDS = [
    # English — cipher/message mechanics
    b'eye message', b'eyemessage', b'eye_message', b'eye messages',
    b'cipher', b'decode', b'decoder', b'decrypt', b'encrypt', b'encode',
    b'glyph', b'rune', b'alphabet', b'pupil', b'orientation', b'trigram',
    b'translat', b'lookup table', b'substitution',
    # English — thematic (puzzle/secret framing)
    b'secret message', b'hidden message', b'ancient text', b'ancient message',
    b'riddle', b'unsolved', b'puzzle text',
    # Finnish equivalents (Nolla Games is a Finnish studio)
    b'silm\xc3\xa4', b'viesti', b'salakirjoitus', b'salainen', b'muinais',
    b'koodi', b'avain', b'tulkkaa', b'ratkaisu',
    # The project's own ASCII rendering of the confirmed shared prefix
    # (chr(32+66), chr(32+5)) = "b%" -- cheap to check, unlikely to hit,
    # but zero-cost and specifically motivated by our own transcription.
    b'b%',
]


def find_all(haystack, needle):
    """All start offsets of `needle` in `haystack` (case-insensitive)."""
    out = []
    start = 0
    hl = haystack.lower()
    nl = needle.lower()
    while True:
        i = hl.find(nl, start)
        if i < 0:
            return out
        out.append(i)
        start = i + 1


def build_offset_index(entries):
    """Sorted list of (start_offset, name, off, length) for bisecting a byte
    position back to its owning file."""
    idx = sorted(entries, key=lambda e: e[1])
    starts = [e[1] for e in idx]
    return starts, idx


def locate(pos, starts, idx):
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0:
        return None
    name, off, length = idx[i][0], idx[i][1], idx[i][2]
    if off <= pos < off + length:
        return name, off, length
    return None


def scan(blob, entries, keywords):
    """Return {keyword: [(name, off, length, local_pos), ...]}"""
    starts, idx = build_offset_index(entries)
    hits = {}
    for kw in keywords:
        positions = find_all(blob, kw)
        located = []
        for p in positions:
            loc = locate(p, starts, idx)
            if loc:
                name, off, length = loc
                located.append((name, off, length, p - off))
        if located:
            hits[kw.decode('utf-8', 'replace')] = located
    return hits


def control_gate(entries, blob):
    """Plant a marker inside a real entry's span; verify the scanner finds it
    and attributes it to the right file. This is the gate: if it fails, the
    byte-offset -> filename mapping (the only non-trivial logic here) is
    broken and the real scan below must not be trusted."""
    rng = random.Random(SEED)
    # pick an entry with enough room for the marker
    candidates = [e for e in entries if e[2] >= 64]
    target = rng.choice(candidates)
    name, off, length = target
    marker = b'PLANTED_NOITEYES_CIPHER_MARKER_XYZ'
    insert_at = off + rng.randrange(0, length - len(marker))
    test_blob = bytearray(blob)
    test_blob[insert_at:insert_at + len(marker)] = marker
    test_blob = bytes(test_blob)

    hits = scan(test_blob, entries, [marker])
    ok = (marker.decode() in hits and
          any(h[0] == name for h in hits[marker.decode()]))
    return ok, name


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else None
    blob, entries = W.read_index(path) if path else W.read_index()
    print(f'index OK: {len(entries)} files, {len(blob)} bytes total')

    # --- control gate ---
    gate_ok, gate_name = control_gate(entries, blob)
    print(f'[gate] planted-marker attribution: {"PASS" if gate_ok else "FAIL"} '
          f'(planted in {gate_name})')
    if not gate_ok:
        print('real verdicts NOT trustworthy')
        sys.exit(1)

    # --- keyword sweep over every byte of every one of the 14,745 files ---
    hits = scan(blob, entries, KEYWORDS)
    print(f'\n[keyword sweep] {len(KEYWORDS)} terms checked against all '
          f'{len(entries)} files (byte-level, case-insensitive)')
    if not hits:
        print('  no keyword matched anywhere in the archive.')
    else:
        for kw, located in hits.items():
            names = sorted(set(n for n, *_ in located))
            print(f'  MATCH {kw!r}: {len(located)} occurrence(s) in {len(names)} file(s)')
            for n in names[:10]:
                print(f'      {n}')

    # PNG pixel content itself is NOT dimension-sniffed here -- an early
    # version of this script flagged PNGs by crude size heuristics (>=100px,
    # a multiple of 20), which produced hundreds of uninformative hits (most
    # game art happens to satisfy that). `tools/eye_mural_scan.py` supersedes
    # this with a gated pixel-level template match (the actual 11x7 eye
    # outline, validated against false positives from solid-fill icons) --
    # run that separately for the image-content half of Step E.
    n_png = sum(1 for e in entries if e[0].lower().endswith('.png'))

    return {
        'gate_ok': gate_ok,
        'keyword_hits': {k: [{'file': n, 'local_offset': p} for n, o, l, p in v]
                          for k, v in hits.items()},
        'n_entries': len(entries),
        'n_png': n_png,
    }


if __name__ == '__main__':
    import json
    report = main()
    with open('data/wak_sweep_report.json', 'w') as f:
        json.dump(report, f, indent=2)
    print('\nwrote data/wak_sweep_report.json')
