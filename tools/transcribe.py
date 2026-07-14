#!/usr/bin/env python3
"""Transcribe the Noita Eye Messages from the 1x pixel-art sheets.

Pipeline (fully self-verifying, no trusted third-party data):
  1. Binarize the 1x sheets (white background / black ink).
  2. Locate every eye by exact match of the 20-pixel almond outline
     (11x7 px glyph, 12 px horizontal pitch, 7 px vertical pitch,
     alternate rows offset by 6 px).
  3. Extract each eye's interior (pupil) pattern; across all 3,108 eyes
     exactly five patterns occur: pupil center/up/right/down/left.
  4. Decode Eyeorientation.png (the community key, 4x nearest-neighbor
     upscale) by majority-pool downscaling; its five glyphs match the
     sheet patterns exactly and give the value mapping 0..4.
  5. Pair rows, group each pair's eyes into alternating down/up triangles
     (down = top,bottom,top in x-order; up = bottom,top,bottom).
  6. Search all 36 digit-significance assignments (6 per triangle type);
     assert exactly one yields trigram values covering 0..82 contiguously:
        down-triangle: msd=top-left,  mid=top-right,   lsd=bottom (lone)
        up-triangle:   msd=bottom-right, mid=bottom-left, lsd=top (lone)
     (both are clockwise walks around the triangle)
  7. Emit per-message letter sequences (0..82), sequenced row-pair by
     row-pair top to bottom, trigrams left to right.
  8. Print a verification report of the community's structural claims.

Outputs: data/messages.json, data/messages.txt
"""
import itertools
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from PIL import Image

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EYE_W, EYE_H, PITCH_X, PITCH_Y = 11, 7, 12, 7

OUTLINE = {0: (4, 5, 6), 1: (2, 3, 7, 8), 2: (1, 9), 3: (0, 10),
           4: (1, 9), 5: (2, 3, 7, 8), 6: (4, 5, 6)}
OUT_PIX = [(r, x) for r, xs in OUTLINE.items() for x in xs]
INSIDE = [(r, x) for r in range(7)
          for x in range(min(OUTLINE[r]) + 1, max(OUTLINE[r]))
          if (r, x) not in OUT_PIX]

NAMES = ['east-1', 'east-2', 'east-3', 'east-4', 'east-5',
         'west-1', 'west-2', 'west-3', 'west-4']


def load(path):
    return np.array(Image.open(path).convert('L')) < 128


def ink_bands(b, axis=1):
    has = b.any(axis=axis)
    out, y = [], 0
    while y < len(has):
        if has[y]:
            y0 = y
            while y < len(has) and has[y]:
                y += 1
            out.append((y0, y))
        else:
            y += 1
    return out


def scan_row(strip):
    """Return x positions of eyes in a 7-row strip (exact outline match)."""
    W = strip.shape[1]
    hits, x = [], 0
    while x <= W - EYE_W:
        if sum(strip[r, x + c] for r, c in OUT_PIX) >= 18:
            hits.append(x)
            x += PITCH_X
        else:
            x += 1
    return hits


def decode_key(path):
    """Map pupil pattern -> value 0..4 from the orientation key image."""
    b = load(path)
    eyes_band = ink_bands(b)[0]
    E = b[eyes_band[0]:eyes_band[1]]
    mapping = {}
    for value, (x0, x1) in enumerate(ink_bands(E.T)):
        eye = E[:, x0:x1]
        s = (x1 - x0) // EYE_W  # integer upscale factor
        assert (x1 - x0) == s * EYE_W, 'key not integer-scaled'
        native = np.zeros((EYE_H, EYE_W), bool)
        for rr in range(EYE_H):
            for cc in range(EYE_W):
                native[rr, cc] = eye[rr * s:(rr + 1) * s, cc * s:(cc + 1) * s].mean() > 0.4
        pat = frozenset((rr, cc) for rr, cc in INSIDE if native[rr, cc])
        outline_ok = all(native[r, c] for r, c in OUT_PIX)
        assert outline_ok, 'key outline mismatch'
        mapping[pat] = value
    assert len(mapping) == 5
    return mapping


def transcribe():
    VAL = decode_key(f'{D}/Eyeorientation.png')
    trigrams = {}   # name -> list of (x-order value triple, 'down'|'up')
    raw = {}
    for sheet, names in [('1x-eyes-east.webp', NAMES[:5]),
                         ('1x-eyes-west.webp', NAMES[5:])]:
        b = load(f'{D}/{sheet}')
        bands = [bd for bd in ink_bands(b) if bd[1] - bd[0] > 20]
        assert len(bands) == len(names)
        for name, (y0, y1) in zip(names, bands):
            assert (y1 - y0) % PITCH_Y == 0
            rows = []
            for r in range((y1 - y0) // PITCH_Y):
                strip = b[y0 + r * PITCH_Y: y0 + (r + 1) * PITCH_Y]
                row = []
                for x in scan_row(strip):
                    pat = frozenset((rr, cc) for rr, cc in INSIDE if strip[rr, x + cc])
                    row.append(((x - 15) // 6, VAL[pat], r % 2))
                rows.append(row)
            raw[name] = [[(c, v) for c, v, _ in row] for row in rows]
            assert len(rows) % 2 == 0
            seq = []
            for p in range(0, len(rows), 2):
                pair = sorted(rows[p] + rows[p + 1])
                cs = [c for c, v, rr in pair]
                assert cs == list(range(len(cs))), f'{name}: pair {p} not contiguous'
                assert len(pair) % 3 == 0
                for i in range(0, len(pair), 3):
                    t = pair[i:i + 3]
                    seq.append((tuple(v for c, v, rr in t),
                                'down' if t[0][2] == 0 else 'up'))
            trigrams[name] = seq
    return raw, trigrams


def find_unique_order(trigrams):
    winners = []
    for pd in itertools.permutations(range(3)):
        for pu in itertools.permutations(range(3)):
            vals = set()
            for n in NAMES:
                for t, typ in trigrams[n]:
                    p = pd if typ == 'down' else pu
                    vals.add(t[p[0]] * 25 + t[p[1]] * 5 + t[p[2]])
            if vals == set(range(83)):
                winners.append((pd, pu))
    assert winners == [((0, 2, 1), (2, 0, 1))], winners
    return winners[0]


def verify(msgs):
    ok = lambda cond, txt: print(('  PASS  ' if cond else '  FAIL  ') + txt)
    all_letters = [v for n in NAMES for v in msgs[n]]
    ok(len(all_letters) == 1036, f'total letters == 1036 (got {len(all_letters)})')
    ok(set(all_letters) == set(range(83)), 'letters cover exactly 0..82')
    firsts = [msgs[n][0] for n in NAMES]
    ok(len(set(firsts)) == 9, f'first letter distinct in all 9: {firsts}')
    ok(len({msgs[n][1] for n in NAMES}) == 1, 'second letter identical in all 9')
    doubles = sum(1 for n in NAMES for a, b in zip(msgs[n], msgs[n][1:]) if a == b)
    ok(doubles == 0, 'no letter ever repeats immediately (expected ~12.4 by chance)')
    d = {k: sum(1 for n in NAMES for i in range(len(msgs[n]) - k) if msgs[n][i] == msgs[n][i + k])
         for k in (2, 4)}
    ok(d[4] >= 20, f'distance-4 repeats elevated: {d[4]} vs ~12 expected')
    ok(d[2] <= 8, f'distance-2 repeats depressed: {d[2]} vs ~12.3 expected')
    occ = defaultdict(list)
    for n in NAMES:
        s = msgs[n]
        for i in range(len(s) - 2):
            occ[tuple(s[i:i + 3])].append((n, i))
    bad = [(g, l) for g, l in occ.items() if len(l) > 1
           and any(a[1] != b[1] or a[0] == b[0]
                   for x, a in enumerate(l) for b in l[x + 1:])]
    ok(not bad, 'no repeated 3-grams outside position-aligned shared sections')


def main():
    raw, trigrams = transcribe()
    pd, pu = find_unique_order(trigrams)
    print(f'unique 0-82 digit order confirmed: down={pd} up={pu}')
    msgs = {}
    for n in NAMES:
        msgs[n] = [t[p[0]] * 25 + t[p[1]] * 5 + t[p[2]]
                   for t, typ in trigrams[n]
                   for p in [pd if typ == 'down' else pu]]
    print('\nverification report:')
    verify(msgs)
    os.makedirs(f'{D}/data', exist_ok=True)
    json.dump({'convention': {
                   'eye_values': 'pupil center=0 up=1 right=2 down=3 left=4 (derived from Eyeorientation.png)',
                   'digit_order': 'down-tri: msd=TL mid=TR lsd=B(lone); up-tri: msd=BR mid=BL lsd=T(lone)',
                   'sequence': 'row pairs top-to-bottom, trigrams left-to-right',
                   'letter': '25*msd + 5*mid + lsd, range 0..82; ASCII rendering = chr(32+letter)'},
               'raw_eye_rows': raw,
               'messages': msgs},
              open(f'{D}/data/messages.json', 'w'), indent=1)
    with open(f'{D}/data/messages.txt', 'w') as f:
        for n in NAMES:
            f.write(f'{n}\t' + ''.join(chr(32 + v) for v in msgs[n]) + '\n')
    print('\nwrote data/messages.json, data/messages.txt')


if __name__ == '__main__':
    sys.exit(main())
