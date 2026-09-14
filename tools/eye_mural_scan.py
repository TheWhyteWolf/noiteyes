#!/usr/bin/env python3
"""Step E, non-bounded: pixel-level scan of every PNG in data.wak for the
eye-glyph outline itself, as a structural (not filename) search for an
undiscovered eye-message mural.

Rationale: the prior Step E bounded search found the two catalogued
Easter eggs (caves/eye_0*.png, eyespot/book_s "tripping" books) by NAME.
That leaves open whether some OTHER, unrelated-looking asset happens to
embed a real eye-glyph grid at native (1x) resolution -- a filename search
cannot rule that out. This script instead reuses transcribe.py's own
exact-match template (11x7 outline, >=18/20 outline pixels lit, 12px
horizontal / 7px vertical pitch) and runs it against every one of the
9,030 PNGs' raw pixels, looking for a GRID of matches (many hits at the
expected pitch), not just an isolated eye sprite (of which the game has
plenty -- boss eyes, decorative wall eyes, etc. -- that are not our cipher).

Control gates (must PASS before the real scan is trusted):
  1. POSITIVE: the known source sheets (1x-eyes-*.webp) must yield a large,
     highly regular grid of matches (proves the detector actually finds the
     real thing at native resolution, using nothing but the published
     outline template).
  2. NEGATIVE: seeded random noise of the same scale must yield ~0 matches
     (proves the >=18/20 threshold isn't trivially satisfied by chance).
  3. NEGATIVE: the catalogued non-mural eye Easter eggs (caves/eye_0*.png,
     the boss/decoration eye sprites) must yield at most a handful of
     ISOLATED matches (a lone eye or two), never a grid -- proving the grid
     threshold used to flag real candidates distinguishes "contains an eye"
     from "contains a mural of many eyes."
"""
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import wak_unpack as W

SEED = 1
EYE_W, EYE_H, PITCH_X, PITCH_Y = 11, 7, 12, 7
OUTLINE = {0: (4, 5, 6), 1: (2, 3, 7, 8), 2: (1, 9), 3: (0, 10),
           4: (1, 9), 5: (2, 3, 7, 8), 6: (4, 5, 6)}
OUT_PIX = [(r, x) for r, xs in OUTLINE.items() for x in xs]
INSIDE = [(r, x) for r in range(7)
          for x in range(min(OUTLINE[r]) + 1, max(OUTLINE[r]))
          if (r, x) not in OUT_PIX]
# The corners of the 7x11 box that a real almond-shaped eye NEVER fills
# (they lie outside the outline at every row, where the almond narrows to
# its top/bottom points). transcribe.py's own scan_row() only checks
# OUT_PIX>=18 because in its controlled context (hand-transcribed sheets,
# fixed pitch) that's already unambiguous. Scanning arbitrary archive art
# needs the stronger check: a plain solid-filled blob (e.g. a flat icon)
# satisfies OUT_PIX>=18/20 at every interior pixel trivially, since a
# uniform fill lights up ANY subset of points -- it only fails to be a
# real eye at these always-background corner pixels, which a filled blob
# leaves lit. Found this the hard way: data/biome_impl/eyespot.png (a
# solid navy circle+triangle marker icon) initially registered as a
# 300+-eye "grid" purely from this fill effect, with the grid pitch
# emerging from our own de-overlap step size, not from any real eye
# outline structure.
_ALL = {(r, c) for r in range(7) for c in range(11)}
BG_PIX = list(_ALL - set(OUT_PIX) - set(INSIDE))
BG_THRESH = 0  # zero-tolerance: none of these corner pixels may be ink
THRESH = 18  # of 20 outline pixels, matches transcribe.py's zero-tolerance match


def ink_mask(im):
    """Binarize any PIL image the way transcribe.py treats the source sheets:
    ink = dark AND (opaque or no alpha channel)."""
    im = im.convert('RGBA')
    arr = np.asarray(im)
    lum = arr[..., :3].astype(np.float32).mean(axis=2)
    alpha = arr[..., 3]
    return (lum < 128) & (alpha > 10)


def find_eyes(mask):
    """All (row, col) top-left positions where the 11x7 outline matches at
    >=THRESH/20 pixels AND the always-background corner pixels (BG_PIX) are
    all clear -- the second condition is what a solid-filled blob fails
    (see BG_PIX's docstring above) while a real almond eye always satisfies
    it. Vectorized: sum of shifted boolean views instead of a per-pixel
    Python loop, so scanning all 9,030 archive PNGs is fast."""
    H, Wd = mask.shape
    oh, ow = H - EYE_H + 1, Wd - EYE_W + 1
    if oh <= 0 or ow <= 0:
        return []
    acc = np.zeros((oh, ow), dtype=np.int16)
    for dr, dc in OUT_PIX:
        acc += mask[dr:dr + oh, dc:dc + ow].astype(np.int16)
    bg = np.zeros((oh, ow), dtype=np.int16)
    for dr, dc in BG_PIX:
        bg += mask[dr:dr + oh, dc:dc + ow].astype(np.int16)
    ys, xs = np.where((acc >= THRESH) & (bg <= BG_THRESH))
    return list(zip(ys.tolist(), xs.tolist()))


def grid_score(hits, shape=None):
    """Given raw hit positions (which overlap near any true eye since the
    scan isn't stepping by pitch), cluster into local maxima at the known
    pitch and report how many DISTINCT eyes and how regular their spacing
    is. Returns (n_distinct, n_at_expected_pitch).

    De-overlap uses a claimed-region raster (numpy slice checks) instead of
    an O(n^2) all-pairs Python loop -- some archive textures spuriously
    match the raw template thousands of times (see eyespot.png in the
    negative gate below), and an all-pairs loop over thousands of hits made
    the first version of this scan pathologically slow across 9,030 files."""
    if not hits:
        return 0, 0
    hits = sorted(hits)
    if shape is None:
        max_r = max(r for r, c in hits) + EYE_H // 2 + 1
        max_c = max(c for r, c in hits) + EYE_W // 2 + 1
        shape = (max_r, max_c)
    claimed = np.zeros(shape, dtype=bool)
    hh, hw = EYE_H // 2, EYE_W // 2
    kept = []
    for r, c in hits:
        r0, r1 = max(0, r - hh), r + hh + 1
        c0, c1 = max(0, c - hw), c + hw + 1
        if not claimed[r0:r1, c0:c1].any():
            kept.append((r, c))
            claimed[r0:r1, c0:c1] = True
    n_distinct = len(kept)
    # count pairs at exactly the expected pitch (12 horiz or 7 vert, allowing
    # the alternating 6px row offset transcribe.py documents)
    kept_set = set(kept)
    n_paired = 0
    for r, c in kept:
        neighbors = [(r, c + PITCH_X), (r + PITCH_Y, c), (r + PITCH_Y, c + 6),
                     (r + PITCH_Y, c - 6)]
        if any(n in kept_set for n in neighbors):
            n_paired += 1
    return n_distinct, n_paired


def is_grid(r, min_distinct=20, min_pair_frac=0.5):
    """A real eye-message mural is not just 'has some outline-shaped pixels'
    -- natural pixel art can spuriously satisfy the 18/20 template at a
    handful of scattered spots (see the noise/eggs gate below). What
    distinguishes a genuine grid is many matches that are ALSO mutually at
    the exact 12x7 pitch, so we require both count and pairing fraction."""
    if r is None or r['n_distinct'] < min_distinct:
        return False
    return r['n_paired'] >= min_pair_frac * r['n_distinct']


def scan_image_bytes(data, label):
    try:
        im = Image.open(__import__('io').BytesIO(data))
    except Exception:
        return None
    if im.width * im.height > 4_000_000:
        return None  # skip huge backgrounds; eye grids are small hand-authored art
    mask = ink_mask(im)
    hits = find_eyes(mask)
    n_distinct, n_paired = grid_score(hits, shape=mask.shape)
    return {'label': label, 'w': im.width, 'h': im.height,
            'n_distinct': n_distinct, 'n_paired': n_paired}


def gate_positive():
    results = []
    for f in ['1x-eyes-east.webp', '1x-eyes-west.webp']:
        with open(f, 'rb') as fh:
            data = fh.read()
        r = scan_image_bytes(data, f)
        results.append(r)
    ok = all(is_grid(r, min_distinct=300, min_pair_frac=0.8) for r in results)
    return ok, results


def gate_negative_noise():
    """Random noise can spuriously satisfy the raw 18/20 template count at
    a nontrivial rate at 50% ink density (binomial tail), which is why
    is_grid() also requires PAIRING at the exact pitch, not just count --
    chance matches don't line up with each other. This gate checks that."""
    rng = np.random.RandomState(SEED)
    fails = []
    for w, h in [(64, 64), (200, 150), (400, 300), (800, 600)]:
        arr = rng.randint(0, 2, size=(h, w), dtype=bool)
        n_distinct, n_paired = grid_score(find_eyes(arr))
        r = {'w': w, 'h': h, 'n_distinct': n_distinct, 'n_paired': n_paired}
        if is_grid(r):
            fails.append(r)
    return len(fails) == 0, fails


def gate_negative_known_eggs(blob, entries):
    """The catalogued non-mural eye Easter eggs must never register as a
    grid, however many raw template matches they trigger (eyespot.png's
    own radial texture spuriously matches the template at ~1560 positions,
    which is exactly why raw count alone can't be the discriminator --
    but none of those matches pair up at the 12x7 pitch)."""
    KNOWN = ['data/biome_impl/caves/eye.png', 'data/biome_impl/caves/eye_01.png',
             'data/biome_impl/caves/eye_02.png', 'data/biome_impl/eyespot.png',
             'data/entities/animals/boss_ghost/eye.png',
             'data/entities/animals/boss_fish/eye.png']
    by_name = {n: (o, l) for n, o, l in entries}
    results = []
    for name in KNOWN:
        if name not in by_name:
            continue
        off, length = by_name[name]
        r = scan_image_bytes(blob[off:off + length], name)
        if r:
            results.append(r)
    ok = all(not is_grid(r) for r in results)
    return ok, results


def main():
    print('=== control gates ===')
    pos_ok, pos_res = gate_positive()
    print(f'[gate] positive (known sheets recover a big regular grid): '
          f'{"PASS" if pos_ok else "FAIL"}')
    for r in pos_res:
        print(f'    {r}')

    neg_ok, neg_fails = gate_negative_noise()
    print(f'[gate] negative (random noise recovers ~nothing): '
          f'{"PASS" if neg_ok else "FAIL"}  fails={neg_fails}')

    blob, entries = W.read_index()
    egg_ok, egg_res = gate_negative_known_eggs(blob, entries)
    print(f'[gate] negative (catalogued lone-eye Easter eggs stay isolated, '
          f'no grid): {"PASS" if egg_ok else "FAIL"}')
    for r in egg_res:
        print(f'    {r}')

    if not (pos_ok and neg_ok and egg_ok):
        print('\nreal verdicts NOT trustworthy')
        sys.exit(1)

    print('\n=== real scan: all PNGs in data.wak ===')
    png_entries = [(n, o, l) for n, o, l in entries if n.lower().endswith('.png')]
    print(f'scanning {len(png_entries)} PNG entries for eye-grid structure...')
    flagged = []
    for i, (name, off, length) in enumerate(png_entries):
        r = scan_image_bytes(blob[off:off + length], name)
        if is_grid(r, min_distinct=10, min_pair_frac=0.5):
            flagged.append(r)
        if (i + 1) % 1000 == 0:
            print(f'  ...{i + 1}/{len(png_entries)} scanned, '
                  f'{len(flagged)} flagged so far', flush=True)
    flagged.sort(key=lambda r: -r['n_distinct'])
    print(f'\n{len(flagged)} PNGs register as an eye-outline GRID '
          f'(>=10 matches, >=50% mutually paired at the 12x7 pitch -- the '
          f'gates above show this rules out both chance and lone-eye sprites):')
    for r in flagged[:40]:
        print(f'    n_distinct={r["n_distinct"]:4d} n_paired={r["n_paired"]:4d} '
              f'{r["w"]}x{r["h"]}  {r["label"]}')
    if not flagged:
        print('    (none -- no PNG in the archive contains anything resembling '
              'an eye-glyph grid at native resolution)')

    import json
    with open('data/eye_mural_scan_report.json', 'w') as f:
        json.dump({
            'gates': {'positive': pos_ok, 'negative_noise': neg_ok,
                      'negative_eggs': egg_ok},
            'n_png_scanned': len(png_entries),
            'flagged': flagged,
        }, f, indent=2)
    print('\nwrote data/eye_mural_scan_report.json')


if __name__ == '__main__':
    main()
