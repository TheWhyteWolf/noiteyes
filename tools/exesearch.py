#!/usr/bin/env python3
"""Step E completion: re-derive the eye messages from the shipped noita.exe.

Community lead (Toboter's progress doc + kaliuresis's decompilation guide,
fetched to community-docs/): the eye data is embedded in the executable as 64-bit
integers whose base-7 digits, minus 1, give the eye-value stream (0-4 = pupil
value, 5 = row separator); the 2021-beta first chunk of east-1 was documented
as 0xacf686745634505c. In the current build (Jan 2025) those u64s are not a
flat array but IMMEDIATE OPERANDS in the spawn function's code:

    c7 45 b4 <lo32>          mov [ebp-0x4c], lo
    ...                      (glue varies: push/lea, sometimes a jmp)
    c7 45 b8 <hi32>          mov [ebp-0x48], hi   -- OMITTED when hi == 0
    e8 <rel32>               call append_base7_digits

Each u64 packs the next <=22 stream characters greedily as base-7 digits
(char+1, so digits are 1..6, 6 = row break) with a final *7 shift (units
digit 0); values that fit in 32 bits (each message's final partial chunk)
get only the lo mov. The messages sit in ONE contiguous stream in the
interleaved order E1 W1 E2 W2 E3 W3 E4 W4 E5, separated by '5' characters
(the same character that separates rows).

The tool extracts every such pair from the binary, decodes runs of
consecutive chunks, and compares them against data/messages.json (our
pixel-transcription). Gates:
  G1 cross-source anchor -- the documented 2021 u64 appears as some run's
     first chunk and decodes to a prefix of OUR east-1 transcription
     (two fully independent derivations of the same data).
  G2 scanner completeness -- a synthetic stream encoded into the same
     instruction pattern and planted in random bytes is recovered exactly.
  G3 re-derivation -- all nine expected streams are found in noita.exe,
     byte-for-byte equal to the transcription.
Also reported (not gated): decodable runs that do NOT match the nine known
messages (candidate extra data), a keyword sweep plus an exact eye-outline
pixel sweep over every PNG inside data.wak (completing the bounded Step E
search), and the same extraction attempted on noita_dev.exe.
"""
import io
import json
import os
import random
import re
import struct
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_solver
import wak_unpack
from transcribe import EYE_H, EYE_W, INSIDE, OUT_PIX

ROOT = pi_solver.D
SEED = 1
ANCHOR_2021 = 0xacf686745634505c

MOV_RE = re.compile(rb'\xc7\x45([\xb4\xb8])(....)', re.DOTALL)


def load_streams():
    d = json.load(open(f'{ROOT}/data/messages.json'))
    out = {}
    for n, rows in d['raw_eye_rows'].items():
        vals = []
        for i, row in enumerate(rows):
            if i:
                vals.append(5)
            vals.extend(v for _x, v in row)
        out[n] = vals
    return out


def decode_u64(n):
    """Base-7 digits (MSD first) -> stream chars, or None if not eye data
    (must end in the *7 shift digit 0, all payload digits 1..6)."""
    ds = []
    while n:
        n, r = divmod(n, 7)
        ds.append(r)
    ds.reverse()
    if len(ds) < 3 or ds[-1] != 0:
        return None
    payload = ds[:-1]
    if any(not 1 <= d <= 6 for d in payload):
        return None
    return [d - 1 for d in payload]


def encode_stream(chars):
    """Greedy inverse of decode_u64: pack chars into u64 chunks."""
    chunks, cur = [], 0
    ncur = 0
    for c in chars:
        nxt = cur * 7 + (c + 1)
        if nxt * 7 >= 1 << 64:
            chunks.append(cur * 7)
            cur, ncur = c + 1, 1
        else:
            cur, ncur = nxt, ncur + 1
    if ncur:
        chunks.append(cur * 7)
    return chunks


def extract_runs(blob, max_gap=64):
    """Chunk stream from the mov-immediate pattern: every `mov [ebp-0x4c],
    imm32` is a chunk's low dword; a following `mov [ebp-0x48], imm32` within
    40 bytes supplies the high dword (omitted by the compiler when zero).
    Decodable chunks are grouped into runs of nearby code."""
    movs = [(m.start(), m.group(1)[0], struct.unpack('<I', m.group(2))[0])
            for m in MOV_RE.finditer(blob)]
    hits = []
    for i, (off, disp, lo) in enumerate(movs):
        if disp != 0xb4:
            continue
        hi = 0
        if (i + 1 < len(movs) and movs[i + 1][1] == 0xb8
                and movs[i + 1][0] - off <= 40):
            hi = movs[i + 1][2]
        n = (hi << 32) | lo
        dec = decode_u64(n)
        if dec is not None:
            hits.append((off, n, dec))
    runs = []
    for off, n, dec in hits:
        if runs and off - runs[-1]['end'] <= max_gap:
            runs[-1]['chunks'].append(n)
            runs[-1]['stream'].extend(dec)
            runs[-1]['end'] = off + 26
        else:
            runs.append({'start': off, 'end': off + 26,
                         'chunks': [n], 'stream': list(dec)})
    return runs


def gate(name, ok):
    print(f'  {"PASS" if ok else "FAIL":4s}  {name}')
    return ok


# ------------------------------------------------------------- wak sweeps

KEYWORDS = [b'lumikki', b'silm', b'eye_message', b'eyemessage', b'trigram',
            b'cipher', b'glyph', b'salaus', b'viesti']

TEXT_EXT = ('.xml', '.lua', '.txt', '.csv', '.frag', '.vert', '.inc')


def wak_keyword_sweep(blob, entries):
    hits = []
    for name, off, length in entries:
        if not name.lower().endswith(TEXT_EXT):
            continue
        data = blob[off:off + length].lower()
        for kw in KEYWORDS:
            if kw in data:
                hits.append((name, kw.decode()))
    return hits


def eye_template_count(img_bool):
    """Positions matching a full eye: all 20 outline pixels ink, the four
    corners of the 7x11 box clear, and a pupil-sized interior (1-5 ink pixels
    inside the outline). The corner/interior constraints reject solid-ink
    regions, which trivially satisfy the outline alone."""
    H, W = img_bool.shape
    if H < EYE_H or W < EYE_W:
        return 0
    h, w = H - EYE_H + 1, W - EYE_W + 1
    out = np.zeros((h, w), np.uint8)
    for r, c in OUT_PIX:
        out += img_bool[r:r + h, c:c + w]
    corners = np.zeros((h, w), np.uint8)
    for r, c in ((0, 0), (0, EYE_W - 1), (EYE_H - 1, 0), (EYE_H - 1, EYE_W - 1)):
        corners += img_bool[r:r + h, c:c + w]
    inner = np.zeros((h, w), np.uint8)
    for r, c in INSIDE:
        inner += img_bool[r:r + h, c:c + w]
    return int(((out == 20) & (corners == 0) & (inner >= 1) & (inner <= 5)).sum())


def wak_png_sweep(blob, entries):
    """Exact eye-outline search over every PNG in the archive. The outline is
    matched in both polarities (ink dark or ink light)."""
    found, n_png, failed = [], 0, 0
    for name, off, length in entries:
        if not name.lower().endswith('.png'):
            continue
        n_png += 1
        try:
            img = Image.open(io.BytesIO(blob[off:off + length])).convert('L')
        except Exception:
            failed += 1
            continue
        a = np.array(img)
        for pol in (a < 128, a >= 128):
            k = eye_template_count(pol)
            if k:
                found.append((name, k))
                break
    return found, n_png, failed


# ------------------------------------------------------------------ main

def main():
    rng = random.Random(SEED)
    streams = load_streams()
    report = {}
    gates = []

    # the install dir is machine-specific; derive it from wak_unpack's probe
    gamedir = os.path.dirname(os.path.dirname(wak_unpack.find_wak()))
    exe = open(f'{gamedir}/noita.exe', 'rb').read()
    runs = extract_runs(exe)
    print(f'noita.exe: {len(exe)} bytes, {len(runs)} decodable chunk-run(s), '
          f'sizes {[len(r["stream"]) for r in runs]}')

    # G1 -- cross-source anchor
    anchor_run = next((r for r in runs if r['chunks'][0] == ANCHOR_2021), None)
    ok1 = (anchor_run is not None
           and anchor_run['stream'][:22] == streams['east-1'][:22])
    gates.append(gate('G1: documented 2021 anchor u64 present and decodes to '
                      'our east-1 prefix', ok1))

    # G2 -- scanner completeness on planted synthetic data (both mov forms)
    synth = [rng.randrange(6) for _ in range(150)]
    code = b''
    for n in encode_stream(synth):
        code += b'\xc7\x45\xb4' + struct.pack('<I', n & 0xffffffff)
        code += b'\x50\x8d\x4d\xa8'
        if n >> 32:
            code += b'\xc7\x45\xb8' + struct.pack('<I', n >> 32)
        code += (b'\xe8' + bytes(rng.randrange(256) for _ in range(4))
                 + b'\x8d\x45\xb4')
    noise = bytes(rng.randrange(256) for _ in range(200000))
    planted = noise[:100000] + code + noise[100000:]
    prec = extract_runs(planted)
    ok2 = any(r['stream'] == synth for r in prec)
    gates.append(gate('G2: planted synthetic stream recovered exactly', ok2))

    # G3 -- all nine messages re-derived, in the interleaved order, from one
    # contiguous '5'-separated stream
    ORDER = ['east-1', 'west-1', 'east-2', 'west-2', 'east-3',
             'west-3', 'east-4', 'west-4', 'east-5']
    expected = []
    for i, n in enumerate(ORDER):
        if i:
            expected.append(5)
        expected.extend(streams[n])
    big = max(runs, key=lambda r: len(r['stream'])) if runs else None
    ok3 = big is not None and (big['stream'] == expected
                               or big['stream'] == expected + [5])
    gates.append(gate('G3: noita.exe eye stream == transcription '
                      '(all 9 messages, interleaved order)', ok3))
    if big:
        print(f'    stream at {big["start"]:#x}: {len(big["stream"])} chars '
              f'in {len(big["chunks"])} chunks (expected {len(expected)})')
    report['stream'] = {'offset': hex(big['start']) if big else None,
                        'chars': len(big['stream']) if big else 0,
                        'order': ORDER, 'matches_transcription': ok3}

    # extras: decodable runs elsewhere in the binary (candidate hidden data)
    extras = [r for r in runs if r is not big and len(r['stream']) >= 40]
    print(f'  other decodable runs with >=40 chars: {len(extras)}')
    for r in extras:
        print(f'    at {r["start"]:#x}: {len(r["stream"])} chars: '
              f'{"".join(map(str, r["stream"][:60]))}...')
    report['extra_runs'] = [{'offset': hex(r['start']),
                             'len': len(r['stream']),
                             'stream': ''.join(map(str, r['stream']))}
                            for r in extras]

    # noita_dev.exe (informational -- different compile may change the pattern)
    dev = open(f'{gamedir}/noita_dev.exe', 'rb').read()
    druns = extract_runs(dev)
    dmatch = sum(1 for r in druns
                 if any(r['stream'] == exp for exp in streams.values()))
    print(f'noita_dev.exe: {len(druns)} decodable run(s), '
          f'{dmatch}/9 messages matched (informational)')
    report['dev_exe'] = {'runs': len(druns), 'matched': dmatch}

    # wak sweeps (completing the bounded Step E search)
    blob, entries = wak_unpack.read_index()
    kw = wak_keyword_sweep(blob, entries)
    print(f'data.wak keyword sweep: {len(kw)} hit(s)')
    for name, k in kw[:20]:
        print(f'    {k:12s} {name}')
    report['wak_keyword_hits'] = kw[:100]

    pngs, n_png, failed = wak_png_sweep(blob, entries)
    print(f'data.wak PNG eye-outline sweep: {n_png} PNGs scanned '
          f'({failed} unreadable), {len(pngs)} with exact outline matches')
    for name, k in pngs[:20]:
        print(f'    {k:4d}  {name}')
    report['wak_png_hits'] = pngs[:100]

    allok = all(gates)
    report['gates_pass'] = allok
    json.dump(report, open(f'{ROOT}/data/exesearch_report.json', 'w'), indent=1)
    print('\nwrote data/exesearch_report.json')
    if not allok:
        print('GATE FAILURE -> verdicts NOT trustworthy')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
