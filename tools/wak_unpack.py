#!/usr/bin/env python3
"""Primary-source instrument: pure-Python reader for a Noita `data.wak` archive.

Used for the game-install / source-verification step (does the shipped game
contain a decoder, plaintext table, or the eye-message murals themselves?).
Reading the format directly avoids launching the Windows `noita.exe -wizard_unpak`
under Proton.

Archive format (decoded from the header of the shipped 2024 build; the leading
u32 is a zero magic, NOT a header size):
  u32 magic (=0)
  u32 num_entries
  u32 data_start   (offset where the first file blob begins; == end of index)
  u32 pad (=0)
  entry* : u32 offset, u32 length, u32 path_len, char[path_len] path

Self-check (the tool's control gate): every offset+length must lie within
[data_start, filesize] and every path must decode as ASCII; otherwise we abort,
because a wrong format guess would surface as an out-of-range offset immediately.
Validated on the 2024 build: 14745 entries, first entry `data/credits.txt`
(path_len 16), and `cat` of a known XML returns readable text.

WAK path is machine-specific (Steam library location varies per machine this
project has run on); `find_wak()` tries known locations, then a bounded glob
under common Steam install roots, before giving up. `read_index()` takes an
explicit path if you have one; this script's own subcommands always call it
with none (autodetect) -- edit CANDIDATES below if a new machine needs a new
path (cheaper than plumbing a flag through every caller).
Usage: wak_unpack.py list [substr] | cat <exact/path> | extract <substr> <dir>
"""
import glob
import os
import struct
import sys

CANDIDATES = [
    '/mnt/hdd1/SteamLibrary/steamapps/common/Noita/data/data.wak',
    os.path.expanduser(
        '~/.local/share/Steam/steamapps/common/Noita/data/data.wak'),
]


def find_wak():
    for p in CANDIDATES:
        if os.path.exists(p):
            return p
    for pattern in (
        os.path.expanduser('~/.local/share/Steam/steamapps/*/common/Noita/data/data.wak'),
        os.path.expanduser('~/.steam/**/steamapps/common/Noita/data/data.wak'),
        '/mnt/**/steamapps/common/Noita/data/data.wak',
    ):
        hits = glob.glob(pattern, recursive=True)
        if hits:
            return hits[0]
    raise SystemExit(
        'data.wak not found in any known location; add its path to '
        'CANDIDATES in tools/wak_unpack.py')


WAK = None  # resolved lazily so importing this module never touches disk


def read_index(path=None):
    if path is None:
        path = find_wak()
    with open(path, 'rb') as f:
        blob = f.read()
    fsz = len(blob)
    magic, num, data_start, _pad = struct.unpack_from('<IIII', blob, 0)
    if not (0 < num < 2000000 and 0 < data_start <= fsz):
        raise SystemExit(f'header sanity fail: num={num} data_start={data_start} fsz={fsz}')
    header_size = data_start
    entries = []
    p = 16
    for _ in range(num):
        off, length, nlen = struct.unpack_from('<III', blob, p)
        p += 12
        name = blob[p:p + nlen]
        p += nlen
        if not (header_size <= off <= fsz and off + length <= fsz):
            raise SystemExit(f'entry offset out of range at {name!r}: off={off} len={length}')
        try:
            name = name.decode('ascii')
        except UnicodeDecodeError:
            raise SystemExit(f'non-ascii path at offset {p}: {name!r}')
        entries.append((name, off, length))
    return blob, entries


def main():
    blob, entries = read_index()
    print(f'index OK: {len(entries)} files')
    if len(sys.argv) > 1 and sys.argv[1] == 'list':
        pat = sys.argv[2] if len(sys.argv) > 2 else ''
        for name, off, length in entries:
            if pat.lower() in name.lower():
                print(f'{length:9d}  {name}')
    elif len(sys.argv) > 1 and sys.argv[1] == 'cat':
        target = sys.argv[2]
        for name, off, length in entries:
            if name == target:
                sys.stdout.buffer.write(blob[off:off + length])
                return
        print(f'not found: {target}', file=sys.stderr)
    elif len(sys.argv) > 1 and sys.argv[1] == 'extract':
        import os
        pat, outdir = sys.argv[2], sys.argv[3]
        n = 0
        for name, off, length in entries:
            if pat.lower() in name.lower():
                dst = os.path.join(outdir, name)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                with open(dst, 'wb') as g:
                    g.write(blob[off:off + length])
                n += 1
        print(f'extracted {n} files matching {pat!r} to {outdir}')


if __name__ == '__main__':
    main()
