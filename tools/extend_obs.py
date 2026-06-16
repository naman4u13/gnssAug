#!/usr/bin/env python3
"""
extend_obs.py
=============
Downloads consecutive 15-minute RINEX 3 high-rate observation blocks and merges
them into a single extended file.

The IGS high-rate archive packages each station's full day into a single
.crx.tar file. This script streams that tar, stops as soon as the N target
files have been extracted, then decompresses (.crx.gz → .rnx) and merges.

Usage
-----
  python3 extend_obs.py <obs_rnx> [--count N] [--out <merged.rnx>] [--dry-run]

  obs_rnx   Path to the existing RINEX 3 observation file.
            Filename must follow the IGS RINEX 3 stream convention:
              {STATION}_S_{YEAR}{DOY}{HHMM}_{DUR}_{RATE}_{TYPE}.rnx
            e.g. AJAC00FRA_S_20242000530_15M_01S_MO.rnx

  --count N  Additional blocks to fetch (default: 6 = 90 min).

  --out      Output merged filename.

  --dry-run  Show what would be done without downloading.

Requires: pip install hatanaka requests
"""

import sys, re, gzip, shutil, datetime, argparse, tarfile
from pathlib import Path

try:
    import requests
    from requests.auth import HTTPBasicAuth
except ImportError:
    sys.exit("requests not installed — run: pip install requests")

try:
    import hatanaka
except ImportError:
    sys.exit("hatanaka not installed — run: pip install hatanaka")

# ── filename parsing ───────────────────────────────────────────────────────────

_OBS_RE = re.compile(
    r'^(?P<station>[A-Z0-9]{9})_(?P<src>[A-Z])_'
    r'(?P<year>\d{4})(?P<doy>\d{3})(?P<hh>\d{2})(?P<mm>\d{2})_'
    r'(?P<dur>\d{2}M)_(?P<rate>\d{2}S)_(?P<ftype>[A-Z]{2})\.rnx$'
)

def parse_obs_name(filename):
    m = _OBS_RE.match(Path(filename).name)
    if not m:
        raise ValueError(
            f"Filename doesn't match IGS RINEX 3 stream convention:\n  {filename}\n"
            "Expected: SSSSSSSSS_S_YYYYDDDHHMM_DDM_RRS_TT.rnx"
        )
    return m.groupdict()


def next_filenames(info, count):
    """Return list of (year, doy, rnx_filename) for `count` successive blocks."""
    dur_min = int(info['dur'][:-1])
    date = datetime.date(int(info['year']), 1, 1) + datetime.timedelta(
        days=int(info['doy']) - 1)
    hh, mm = int(info['hh']), int(info['mm'])

    results = []
    for _ in range(count):
        mm += dur_min
        if mm >= 60:
            hh += mm // 60
            mm %= 60
        if hh >= 24:
            date += datetime.timedelta(days=hh // 24)
            hh %= 24
        doy = date.timetuple().tm_yday
        fname = (
            f"{info['station']}_{info['src']}_{date.year}{doy:03d}"
            f"{hh:02d}{mm:02d}_{info['dur']}_{info['rate']}_{info['ftype']}.rnx"
        )
        results.append((date.year, doy, fname))
    return results


# ── archive URL ────────────────────────────────────────────────────────────────

def full_day_tar_url(station, src, year, doy, rate, ftype):
    """
    CDDIS full-day .crx.tar for a station.
    e.g. AJAC00FRA_S_20242000000_01D_01S_MO.crx.tar
    """
    yy = str(year)[2:]
    dd = f"{doy:03d}"
    tar_name = f"{station}_{src}_{year}{dd}0000_01D_{rate}_{ftype}.crx.tar"
    url = f"https://cddis.nasa.gov/archive/gnss/data/highrate/{year}/{dd}/{tar_name}"
    return url, tar_name


# ── auth ───────────────────────────────────────────────────────────────────────

def _netrc_creds(host):
    try:
        import netrc as _netrc
        n = _netrc.netrc()
        e = n.authenticators(host)
        if e:
            return e[0], e[2]
    except Exception:
        pass
    return None, None


# ── streaming tar extraction ───────────────────────────────────────────────────

def stream_extract(tar_url, target_crx_gz_names, out_dir, auth=None):
    """
    Stream the full-day .crx.tar from CDDIS, extracting only the files whose
    basenames are in target_crx_gz_names. Stops as soon as all targets found.

    Returns list of (rnx_filename, rnx_bytes) tuples.
    """
    remaining = set(target_crx_gz_names)
    extracted = []
    bytes_read = 0

    with requests.get(tar_url, auth=auth, stream=True,
                      timeout=120, allow_redirects=True) as r:
        r.raise_for_status()
        with tarfile.open(fileobj=r.raw, mode='r|') as tar:
            for member in tar:
                bytes_read += member.size
                basename = Path(member.name).name
                if basename in remaining:
                    f_obj = tar.extractfile(member)
                    if f_obj:
                        crx_gz_data = f_obj.read()
                        crx_data    = gzip.decompress(crx_gz_data)
                        rnx_bytes   = hatanaka.decompress(crx_data)
                        rnx_name    = basename.replace('.crx.gz', '.rnx')
                        extracted.append((rnx_name, rnx_bytes))
                        remaining.discard(basename)
                        print(f"           ✓  extracted {rnx_name}  "
                              f"({len(rnx_bytes)/1024:.0f} KB)  "
                              f"[{bytes_read/1e6:.1f} MB streamed so far]")
                if not remaining:
                    break   # all targets found — stop streaming

    return extracted, bytes_read


# ── RINEX merge ────────────────────────────────────────────────────────────────

def after_header_offset(path):
    with open(path, 'rb') as f:
        for line in f:
            if b'END OF HEADER' in line:
                return f.tell()
    raise ValueError(f"No END OF HEADER found in {path}")


def merge_rinex(base, extras, out_path):
    with open(out_path, 'wb') as out:
        with open(base, 'rb') as f:
            shutil.copyfileobj(f, out)
        for p in extras:
            offset = after_header_offset(p)
            with open(p, 'rb') as f:
                f.seek(offset)
                shutil.copyfileobj(f, out)


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("obs_rnx", help="Existing RINEX 3 high-rate obs file")
    ap.add_argument("--count", type=int, default=6,
                    help="Additional blocks to fetch (default: 6 → 90 min)")
    ap.add_argument("--out", default=None, help="Output merged file path")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    obs_path = Path(args.obs_rnx).resolve()
    if not obs_path.exists():
        sys.exit(f"File not found: {obs_path}")

    info = parse_obs_name(obs_path.name)
    dur_min = int(info['dur'][:-1])

    print(f"Base file  : {obs_path.name}")
    print(f"Station    : {info['station']}")
    print(f"Start      : {info['year']}-DOY{info['doy']}  {info['hh']}:{info['mm']} UTC")
    print(f"Fetching   : {args.count} × {info['dur']} = {args.count * dur_min} min")
    print()

    out_dir  = obs_path.parent
    default_out = obs_path.name.replace('.rnx', f'_plus{args.count}x{info["dur"]}.rnx')
    out_path = Path(args.out).resolve() if args.out else out_dir / default_out

    blocks = next_filenames(info, args.count)

    # Group blocks by (year, doy) — normally all same day but handle midnight rollover
    from itertools import groupby
    from operator import itemgetter

    day_groups = {}
    for year, doy, fname in blocks:
        day_groups.setdefault((year, doy), []).append(fname)

    user_creds = _netrc_creds("urs.earthdata.nasa.gov")
    if not user_creds[0]:
        sys.exit("No CDDIS credentials found in ~/.netrc — cannot download high-rate data.\n"
                 "Add:  machine urs.earthdata.nasa.gov login <user> password <pass>")

    auth = HTTPBasicAuth(*user_creds)

    all_downloaded = []   # ordered list of Path objects

    for (year, doy), day_fnames in day_groups.items():
        tar_url, tar_name = full_day_tar_url(
            info['station'], info['src'], year, doy, info['rate'], info['ftype'])

        # Check which files already exist
        needed, already = [], []
        for fname in day_fnames:
            dest = out_dir / fname
            if dest.exists():
                print(f"  [SKIP]  {fname}")
                already.append(dest)
            else:
                needed.append(fname)

        all_downloaded.extend(already)

        if not needed:
            continue

        if args.dry_run:
            print(f"  [DRY]   would stream {tar_url}")
            for f in needed:
                print(f"           extract {f}")
            continue

        target_crx_gz = {f.replace('.rnx', '.crx.gz') for f in needed}
        print(f"  [TAR]   streaming {tar_name}")
        print(f"           extracting {len(needed)} of {(24*60)//dur_min} files in archive")

        try:
            extracted, bytes_read = stream_extract(tar_url, target_crx_gz, out_dir, auth=auth)
        except Exception as e:
            print(f"  [FAIL]  {e}")
            continue

        # Save extracted files in the correct order
        extracted_map = {name: data for name, data in extracted}
        for fname in needed:
            if fname in extracted_map:
                dest = out_dir / fname
                dest.write_bytes(extracted_map[fname])
                all_downloaded.append(dest)
            else:
                print(f"  [MISS]  {fname} not found in tar")

        print(f"           {bytes_read/1e6:.1f} MB streamed total")

    print()
    if args.dry_run:
        return

    if not all_downloaded:
        print("No files available — nothing to merge.")
        return

    # Sort by timestamp embedded in filename before merging
    all_downloaded.sort(key=lambda p: p.name)

    merge_rinex(obs_path, all_downloaded, out_path)
    total_min = (1 + len(all_downloaded)) * dur_min
    print(f"Total coverage : {total_min} min  ({total_min / 60:.2f} h)")
    print(f"Merged file    : {out_path}")


if __name__ == "__main__":
    main()
