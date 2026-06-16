#!/usr/bin/env python3
"""
IGS Product Downloader
======================
Parses a RINEX 3 observation filename, derives year + DOY, then downloads
all IGS/MGEX product files needed by IGS.posEstimate() into the correct
directory alongside existing products.

Usage
-----
  python3 download_igs_products.py <obs_rnx_path> [--base <products_base_dir>] [--ac <AC>] [--dry-run]

Arguments
---------
  obs_rnx_path        Path to the RINEX 3 observation file (or just its filename).
                      The year and DOY are extracted from the embedded timestamp
                      (e.g. AJAC00FRA_S_20242000530_15M_01S_MO.rnx → year=2024, doy=200).

  --base DIR          Root products directory.  Products are stored in
                      <base>/<year>_<doy>/  to match the pDir() convention.
                      Default: the sibling directory of the obs file's parent
                      named "input_files" (or the same dir as the obs file).

  --ac AC             Analysis centre for orbit/clock/OSB products.
                      Choices: WUM (default), COD, GFZ, JAX
                      DCB: always CAS0OPSRAP; IONEX/SINEX: always IGS0OPS*.

  --dry-run           Print URLs and target paths without downloading.

Mirror priority
---------------
  1. IGN   (https://igs.ign.fr)          — public, no credentials
  2. BKG   (https://igs.bkg.bund.de)     — public, no credentials
  3. CDDIS (https://cddis.nasa.gov)      — NASA Earthdata login required

CDDIS credentials (never put these in code or on the command line):
  Add one line to ~/.netrc  (create the file if it doesn't exist):

      machine urs.earthdata.nasa.gov login YOUR_USERNAME password YOUR_PASSWORD

  Then lock it down:  chmod 600 ~/.netrc
  The script reads ~/.netrc automatically — no env vars, no flags needed.

Products downloaded
-------------------
  - {AC}0MGXFIN_{year}{doy}0000_01D_05M_ORB.SP3
  - {AC}0MGXFIN_{year}{doy}0000_01D_30S_CLK.CLK
  - {AC}0MGXFIN_{year}{doy}0000_01D_01D_OSB.BIA
  - CAS0OPSRAP_{year}{doy}0000_01D_01D_DCB.BIA
  - IGS0OPSFIN_{year}{doy}0000_01D_02H_GIM.INX
  - IGS0OPSSNX_{year}{doy}0000_01D_01D_SOL.SNX
  - BRDC00IGS_R_{year}{doy}0000_01D_MN.rnx
"""

import sys, os, re, gzip, shutil, datetime, argparse
from pathlib import Path

try:
    import requests
    from requests.auth import HTTPBasicAuth
    HAS_REQUESTS = True
except ImportError:
    import urllib.request
    HAS_REQUESTS = False

# ─── GPS week ─────────────────────────────────────────────────────────────────

GPS_EPOCH = datetime.date(1980, 1, 6)

def gps_week_dow(year: int, doy: int):
    d = datetime.date(year, 1, 1) + datetime.timedelta(days=doy - 1)
    delta = (d - GPS_EPOCH).days
    return delta // 7, delta % 7

# ─── Filename parsing ─────────────────────────────────────────────────────────

def parse_obs_filename(obs_path: str):
    """
    Extracts (year, doy) from a RINEX 3 filename like:
      AJAC00FRA_S_20242000530_15M_01S_MO.rnx
                   ^^^^---  year=2024
                       ^^^  doy=200
    """
    name = Path(obs_path).name
    m = re.search(r'_[A-Z_]_(\d{4})(\d{3})\d{4}_', name)
    if not m:
        # try bare timestamp (older style)
        m = re.search(r'(\d{4})(\d{3})', name)
    if not m:
        raise ValueError(f"Cannot parse year/DOY from filename: {name}")
    year, doy = int(m.group(1)), int(m.group(2))
    return year, doy

# ─── URL catalogue ────────────────────────────────────────────────────────────

def build_file_list(year: int, doy: int, ac: str):
    """
    Returns a list of (local_filename, [(mirror_url, auth_needed), ...]) tuples.
    Files are served compressed (.gz) on IGS mirrors; we decompress after download.
    """
    yy   = str(year)
    dd   = f"{doy:03d}"
    wk, dow = gps_week_dow(year, doy)
    wk_s = f"{wk:04d}"

    # centres
    AC = ac.upper()
    ac_prod = f"{AC}0MGXFIN"

    # --- IGS MGEX product mirrors for orbit/clock/bias products
    # IGN  https://igs.ign.fr/pub/igs/products/mgex/{week}/
    # BKG  https://igs.bkg.bund.de/root_ftp/MGEX/products/{week}/
    # CDDIS https://cddis.nasa.gov/archive/gnss/products/mgex/{week}/   (auth)

    def mgex_urls(filename):
        f = filename + ".gz"
        return [
            (f"https://igs.ign.fr/pub/igs/products/mgex/{wk_s}/{f}",    False),
            (f"https://igs.bkg.bund.de/root_ftp/MGEX/products/{wk_s}/{f}", False),
            (f"https://cddis.nasa.gov/archive/gnss/products/mgex/{wk_s}/{f}", True),
        ]

    def ionex_urls(filename):
        f = filename + ".gz"
        return [
            (f"https://igs.ign.fr/pub/igs/products/ionos/{yy}/{dd}/{f}", False),
            (f"https://igs.bkg.bund.de/root_ftp/IGS/products/ionos/{yy}/{dd}/{f}", False),
            (f"https://cddis.nasa.gov/archive/gnss/products/ionex/{yy}/{dd}/{f}", True),
        ]

    def nav_urls(filename):
        f = filename + ".gz"
        yy2 = yy[2:]    # 2-digit year for old-style path component
        return [
            # CDDIS is the canonical source for BRDC combined nav
            (f"https://cddis.nasa.gov/archive/gnss/data/daily/{yy}/{dd}/{yy2}p/{f}", True),
            # BKG uses slightly different subpath structure
            (f"https://igs.bkg.bund.de/root_ftp/IGS/data/daily/{yy}/{dd}/{yy2}p/{f}", False),
            # IGN
            (f"https://igs.ign.fr/pub/igs/data/{yy}/{dd}/{yy2}p/{f}",   False),
            # WHU as extra fallback
            (f"https://ftp.gnsswhu.cn/pub/gps/data/daily/{yy}/{dd}/{yy2}p/{f}", False),
        ]

    def sinex_urls(filename):
        f = filename + ".gz"
        return [
            (f"https://igs.ign.fr/pub/igs/products/{wk_s}/{f}",         False),
            (f"https://igs.bkg.bund.de/root_ftp/IGS/products/{wk_s}/{f}", False),
            (f"https://cddis.nasa.gov/archive/gnss/products/{wk_s}/{f}", True),
        ]

    orb_file  = f"{ac_prod}_{yy}{dd}0000_01D_05M_ORB.SP3"
    clk_file  = f"{ac_prod}_{yy}{dd}0000_01D_30S_CLK.CLK"
    osb_file  = f"{ac_prod}_{yy}{dd}0000_01D_01D_OSB.BIA"
    dcb_file  = f"CAS0OPSRAP_{yy}{dd}0000_01D_01D_DCB.BIA"
    gim_file  = f"IGS0OPSFIN_{yy}{dd}0000_01D_02H_GIM.INX"
    snx_file  = f"IGS0OPSSNX_{yy}{dd}0000_01D_01D_SOL.SNX"
    nav_file  = f"BRDC00IGS_R_{yy}{dd}0000_01D_MN.rnx"

    return [
        (orb_file,  mgex_urls(orb_file)),
        (clk_file,  mgex_urls(clk_file)),
        (osb_file,  mgex_urls(osb_file)),
        (dcb_file,  mgex_urls(dcb_file)),
        (gim_file,  ionex_urls(gim_file)),
        (snx_file,  sinex_urls(snx_file)),
        (nav_file,  nav_urls(nav_file)),
    ]

# ─── Download helpers ─────────────────────────────────────────────────────────

def _download_requests(url: str, dest: Path, auth=None):
    r = requests.get(url, timeout=60, stream=True, auth=auth,
                     allow_redirects=True, verify=True)
    r.raise_for_status()
    with open(dest, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1 << 16):
            f.write(chunk)

def _download_urllib(url: str, dest: Path, user=None, pwd=None):
    if user:
        mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        mgr.add_password(None, url, user, pwd)
        handler = urllib.request.HTTPBasicAuthHandler(mgr)
        opener = urllib.request.build_opener(handler)
    else:
        opener = urllib.request.build_opener()
    with opener.open(url, timeout=60) as resp, open(dest, 'wb') as f:
        shutil.copyfileobj(resp, f)

def decompress_gz(gz_path: Path, out_path: Path):
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()

def _netrc_creds(host: str):
    """Return (user, password) from ~/.netrc for host, or (None, None)."""
    netrc_path = Path.home() / ".netrc"
    if not netrc_path.exists():
        return None, None
    try:
        import netrc as _netrc
        n = _netrc.netrc(str(netrc_path))
        entry = n.authenticators(host)
        if entry:
            return entry[0], entry[2]  # (login, account, password)
    except Exception:
        pass
    return None, None

def try_download(url: str, dest_gz: Path, needs_auth: bool):
    user, pwd = None, None
    if needs_auth:
        user, pwd = _netrc_creds("urs.earthdata.nasa.gov")
        if not user:
            return False, "no CDDIS credentials in ~/.netrc — skipping"
    try:
        if HAS_REQUESTS:
            auth = HTTPBasicAuth(user, pwd) if user else None
            _download_requests(url, dest_gz, auth=auth)
        else:
            _download_urllib(url, dest_gz, user, pwd)
        return True, None
    except Exception as e:
        return False, str(e)

# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("obs_path", help="RINEX 3 observation file path or filename")
    ap.add_argument("--base",   default=None,
                    help="Root products directory (default: auto-detected from obs path)")
    ap.add_argument("--ac",    default="WUM",
                    choices=["WUM", "COD", "GFZ", "JAX"],
                    help="Analysis centre for orbit/clock/OSB (default: WUM)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print URLs and paths without downloading")
    args = ap.parse_args()

    # ── 1. Parse year / DOY from obs filename ─────────────────────────────────
    year, doy = parse_obs_filename(args.obs_path)
    wk, dow = gps_week_dow(year, doy)
    print(f"Observation file : {Path(args.obs_path).name}")
    print(f"Year={year}  DOY={doy:03d}  GPS week={wk}  DOW={dow}")

    # ── 2. Resolve output directory ───────────────────────────────────────────
    if args.base:
        base = Path(args.base)
    else:
        # Walk upward from the obs file looking for a known base directory name
        obs_p = Path(args.obs_path).resolve()
        candidates = [
            obs_p.parent.parent / "input_files",  # sibling of obs parent
            obs_p.parent,                          # same dir as obs
        ]
        base = next((c for c in candidates if c.is_dir()), obs_p.parent)

    out_dir = base / f"{year}_{doy:03d}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory : {out_dir}")
    print()

    # ── 3. Build file list ────────────────────────────────────────────────────
    file_list = build_file_list(year, doy, args.ac)

    # ── 4. Download ───────────────────────────────────────────────────────────
    ok_count = 0
    for local_name, mirror_urls in file_list:
        dest = out_dir / local_name
        if dest.exists():
            print(f"  [SKIP]  {local_name}  (already exists)")
            ok_count += 1
            continue

        if args.dry_run:
            print(f"  [DRY]   {local_name}")
            for url, auth in mirror_urls:
                print(f"          {url}")
            continue

        print(f"  [GET]   {local_name}")
        success = False
        for url, needs_auth in mirror_urls:
            dest_gz = dest.with_suffix(dest.suffix + ".gz")
            ok, err = try_download(url, dest_gz, needs_auth)
            if ok:
                if dest_gz.suffix == ".gz" and local_name.endswith(".gz") is False:
                    # decompress
                    try:
                        decompress_gz(dest_gz, dest)
                        print(f"          ✓  {url}")
                    except Exception as e:
                        dest_gz.unlink(missing_ok=True)
                        print(f"          ✗  decompress failed: {e}")
                        continue
                else:
                    dest_gz.rename(dest)
                    print(f"          ✓  {url}")
                success = True
                ok_count += 1
                break
            else:
                print(f"          ✗  {url}  ({err})")

        if not success:
            print(f"  [FAIL]  {local_name} — all mirrors failed")

    print()
    if not args.dry_run:
        print(f"Done: {ok_count}/{len(file_list)} files available in {out_dir}")
        print()
        print("── MainApp.java snippet ──────────────────────────────────────────")
        dd = f"{doy:03d}"
        yy = str(year)
        ac = args.ac.upper()
        ac_prod = f"{ac}0MGXFIN"
        print(f'  String year = "{year}", doy = "{dd}";')
        print(f'  String dir  = "{out_dir}/";')
        print(f'  String dcb_bias_path = dir + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";')
        print(f'  String osb_bias_path = dir + "{ac_prod}_"  + year + doy + "0000_01D_01D_OSB.BIA";')
        print(f'  String clock_path    = dir + "{ac_prod}_"  + year + doy + "0000_01D_30S_CLK.CLK";')
        print(f'  String orbit_path    = dir + "{ac_prod}_"  + year + doy + "0000_01D_05M_ORB.SP3";')
        print(f'  String ionex_path    = dir + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";')
        print(f'  String sinex_path    = dir + "IGS0OPSSNX_" + year + doy + "0000_01D_01D_SOL.SNX";')
        nav = f"BRDC00IGS_R_{yy}{dd}0000_01D_MN.rnx"
        print(f'  String nav_path      = dir + "{nav}";')

if __name__ == "__main__":
    main()
