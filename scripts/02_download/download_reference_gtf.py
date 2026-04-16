#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download a reference GTF annotation file into the project structure.

Default target:
- GENCODE human Release 49 (GRCh38.p14)
- Basic gene annotation CHR
- File:
  gencode.v49.basic.annotation.gtf.gz

Output directory:
- raw/annotations/

Recommended usage:
    cd /data15/data15_5/junguang/wangshuo/rna_dynamics
    /usr/bin/python3 scripts/02_download/download_reference_gtf.py

Optional:
    /usr/bin/python3 scripts/02_download/download_reference_gtf.py --overwrite
    /usr/bin/python3 scripts/02_download/download_reference_gtf.py --url <custom_url>
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Optional

import requests


DEFAULT_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_49/gencode.v49.basic.annotation.gtf.gz"
)
DEFAULT_FILENAME = "gencode.v49.basic.annotation.gtf.gz"
CHUNK_SIZE = 1024 * 1024  # 1 MB
TIMEOUT = 120


def ensure_project_root() -> Path:
    script_path = Path(__file__).resolve()
    return script_path.parents[2]


def log(msg: str) -> None:
    print(msg, flush=True)


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def is_gzip_file(path: Path) -> bool:
    with path.open("rb") as f:
        magic = f.read(2)
    return magic == b"\x1f\x8b"


def download_file(url: str, out_path: Path, overwrite: bool = False) -> None:
    if out_path.exists() and not overwrite:
        log(f"[SKIP] File already exists: {out_path}")
        log(f"[INFO] SHA256: {sha256_of_file(out_path)}")
        log(f"[INFO] Gzip magic check: {'yes' if is_gzip_file(out_path) else 'no'}")
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")

    if tmp_path.exists():
        tmp_path.unlink()

    log(f"[INFO] Downloading: {url}")
    log(f"[INFO] Output path: {out_path}")

    bytes_written = 0
    with requests.get(url, stream=True, timeout=TIMEOUT) as resp:
        resp.raise_for_status()
        total_size = int(resp.headers.get("Content-Length", 0))

        with tmp_path.open("wb") as f:
            for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
                if not chunk:
                    continue
                f.write(chunk)
                bytes_written += len(chunk)

                if total_size > 0:
                    pct = bytes_written / total_size * 100
                    log(
                        f"[INFO] Downloaded {bytes_written:,} / {total_size:,} bytes "
                        f"({pct:.2f}%)"
                    )
                else:
                    log(f"[INFO] Downloaded {bytes_written:,} bytes")

    tmp_path.replace(out_path)

    log(f"[DONE] Download completed: {out_path}")
    log(f"[INFO] File size: {out_path.stat().st_size:,} bytes")
    log(f"[INFO] SHA256: {sha256_of_file(out_path)}")
    log(f"[INFO] Gzip magic check: {'yes' if is_gzip_file(out_path) else 'no'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download a reference GTF into the project raw/annotations directory."
    )
    parser.add_argument(
        "--url",
        default=DEFAULT_URL,
        help="GTF download URL",
    )
    parser.add_argument(
        "--filename",
        default=DEFAULT_FILENAME,
        help="Output filename under raw/annotations/",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing file",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    out_dir = root / "raw" / "annotations"
    out_path = out_dir / args.filename

    log(f"[INFO] Project root: {root}")
    log(f"[INFO] Target directory: {out_dir}")

    try:
        download_file(
            url=args.url,
            out_path=out_path,
            overwrite=args.overwrite,
        )
    except Exception as e:
        log(f"[ERROR] Download failed: {e}")
        return 1

    log("")
    log("[INFO] Suggested next step:")
    log(
        f"  /usr/bin/python3 {root}/scripts/04_analysis/01_annotate_peaks_with_gtf_regions.py "
        f"--gtf {out_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())