#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download ENCODE files listed in metadata/manifest/download_manifest.tsv.

Input:
- metadata/manifest/download_manifest.tsv

Expected columns:
- file_accession
- experiment_accession
- target_label
- biosample_term_name
- file_format
- output_type
- assembly
- download_url
- local_path
- download_flag

Features:
- Skip rows where download_flag != yes
- Skip existing files by default
- Create parent directories automatically
- Retry failed downloads
- Write download log to logs/download_encode_data.log
- Write per-file status table to logs/download_status.tsv

Recommended usage:
    cd ~/wangshuo/rna_dynamics
    python scripts/02_download/download_encode_data.py

Useful options:
    python scripts/02_download/download_encode_data.py --dry-run
    python scripts/02_download/download_encode_data.py --limit 3
    python scripts/02_download/download_encode_data.py --overwrite
    python scripts/02_download/download_encode_data.py --targets TARDBP
    python scripts/02_download/download_encode_data.py --biosamples HepG2
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests


CHUNK_SIZE = 1024 * 1024  # 1 MB
DEFAULT_TIMEOUT = 120
DEFAULT_RETRY = 3
DEFAULT_SLEEP = 1.5


STATUS_COLUMNS = [
    "timestamp",
    "file_accession",
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "download_url",
    "local_path",
    "status",
    "bytes_downloaded",
    "message",
]


def ensure_project_root() -> Path:
    script_path = Path(__file__).resolve()
    return script_path.parents[2]


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, log_file: Optional[Path] = None) -> None:
    line = f"{msg}"
    print(line, flush=True)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a", encoding="utf-8") as f:
            f.write(line + "\n")


def read_tsv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def append_tsv_row(path: Path, row: Dict[str, Any], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    file_exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        if not file_exists:
            writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in fieldnames})


def normalize_text(x: Any) -> str:
    if x is None:
        return ""
    return str(x).strip()


def normalize_lower(x: Any) -> str:
    return normalize_text(x).lower()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download ENCODE files from download_manifest.tsv"
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Path to download manifest TSV. Default: metadata/manifest/download_manifest.tsv",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N eligible rows",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing local files",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without downloading",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"HTTP timeout in seconds (default: {DEFAULT_TIMEOUT})",
    )
    parser.add_argument(
        "--retry",
        type=int,
        default=DEFAULT_RETRY,
        help=f"Retry count for failed downloads (default: {DEFAULT_RETRY})",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="Sleep seconds between files (default: 0.2)",
    )
    parser.add_argument(
        "--targets",
        nargs="*",
        default=None,
        help="Restrict to these target_label values",
    )
    parser.add_argument(
        "--biosamples",
        nargs="*",
        default=None,
        help="Restrict to these biosample_term_name values",
    )
    return parser.parse_args()


def row_is_eligible(
    row: Dict[str, str],
    target_filter: Optional[set],
    biosample_filter: Optional[set],
) -> bool:
    if normalize_lower(row.get("download_flag")) != "yes":
        return False

    target_label = normalize_lower(row.get("target_label"))
    biosample = normalize_lower(row.get("biosample_term_name"))

    if target_filter is not None and target_label not in target_filter:
        return False
    if biosample_filter is not None and biosample not in biosample_filter:
        return False

    return True


def prepare_rows(
    rows: List[Dict[str, str]],
    target_filter: Optional[set],
    biosample_filter: Optional[set],
    limit: Optional[int],
) -> List[Dict[str, str]]:
    eligible = [r for r in rows if row_is_eligible(r, target_filter, biosample_filter)]
    eligible = sorted(
        eligible,
        key=lambda x: (
            normalize_text(x.get("target_label")),
            normalize_text(x.get("biosample_term_name")),
            normalize_text(x.get("experiment_accession")),
            normalize_text(x.get("file_accession")),
        ),
    )
    if limit is not None:
        eligible = eligible[:limit]
    return eligible


def write_status(
    status_path: Path,
    row: Dict[str, str],
    status: str,
    bytes_downloaded: int,
    message: str,
) -> None:
    status_row = {
        "timestamp": now_str(),
        "file_accession": normalize_text(row.get("file_accession")),
        "experiment_accession": normalize_text(row.get("experiment_accession")),
        "target_label": normalize_text(row.get("target_label")),
        "biosample_term_name": normalize_text(row.get("biosample_term_name")),
        "download_url": normalize_text(row.get("download_url")),
        "local_path": normalize_text(row.get("local_path")),
        "status": status,
        "bytes_downloaded": bytes_downloaded,
        "message": message,
    }
    append_tsv_row(status_path, status_row, STATUS_COLUMNS)


def stream_download(
    url: str,
    out_path: Path,
    timeout: int,
) -> int:
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    bytes_written = 0

    with requests.get(url, stream=True, timeout=timeout) as resp:
        resp.raise_for_status()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with tmp_path.open("wb") as f:
            for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
                if not chunk:
                    continue
                f.write(chunk)
                bytes_written += len(chunk)

    tmp_path.replace(out_path)
    return bytes_written


def download_one(
    row: Dict[str, str],
    root: Path,
    timeout: int,
    retry: int,
    overwrite: bool,
    log_file: Path,
    status_path: Path,
) -> None:
    file_accession = normalize_text(row.get("file_accession"))
    url = normalize_text(row.get("download_url"))
    local_path_str = normalize_text(row.get("local_path"))

    if not url:
        msg = f"[ERROR] Missing download_url for {file_accession}"
        log(msg, log_file)
        write_status(status_path, row, "error", 0, "missing_download_url")
        return

    if not local_path_str:
        msg = f"[ERROR] Missing local_path for {file_accession}"
        log(msg, log_file)
        write_status(status_path, row, "error", 0, "missing_local_path")
        return

    out_path = root / local_path_str

    if out_path.exists() and not overwrite:
        msg = f"[SKIP] Exists: {out_path}"
        log(msg, log_file)
        write_status(status_path, row, "skipped_exists", out_path.stat().st_size, "file_exists")
        return

    last_err = None
    for attempt in range(1, retry + 1):
        try:
            msg = f"[INFO] Downloading ({attempt}/{retry}): {file_accession} -> {out_path}"
            log(msg, log_file)

            bytes_written = stream_download(url, out_path, timeout)

            msg = f"[DONE] Downloaded: {file_accession} | {bytes_written} bytes"
            log(msg, log_file)
            write_status(status_path, row, "downloaded", bytes_written, "ok")
            return

        except Exception as e:
            last_err = e
            msg = f"[WARN] Download failed ({attempt}/{retry}) for {file_accession}: {e}"
            log(msg, log_file)

            tmp_path = out_path.with_suffix(out_path.suffix + ".part")
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except Exception:
                    pass

            if attempt < retry:
                time.sleep(DEFAULT_SLEEP * attempt)

    msg = f"[ERROR] Failed after {retry} attempts: {file_accession} | {last_err}"
    log(msg, log_file)
    write_status(status_path, row, "error", 0, str(last_err))


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    manifest_path = (
        Path(args.manifest).resolve()
        if args.manifest
        else root / "metadata" / "manifest" / "download_manifest.tsv"
    )
    log_file = root / "logs" / "download_encode_data.log"
    status_path = root / "logs" / "download_status.tsv"

    rows = read_tsv_rows(manifest_path)

    target_filter = {normalize_lower(x) for x in args.targets} if args.targets else None
    biosample_filter = {normalize_lower(x) for x in args.biosamples} if args.biosamples else None

    eligible_rows = prepare_rows(
        rows=rows,
        target_filter=target_filter,
        biosample_filter=biosample_filter,
        limit=args.limit,
    )

    log(f"[INFO] Project root: {root}", log_file)
    log(f"[INFO] Manifest path: {manifest_path}", log_file)
    log(f"[INFO] Manifest rows total: {len(rows)}", log_file)
    log(f"[INFO] Eligible rows: {len(eligible_rows)}", log_file)
    log(f"[INFO] Overwrite: {args.overwrite}", log_file)
    log(f"[INFO] Dry run: {args.dry_run}", log_file)

    if args.targets:
        log(f"[INFO] Target filter: {', '.join(args.targets)}", log_file)
    if args.biosamples:
        log(f"[INFO] Biosample filter: {', '.join(args.biosamples)}", log_file)

    if not eligible_rows:
        log("[WARN] No eligible rows found. Nothing to do.", log_file)
        return 0

    log("", log_file)
    log("[INFO] First 5 eligible rows preview:", log_file)
    for row in eligible_rows[:5]:
        log(
            "  "
            + " | ".join(
                [
                    normalize_text(row.get("file_accession")),
                    normalize_text(row.get("experiment_accession")),
                    normalize_text(row.get("target_label")),
                    normalize_text(row.get("biosample_term_name")),
                    normalize_text(row.get("local_path")),
                ]
            ),
            log_file,
        )

    if args.dry_run:
        log("", log_file)
        log("[DONE] Dry run only. No files were downloaded.", log_file)
        return 0

    downloaded = 0
    skipped = 0
    failed = 0

    for i, row in enumerate(eligible_rows, start=1):
        file_accession = normalize_text(row.get("file_accession"))
        local_path_str = normalize_text(row.get("local_path"))
        out_path = root / local_path_str

        log("", log_file)
        log(f"[INFO] [{i}/{len(eligible_rows)}] Processing {file_accession}", log_file)

        before_exists = out_path.exists()
        before_size = out_path.stat().st_size if before_exists else 0

        download_one(
            row=row,
            root=root,
            timeout=args.timeout,
            retry=args.retry,
            overwrite=args.overwrite,
            log_file=log_file,
            status_path=status_path,
        )

        after_exists = out_path.exists()
        after_size = out_path.stat().st_size if after_exists else 0

        if before_exists and not args.overwrite:
            skipped += 1
        elif after_exists and after_size > 0:
            downloaded += 1
        else:
            failed += 1

        time.sleep(args.sleep)

    log("", log_file)
    log("[DONE] Download summary:", log_file)
    log(f"  Downloaded : {downloaded}", log_file)
    log(f"  Skipped    : {skipped}", log_file)
    log(f"  Failed     : {failed}", log_file)
    log(f"  Log file   : {log_file}", log_file)
    log(f"  Status TSV : {status_path}", log_file)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())