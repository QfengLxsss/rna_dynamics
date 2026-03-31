#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Inventory downloaded ENCODE peak BED files.

Inputs:
- metadata/manifest/download_manifest.tsv

Outputs:
- processed/peaks/inventory/peak_file_inventory.tsv
- processed/peaks/inventory/experiment_peak_summary.tsv
- logs/peak_inventory.log

Main tasks:
1. Check whether each manifest local_path exists
2. Record file size
3. Read first N non-empty, non-comment lines
4. Infer BED-like column count consistency
5. Summarize per experiment how many peak files were downloaded
6. Flag obviously suspicious files (missing/empty/inconsistent)

Important robustness updates:
- Detect gzip by magic bytes, not filename suffix alone
- Split columns on arbitrary whitespace, not only tabs
- Better BED-like status classification

Recommended usage:
    cd /data15/data15_5/junguang/wangshuo/rna_dynamics
    python scripts/03_processing/01_inventory_peak_files.py

Optional:
    python scripts/03_processing/01_inventory_peak_files.py --preview-lines 3
    python scripts/03_processing/01_inventory_peak_files.py --limit 5
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


INVENTORY_COLUMNS = [
    "file_accession",
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "file_format",
    "output_type",
    "assembly",
    "download_url",
    "local_path",
    "exists",
    "file_size_bytes",
    "preview_line_count",
    "first_data_line",
    "inferred_column_count",
    "column_count_consistent",
    "status",
    "note",
]

EXPERIMENT_SUMMARY_COLUMNS = [
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "n_manifest_files",
    "n_existing_files",
    "n_missing_files",
    "n_empty_files",
    "n_consistent_files",
    "observed_column_counts",
    "status",
    "note",
]


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, log_path: Optional[Path] = None) -> None:
    print(msg, flush=True)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(f"[{now_str()}] {msg}\n")


def ensure_project_root() -> Path:
    script_path = Path(__file__).resolve()
    return script_path.parents[2]


def read_tsv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def normalize_text(x: Any) -> str:
    if x is None:
        return ""
    return str(x).strip()


def open_maybe_gzip(path: Path):
    """
    Open plain text or gzip-compressed text file, regardless of filename suffix.
    Detect gzip by magic bytes, not by extension alone.
    """
    with path.open("rb") as f:
        magic = f.read(2)

    # gzip magic number: 1f 8b
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")

    return path.open("r", encoding="utf-8", errors="replace")


def inspect_bed_like_file(path: Path, preview_lines: int = 5) -> Tuple[int, str, str, str]:
    """
    Return:
    - preview_line_count
    - first_data_line
    - inferred_column_count
    - column_count_consistent

    Rules:
    - skip empty lines
    - skip comment/header-like lines
    - split on arbitrary whitespace, not only tabs
    """
    kept_lines: List[str] = []
    col_counts: List[int] = []

    with open_maybe_gzip(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            lower = stripped.lower()
            if (
                stripped.startswith("#")
                or lower.startswith("track")
                or lower.startswith("browser")
            ):
                continue

            fields = stripped.split()
            if not fields:
                continue

            kept_lines.append(stripped)
            col_counts.append(len(fields))

            if len(kept_lines) >= preview_lines:
                break

    if not kept_lines:
        return 0, "", "", ""

    first_data_line = kept_lines[0]
    unique_counts = sorted(set(col_counts))
    inferred_column_count = ";".join(str(x) for x in unique_counts)
    column_count_consistent = "yes" if len(unique_counts) == 1 else "no"

    return len(kept_lines), first_data_line, inferred_column_count, column_count_consistent


def classify_file_status(
    exists: str,
    file_size_bytes: int,
    preview_line_count: int,
    column_count_consistent: str,
    inferred_column_count: str,
) -> Tuple[str, str]:
    if exists != "yes":
        return "missing", "file_not_found"
    if file_size_bytes <= 0:
        return "empty", "zero_byte_file"
    if preview_line_count == 0:
        return "empty_or_no_data", "no_nonempty_data_lines_found"
    if column_count_consistent != "yes":
        return "inconsistent_columns", f"multiple_column_counts:{inferred_column_count}"

    try:
        first_ncol = int(inferred_column_count.split(";")[0])
    except Exception:
        return "manual_check", "cannot_parse_column_count"

    if first_ncol < 3:
        return "invalid_bed_like", f"too_few_columns:{first_ncol}"
    elif first_ncol in {3, 4, 5, 6, 9, 10, 12}:
        return "ok", f"standard_bed_like_ncol:{first_ncol}"
    else:
        return "ok", f"nonstandard_bed_like_ncol:{first_ncol}"


def build_inventory_row(
    manifest_row: Dict[str, str],
    root: Path,
    preview_lines: int,
) -> Dict[str, Any]:
    local_path_str = normalize_text(manifest_row.get("local_path"))
    abs_path = root / local_path_str

    exists = "yes" if abs_path.exists() else "no"
    file_size_bytes = abs_path.stat().st_size if abs_path.exists() else 0

    preview_line_count = 0
    first_data_line = ""
    inferred_column_count = ""
    column_count_consistent = ""

    if abs_path.exists() and file_size_bytes > 0:
        try:
            (
                preview_line_count,
                first_data_line,
                inferred_column_count,
                column_count_consistent,
            ) = inspect_bed_like_file(
                abs_path,
                preview_lines=preview_lines,
            )
        except Exception as e:
            return {
                "file_accession": normalize_text(manifest_row.get("file_accession")),
                "experiment_accession": normalize_text(manifest_row.get("experiment_accession")),
                "target_label": normalize_text(manifest_row.get("target_label")),
                "biosample_term_name": normalize_text(manifest_row.get("biosample_term_name")),
                "file_format": normalize_text(manifest_row.get("file_format")),
                "output_type": normalize_text(manifest_row.get("output_type")),
                "assembly": normalize_text(manifest_row.get("assembly")),
                "download_url": normalize_text(manifest_row.get("download_url")),
                "local_path": local_path_str,
                "exists": exists,
                "file_size_bytes": file_size_bytes,
                "preview_line_count": 0,
                "first_data_line": "",
                "inferred_column_count": "",
                "column_count_consistent": "",
                "status": "read_error",
                "note": str(e),
            }

    status, note = classify_file_status(
        exists=exists,
        file_size_bytes=file_size_bytes,
        preview_line_count=preview_line_count,
        column_count_consistent=column_count_consistent,
        inferred_column_count=inferred_column_count,
    )

    return {
        "file_accession": normalize_text(manifest_row.get("file_accession")),
        "experiment_accession": normalize_text(manifest_row.get("experiment_accession")),
        "target_label": normalize_text(manifest_row.get("target_label")),
        "biosample_term_name": normalize_text(manifest_row.get("biosample_term_name")),
        "file_format": normalize_text(manifest_row.get("file_format")),
        "output_type": normalize_text(manifest_row.get("output_type")),
        "assembly": normalize_text(manifest_row.get("assembly")),
        "download_url": normalize_text(manifest_row.get("download_url")),
        "local_path": local_path_str,
        "exists": exists,
        "file_size_bytes": file_size_bytes,
        "preview_line_count": preview_line_count,
        "first_data_line": first_data_line,
        "inferred_column_count": inferred_column_count,
        "column_count_consistent": column_count_consistent,
        "status": status,
        "note": note,
    }


def summarize_by_experiment(inventory_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in inventory_rows:
        exp = normalize_text(row.get("experiment_accession"))
        grouped[exp].append(row)

    summary_rows: List[Dict[str, Any]] = []

    for exp_acc, rows in grouped.items():
        target_label = normalize_text(rows[0].get("target_label"))
        biosample_term_name = normalize_text(rows[0].get("biosample_term_name"))

        n_manifest_files = len(rows)
        n_existing_files = sum(1 for r in rows if normalize_text(r.get("exists")) == "yes")
        n_missing_files = sum(1 for r in rows if normalize_text(r.get("status")) == "missing")
        n_empty_files = sum(
            1 for r in rows if normalize_text(r.get("status")) in {"empty", "empty_or_no_data"}
        )
        n_consistent_files = sum(
            1 for r in rows if normalize_text(r.get("column_count_consistent")) == "yes"
        )

        observed_column_counts = sorted(
            {
                normalize_text(r.get("inferred_column_count"))
                for r in rows
                if normalize_text(r.get("inferred_column_count"))
            }
        )
        observed_column_counts_text = ";".join(observed_column_counts)

        if n_missing_files > 0:
            status = "incomplete"
            note = "contains_missing_files"
        elif n_empty_files > 0:
            status = "manual_check"
            note = "contains_empty_files"
        else:
            inconsistent = any(normalize_text(r.get("column_count_consistent")) == "no" for r in rows)
            read_error = any(normalize_text(r.get("status")) == "read_error" for r in rows)
            invalid_bed = any(normalize_text(r.get("status")) == "invalid_bed_like" for r in rows)

            if read_error:
                status = "manual_check"
                note = "contains_read_error"
            elif invalid_bed:
                status = "manual_check"
                note = "contains_invalid_bed_like_file"
            elif inconsistent:
                status = "manual_check"
                note = "contains_inconsistent_column_counts"
            else:
                status = "ok"
                note = "all_files_present_and_bed_like"

        summary_rows.append(
            {
                "experiment_accession": exp_acc,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "n_manifest_files": n_manifest_files,
                "n_existing_files": n_existing_files,
                "n_missing_files": n_missing_files,
                "n_empty_files": n_empty_files,
                "n_consistent_files": n_consistent_files,
                "observed_column_counts": observed_column_counts_text,
                "status": status,
                "note": note,
            }
        )

    summary_rows = sorted(
        summary_rows,
        key=lambda x: (
            normalize_text(x.get("target_label")),
            normalize_text(x.get("biosample_term_name")),
            normalize_text(x.get("experiment_accession")),
        ),
    )
    return summary_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inventory downloaded ENCODE peak files."
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Path to download_manifest.tsv; default: metadata/manifest/download_manifest.tsv",
    )
    parser.add_argument(
        "--preview-lines",
        type=int,
        default=5,
        help="Number of non-empty data lines to preview per file (default: 5)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only inspect the first N manifest rows",
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


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    manifest_path = (
        Path(args.manifest).resolve()
        if args.manifest
        else root / "metadata" / "manifest" / "download_manifest.tsv"
    )

    out_dir = root / "processed" / "peaks" / "inventory"
    inventory_path = out_dir / "peak_file_inventory.tsv"
    summary_path = out_dir / "experiment_peak_summary.tsv"
    log_path = root / "logs" / "peak_inventory.log"

    manifest_rows = read_tsv_rows(manifest_path)

    if args.targets:
        target_filter = {x.strip().lower() for x in args.targets if x.strip()}
        manifest_rows = [
            r for r in manifest_rows
            if normalize_text(r.get("target_label")).lower() in target_filter
        ]

    if args.biosamples:
        biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()}
        manifest_rows = [
            r for r in manifest_rows
            if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter
        ]

    if args.limit is not None:
        manifest_rows = manifest_rows[:args.limit]

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Manifest path: {manifest_path}", log_path)
    log(f"[INFO] Manifest rows to inspect: {len(manifest_rows)}", log_path)
    log(f"[INFO] Preview lines per file: {args.preview_lines}", log_path)

    inventory_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(manifest_rows, start=1):
        file_accession = normalize_text(row.get("file_accession"))
        log(f"[INFO] [{i}/{len(manifest_rows)}] Inspecting {file_accession}", log_path)
        inv = build_inventory_row(
            manifest_row=row,
            root=root,
            preview_lines=args.preview_lines,
        )
        inventory_rows.append(inv)

    inventory_rows = sorted(
        inventory_rows,
        key=lambda x: (
            normalize_text(x.get("target_label")),
            normalize_text(x.get("biosample_term_name")),
            normalize_text(x.get("experiment_accession")),
            normalize_text(x.get("file_accession")),
        ),
    )

    summary_rows = summarize_by_experiment(inventory_rows)

    write_tsv(inventory_path, inventory_rows, INVENTORY_COLUMNS)
    write_tsv(summary_path, summary_rows, EXPERIMENT_SUMMARY_COLUMNS)

    n_ok = sum(1 for r in inventory_rows if normalize_text(r.get("status")) == "ok")
    n_missing = sum(1 for r in inventory_rows if normalize_text(r.get("status")) == "missing")
    n_empty = sum(1 for r in inventory_rows if normalize_text(r.get("status")) in {"empty", "empty_or_no_data"})
    n_inconsistent = sum(1 for r in inventory_rows if normalize_text(r.get("status")) == "inconsistent_columns")
    n_invalid = sum(1 for r in inventory_rows if normalize_text(r.get("status")) == "invalid_bed_like")
    n_read_error = sum(1 for r in inventory_rows if normalize_text(r.get("status")) == "read_error")

    log("", log_path)
    log("[DONE] Peak inventory completed.", log_path)
    log(f"[DONE] Inventory table: {inventory_path}", log_path)
    log(f"[DONE] Experiment summary: {summary_path}", log_path)
    log("", log_path)
    log("[INFO] File-level summary:", log_path)
    log(f"  OK files               : {n_ok}", log_path)
    log(f"  Missing files          : {n_missing}", log_path)
    log(f"  Empty/no-data files    : {n_empty}", log_path)
    log(f"  Inconsistent columns   : {n_inconsistent}", log_path)
    log(f"  Invalid BED-like files : {n_invalid}", log_path)
    log(f"  Read errors            : {n_read_error}", log_path)
    log("", log_path)
    log("[INFO] Experiment-level summary:", log_path)

    for row in summary_rows:
        log(
            "  "
            + " | ".join(
                [
                    normalize_text(row.get("experiment_accession")),
                    normalize_text(row.get("target_label")),
                    normalize_text(row.get("biosample_term_name")),
                    f"n={normalize_text(row.get('n_manifest_files'))}",
                    normalize_text(row.get("observed_column_counts")),
                    normalize_text(row.get("status")),
                ]
            ),
            log_path,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())