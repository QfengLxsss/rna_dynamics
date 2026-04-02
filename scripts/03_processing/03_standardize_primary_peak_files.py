#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Standardize selected primary ENCODE peak files into a unified narrowPeak-like BED10 format.

Input:
- processed/peaks/selected/primary_peak_files.tsv

Outputs:
- processed/peaks/standardized/<target>/<biosample>/<experiment>/<file>.standardized.bed
- processed/peaks/standardized/standardized_peak_file_summary.tsv
- logs/standardize_primary_peak_files.log

Standard output columns:
1. chrom
2. start
3. end
4. name
5. score
6. strand
7. signalValue
8. pValue
9. qValue
10. peak

Rules:
- keep non-empty, non-comment lines only
- split on arbitrary whitespace
- require at least 3 columns
- if <10 columns, pad to BED10-compatible defaults
- if >10 columns, keep first 10 columns
- enforce integer start/end
- require end > start and start >= 0
- de-duplicate by (chrom, start, end, name)
- sort by chrom, start, end, name
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


SUMMARY_COLUMNS = [
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "selected_file_accession",
    "input_local_path",
    "output_standardized_path",
    "input_peak_count_expected",
    "raw_data_lines_seen",
    "written_peak_count",
    "dropped_too_few_columns",
    "dropped_bad_coordinates",
    "dropped_end_le_start",
    "dropped_negative_start",
    "dropped_duplicate_rows",
    "status",
    "note",
]

STANDARD_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak",
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


def safe_int(x: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(str(x).strip())
    except Exception:
        return default


def safe_filename_part(text: str) -> str:
    s = normalize_text(text)
    if not s:
        return "NA"
    out = []
    for ch in s:
        if ch.isalnum() or ch in {"-", "_", "."}:
            out.append(ch)
        else:
            out.append("_")
    x = "".join(out)
    while "__" in x:
        x = x.replace("__", "_")
    return x.strip("_") or "NA"


def open_maybe_gzip(path: Path):
    with path.open("rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def chrom_sort_key(chrom: str) -> Tuple[int, Any]:
    c = chrom.strip()
    if c.startswith("chr"):
        tail = c[3:]
    else:
        tail = c

    # Common human chromosome ordering
    if tail.isdigit():
        return (0, int(tail))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if tail in special:
        return (1, special[tail])
    return (2, tail)


def normalize_bed10(fields: List[str], row_index: int) -> Optional[List[str]]:
    """
    Convert an arbitrary BED-like row to BED10.
    Return None if fewer than 3 columns.
    """
    if len(fields) < 3:
        return None

    fields = list(fields[:10])

    while len(fields) < 10:
        if len(fields) == 3:
            fields.append(f"peak_{row_index}")
        elif len(fields) == 4:
            fields.append("0")      # score
        elif len(fields) == 5:
            fields.append(".")      # strand
        elif len(fields) == 6:
            fields.append("0")      # signalValue
        elif len(fields) == 7:
            fields.append("-1")     # pValue
        elif len(fields) == 8:
            fields.append("-1")     # qValue
        elif len(fields) == 9:
            fields.append("-1")     # peak
        else:
            fields.append(".")

    chrom = fields[0]
    start = fields[1]
    end = fields[2]
    name = fields[3] if fields[3] not in {"", "."} else f"peak_{row_index}"
    score = fields[4] if fields[4] != "" else "0"
    strand = fields[5] if fields[5] in {"+", "-", "."} else "."
    signalValue = fields[6] if fields[6] != "" else "0"
    pValue = fields[7] if fields[7] != "" else "-1"
    qValue = fields[8] if fields[8] != "" else "-1"
    peak = fields[9] if fields[9] != "" else "-1"

    return [chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak]


def make_output_path(
    root: Path,
    target_label: str,
    biosample_term_name: str,
    experiment_accession: str,
    file_accession: str,
) -> Path:
    return (
        root
        / "processed"
        / "peaks"
        / "standardized"
        / safe_filename_part(target_label)
        / safe_filename_part(biosample_term_name)
        / safe_filename_part(experiment_accession)
        / f"{safe_filename_part(file_accession)}.standardized.bed"
    )


def standardize_one_file(
    input_path: Path,
    output_path: Path,
    row_label: str,
) -> Dict[str, Any]:
    raw_data_lines_seen = 0
    dropped_too_few_columns = 0
    dropped_bad_coordinates = 0
    dropped_end_le_start = 0
    dropped_negative_start = 0
    dropped_duplicate_rows = 0

    seen_keys = set()
    standardized_rows: List[List[str]] = []

    with open_maybe_gzip(input_path) as f:
        for i, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            lower = line.lower()
            if line.startswith("#") or lower.startswith("track") or lower.startswith("browser"):
                continue

            raw_data_lines_seen += 1
            fields = line.split()

            bed10 = normalize_bed10(fields, row_index=raw_data_lines_seen)
            if bed10 is None:
                dropped_too_few_columns += 1
                continue

            chrom, start_s, end_s, name, score, strand, signalValue, pValue, qValue, peak = bed10

            start = safe_int(start_s, None)
            end = safe_int(end_s, None)
            if start is None or end is None:
                dropped_bad_coordinates += 1
                continue
            if start < 0:
                dropped_negative_start += 1
                continue
            if end <= start:
                dropped_end_le_start += 1
                continue

            key = (chrom, start, end, name)
            if key in seen_keys:
                dropped_duplicate_rows += 1
                continue
            seen_keys.add(key)

            # write normalized values back
            standardized_rows.append([
                chrom,
                str(start),
                str(end),
                name,
                str(score),
                strand if strand in {"+", "-", "."} else ".",
                str(signalValue),
                str(pValue),
                str(qValue),
                str(peak),
            ])

    standardized_rows.sort(
        key=lambda r: (
            chrom_sort_key(r[0]),
            int(r[1]),
            int(r[2]),
            r[3],
        )
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as out:
        for row in standardized_rows:
            out.write("\t".join(row) + "\n")

    return {
        "raw_data_lines_seen": raw_data_lines_seen,
        "written_peak_count": len(standardized_rows),
        "dropped_too_few_columns": dropped_too_few_columns,
        "dropped_bad_coordinates": dropped_bad_coordinates,
        "dropped_end_le_start": dropped_end_le_start,
        "dropped_negative_start": dropped_negative_start,
        "dropped_duplicate_rows": dropped_duplicate_rows,
        "status": "ok",
        "note": row_label,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Standardize selected primary peak files to BED10."
    )
    parser.add_argument(
        "--primary-tsv",
        default=None,
        help="Override processed/peaks/selected/primary_peak_files.tsv",
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
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process first N selected rows",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    primary_tsv = (
        Path(args.primary_tsv).resolve()
        if args.primary_tsv
        else root / "processed" / "peaks" / "selected" / "primary_peak_files.tsv"
    )

    summary_out = root / "processed" / "peaks" / "standardized" / "standardized_peak_file_summary.tsv"
    log_path = root / "logs" / "standardize_primary_peak_files.log"

    rows = read_tsv_rows(primary_tsv)

    if args.targets:
        target_filter = {x.strip().lower() for x in args.targets if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("target_label")).lower() in target_filter]

    if args.biosamples:
        biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter]

    if args.limit is not None:
        rows = rows[:args.limit]

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Primary TSV: {primary_tsv}", log_path)
    log(f"[INFO] Selected files to standardize: {len(rows)}", log_path)

    summary_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(rows, start=1):
        experiment_accession = normalize_text(row.get("experiment_accession"))
        target_label = normalize_text(row.get("target_label"))
        biosample_term_name = normalize_text(row.get("biosample_term_name"))
        file_accession = normalize_text(row.get("selected_file_accession"))
        input_local_path = normalize_text(row.get("selected_local_path"))
        expected_peak_count = normalize_text(row.get("selected_peak_count"))

        input_path = root / input_local_path
        output_path = make_output_path(
            root=root,
            target_label=target_label,
            biosample_term_name=biosample_term_name,
            experiment_accession=experiment_accession,
            file_accession=file_accession,
        )

        log(
            f"[INFO] [{i}/{len(rows)}] Standardizing {file_accession} | "
            f"{target_label} | {biosample_term_name} | {experiment_accession}",
            log_path,
        )

        if not input_path.exists():
            summary_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": file_accession,
                "input_local_path": input_local_path,
                "output_standardized_path": str(output_path.relative_to(root)),
                "input_peak_count_expected": expected_peak_count,
                "raw_data_lines_seen": 0,
                "written_peak_count": 0,
                "dropped_too_few_columns": 0,
                "dropped_bad_coordinates": 0,
                "dropped_end_le_start": 0,
                "dropped_negative_start": 0,
                "dropped_duplicate_rows": 0,
                "status": "missing_input",
                "note": "selected_local_path_not_found",
            })
            log(f"[WARN] Missing input file: {input_path}", log_path)
            continue

        try:
            stats = standardize_one_file(
                input_path=input_path,
                output_path=output_path,
                row_label=file_accession,
            )
            summary_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": file_accession,
                "input_local_path": input_local_path,
                "output_standardized_path": str(output_path.relative_to(root)),
                "input_peak_count_expected": expected_peak_count,
                **stats,
            })
        except Exception as e:
            summary_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": file_accession,
                "input_local_path": input_local_path,
                "output_standardized_path": str(output_path.relative_to(root)),
                "input_peak_count_expected": expected_peak_count,
                "raw_data_lines_seen": 0,
                "written_peak_count": 0,
                "dropped_too_few_columns": 0,
                "dropped_bad_coordinates": 0,
                "dropped_end_le_start": 0,
                "dropped_negative_start": 0,
                "dropped_duplicate_rows": 0,
                "status": "error",
                "note": str(e),
            })
            log(f"[ERROR] Failed to standardize {file_accession}: {e}", log_path)

    summary_rows = sorted(
        summary_rows,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
        ),
    )

    write_tsv(summary_out, summary_rows, SUMMARY_COLUMNS)

    n_ok = sum(1 for r in summary_rows if normalize_text(r.get("status")) == "ok")
    n_err = sum(1 for r in summary_rows if normalize_text(r.get("status")) == "error")
    n_missing = sum(1 for r in summary_rows if normalize_text(r.get("status")) == "missing_input")

    log("", log_path)
    log("[DONE] Standardization completed.", log_path)
    log(f"[DONE] Summary TSV: {summary_out}", log_path)
    log(f"[INFO] OK files      : {n_ok}", log_path)
    log(f"[INFO] Error files   : {n_err}", log_path)
    log(f"[INFO] Missing files : {n_missing}", log_path)

    for r in summary_rows:
        log(
            "  "
            + " | ".join(
                [
                    normalize_text(r.get("experiment_accession")),
                    normalize_text(r.get("target_label")),
                    normalize_text(r.get("biosample_term_name")),
                    normalize_text(r.get("selected_file_accession")),
                    f"written={normalize_text(r.get('written_peak_count'))}",
                    normalize_text(r.get("status")),
                ]
            ),
            log_path,
        )

    return 0 if n_err == 0 and n_missing == 0 else 1


if __name__ == "__main__":
    sys.exit(main())