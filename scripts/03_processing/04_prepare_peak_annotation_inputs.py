#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare annotation-ready inputs from standardized primary peak files.

Purpose
-------
Convert standardized BED10 peak files into annotation-ready BED6 files and
build a manifest for downstream peak annotation.

Typical downstream tools
------------------------
- bedtools intersect
- bedtools closest
- HOMER annotatePeaks.pl
- custom Python/R annotation scripts

Inputs
------
- processed/peaks/selected/primary_peak_files.tsv
- processed/peaks/standardized/standardized_peak_file_summary.tsv
- processed/peaks/standardized/<...>/*.standardized.bed

Optional input
--------------
- GTF/GFF annotation file (--gtf)

Outputs
-------
- processed/peaks/annotation_inputs/annotation_input_manifest.tsv
- processed/peaks/annotation_inputs/annotation_input_summary.tsv
- processed/peaks/annotation_inputs/bed6/<target>/<biosample>/<experiment>/<file>.annotation_ready.bed
- logs/prepare_peak_annotation_inputs.log

What this script does
---------------------
1. Read selected primary peak files and their standardized outputs
2. Convert BED10 standardized peaks to annotation-ready BED6
3. Rebuild peak names into stable unique IDs
4. Keep coordinate sorting
5. Summarize chromosome naming style
6. Optionally inspect GTF chromosome naming style and report compatibility

BED6 output columns
-------------------
1. chrom
2. start
3. end
4. peak_id
5. score
6. strand
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


MANIFEST_COLUMNS = [
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "selected_file_accession",
    "standardized_input_path",
    "annotation_ready_bed6_path",
    "n_peaks",
    "chrom_style",
    "has_chr_prefix",
    "gtf_path",
    "gtf_chrom_style",
    "gtf_has_chr_prefix",
    "chrom_compatible_with_gtf",
    "status",
    "note",
]

SUMMARY_COLUMNS = [
    "target_label",
    "biosample_term_name",
    "experiment_accession",
    "selected_file_accession",
    "n_peaks",
    "n_unique_chromosomes",
    "top_chromosomes",
    "chrom_style",
    "has_chr_prefix",
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
    tail = c[3:] if c.startswith("chr") else c
    if tail.isdigit():
        return (0, int(tail))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if tail in special:
        return (1, special[tail])
    return (2, tail)


def infer_chrom_style(chroms: Iterable[str]) -> Tuple[str, str]:
    chrom_list = [normalize_text(c) for c in chroms if normalize_text(c)]
    if not chrom_list:
        return "unknown", "unknown"

    has_chr = sum(1 for c in chrom_list if c.startswith("chr"))
    has_no_chr = sum(1 for c in chrom_list if not c.startswith("chr"))

    if has_chr > 0 and has_no_chr == 0:
        return "chr_prefixed", "yes"
    if has_no_chr > 0 and has_chr == 0:
        return "non_chr_prefixed", "no"
    return "mixed", "mixed"


def inspect_gtf_style(gtf_path: Path, max_lines: int = 50000) -> Tuple[str, str]:
    chroms: List[str] = []
    seen = 0
    with open_maybe_gzip(gtf_path) as f:
        for raw in f:
            if seen >= max_lines:
                break
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 1:
                continue
            chroms.append(fields[0])
            seen += 1
            if len(chroms) >= 200:
                break
    return infer_chrom_style(chroms)


def make_output_bed6_path(
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
        / "annotation_inputs"
        / "bed6"
        / safe_filename_part(target_label)
        / safe_filename_part(biosample_term_name)
        / safe_filename_part(experiment_accession)
        / f"{safe_filename_part(file_accession)}.annotation_ready.bed"
    )


def convert_standardized_bed10_to_bed6(
    input_path: Path,
    output_path: Path,
    target_label: str,
    biosample_term_name: str,
    experiment_accession: str,
    file_accession: str,
) -> Dict[str, Any]:
    rows: List[List[str]] = []
    chrom_counter: Counter[str] = Counter()
    peak_idx = 0

    with open_maybe_gzip(input_path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 6:
                continue

            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            score = fields[4]
            strand = fields[5] if fields[5] in {"+", "-", "."} else "."

            start_i = safe_int(start, None)
            end_i = safe_int(end, None)
            score_i = safe_int(score, 0)

            if start_i is None or end_i is None:
                continue
            if start_i < 0 or end_i <= start_i:
                continue

            # BED score convention is usually 0-1000; clip for safety
            if score_i is None:
                score_i = 0
            score_i = max(0, min(1000, score_i))

            peak_idx += 1
            peak_id = (
                f"{safe_filename_part(target_label)}"
                f"__{safe_filename_part(biosample_term_name)}"
                f"__{safe_filename_part(experiment_accession)}"
                f"__{safe_filename_part(file_accession)}"
                f"__peak{peak_idx:07d}"
            )

            rows.append([
                chrom,
                str(start_i),
                str(end_i),
                peak_id,
                str(score_i),
                strand,
            ])
            chrom_counter[chrom] += 1

    rows.sort(key=lambda r: (chrom_sort_key(r[0]), int(r[1]), int(r[2]), r[3]))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")

    chrom_style, has_chr_prefix = infer_chrom_style(chrom_counter.keys())
    top_chromosomes = ";".join([f"{k}:{v}" for k, v in chrom_counter.most_common(10)])

    return {
        "n_peaks": len(rows),
        "n_unique_chromosomes": len(chrom_counter),
        "top_chromosomes": top_chromosomes,
        "chrom_style": chrom_style,
        "has_chr_prefix": has_chr_prefix,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare annotation-ready BED6 files from standardized primary peaks."
    )
    parser.add_argument(
        "--primary-tsv",
        default=None,
        help="Override processed/peaks/selected/primary_peak_files.tsv",
    )
    parser.add_argument(
        "--standardized-summary-tsv",
        default=None,
        help="Override processed/peaks/standardized/standardized_peak_file_summary.tsv",
    )
    parser.add_argument(
        "--gtf",
        default=None,
        help="Optional GTF/GFF file for chromosome naming compatibility check",
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
    standardized_summary_tsv = (
        Path(args.standardized_summary_tsv).resolve()
        if args.standardized_summary_tsv
        else root / "processed" / "peaks" / "standardized" / "standardized_peak_file_summary.tsv"
    )
    gtf_path = Path(args.gtf).resolve() if args.gtf else None

    out_dir = root / "processed" / "peaks" / "annotation_inputs"
    manifest_out = out_dir / "annotation_input_manifest.tsv"
    summary_out = out_dir / "annotation_input_summary.tsv"
    log_path = root / "logs" / "prepare_peak_annotation_inputs.log"

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Primary TSV: {primary_tsv}", log_path)
    log(f"[INFO] Standardized summary TSV: {standardized_summary_tsv}", log_path)

    primary_rows = read_tsv_rows(primary_tsv)
    standardized_rows = read_tsv_rows(standardized_summary_tsv)

    standardized_idx = {
        normalize_text(r.get("selected_file_accession")): r
        for r in standardized_rows
        if normalize_text(r.get("selected_file_accession"))
    }

    if args.targets:
        target_filter = {x.strip().lower() for x in args.targets if x.strip()}
        primary_rows = [
            r for r in primary_rows
            if normalize_text(r.get("target_label")).lower() in target_filter
        ]

    if args.biosamples:
        biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()}
        primary_rows = [
            r for r in primary_rows
            if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter
        ]

    if args.limit is not None:
        primary_rows = primary_rows[:args.limit]

    gtf_style = "not_provided"
    gtf_has_chr_prefix = "not_provided"
    if gtf_path is not None:
        if not gtf_path.exists():
            log(f"[ERROR] GTF file not found: {gtf_path}", log_path)
            return 1
        gtf_style, gtf_has_chr_prefix = inspect_gtf_style(gtf_path)
        log(
            f"[INFO] GTF style: {gtf_style} | has_chr_prefix={gtf_has_chr_prefix} | path={gtf_path}",
            log_path,
        )

    manifest_rows: List[Dict[str, Any]] = []
    summary_rows: List[Dict[str, Any]] = []

    log(f"[INFO] Files to prepare: {len(primary_rows)}", log_path)

    for i, row in enumerate(primary_rows, start=1):
        experiment_accession = normalize_text(row.get("experiment_accession"))
        target_label = normalize_text(row.get("target_label"))
        biosample_term_name = normalize_text(row.get("biosample_term_name"))
        file_accession = normalize_text(row.get("selected_file_accession"))

        std_row = standardized_idx.get(file_accession)
        if std_row is None:
            manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": file_accession,
                "standardized_input_path": "",
                "annotation_ready_bed6_path": "",
                "n_peaks": 0,
                "chrom_style": "unknown",
                "has_chr_prefix": "unknown",
                "gtf_path": str(gtf_path) if gtf_path else "",
                "gtf_chrom_style": gtf_style,
                "gtf_has_chr_prefix": gtf_has_chr_prefix,
                "chrom_compatible_with_gtf": "unknown",
                "status": "missing_standardized_summary",
                "note": "selected_file_accession_not_found_in_standardized_summary",
            })
            continue

        standardized_input_rel = normalize_text(std_row.get("output_standardized_path"))
        standardized_input_path = root / standardized_input_rel

        annotation_ready_path = make_output_bed6_path(
            root=root,
            target_label=target_label,
            biosample_term_name=biosample_term_name,
            experiment_accession=experiment_accession,
            file_accession=file_accession,
        )

        log(
            f"[INFO] [{i}/{len(primary_rows)}] Preparing annotation input for "
            f"{file_accession} | {target_label} | {biosample_term_name} | {experiment_accession}",
            log_path,
        )

        if not standardized_input_path.exists():
            manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": file_accession,
                "standardized_input_path": standardized_input_rel,
                "annotation_ready_bed6_path": str(annotation_ready_path.relative_to(root)),
                "n_peaks": 0,
                "chrom_style": "unknown",
                "has_chr_prefix": "unknown",
                "gtf_path": str(gtf_path) if gtf_path else "",
                "gtf_chrom_style": gtf_style,
                "gtf_has_chr_prefix": gtf_has_chr_prefix,
                "chrom_compatible_with_gtf": "unknown",
                "status": "missing_standardized_input",
                "note": "standardized_input_file_not_found",
            })
            continue

        try:
            stats = convert_standardized_bed10_to_bed6(
                input_path=standardized_input_path,
                output_path=annotation_ready_path,
                target_label=target_label,
                biosample_term_name=biosample_term_name,
                experiment_accession=experiment_accession,
                file_accession=file_accession,
            )

            chrom_compatible = "not_checked"
            if gtf_path is not None:
                chrom_compatible = (
                    "yes" if stats["has_chr_prefix"] == gtf_has_chr_prefix else "no"
                )

            manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": file_accession,
                "standardized_input_path": standardized_input_rel,
                "annotation_ready_bed6_path": str(annotation_ready_path.relative_to(root)),
                "n_peaks": stats["n_peaks"],
                "chrom_style": stats["chrom_style"],
                "has_chr_prefix": stats["has_chr_prefix"],
                "gtf_path": str(gtf_path) if gtf_path else "",
                "gtf_chrom_style": gtf_style,
                "gtf_has_chr_prefix": gtf_has_chr_prefix,
                "chrom_compatible_with_gtf": chrom_compatible,
                "status": "ok",
                "note": "annotation_ready_bed6_created",
            })

            summary_rows.append({
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "experiment_accession": experiment_accession,
                "selected_file_accession": file_accession,
                "n_peaks": stats["n_peaks"],
                "n_unique_chromosomes": stats["n_unique_chromosomes"],
                "top_chromosomes": stats["top_chromosomes"],
                "chrom_style": stats["chrom_style"],
                "has_chr_prefix": stats["has_chr_prefix"],
                "status": "ok",
                "note": "annotation_ready_bed6_created",
            })

        except Exception as e:
            manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": file_accession,
                "standardized_input_path": standardized_input_rel,
                "annotation_ready_bed6_path": str(annotation_ready_path.relative_to(root)),
                "n_peaks": 0,
                "chrom_style": "unknown",
                "has_chr_prefix": "unknown",
                "gtf_path": str(gtf_path) if gtf_path else "",
                "gtf_chrom_style": gtf_style,
                "gtf_has_chr_prefix": gtf_has_chr_prefix,
                "chrom_compatible_with_gtf": "unknown",
                "status": "error",
                "note": str(e),
            })

    manifest_rows = sorted(
        manifest_rows,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
        ),
    )
    summary_rows = sorted(
        summary_rows,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
        ),
    )

    write_tsv(manifest_out, manifest_rows, MANIFEST_COLUMNS)
    write_tsv(summary_out, summary_rows, SUMMARY_COLUMNS)

    n_ok = sum(1 for r in manifest_rows if normalize_text(r.get("status")) == "ok")
    n_err = sum(1 for r in manifest_rows if normalize_text(r.get("status")) == "error")
    n_missing = sum(
        1
        for r in manifest_rows
        if normalize_text(r.get("status")) in {
            "missing_standardized_summary",
            "missing_standardized_input",
        }
    )

    log("", log_path)
    log("[DONE] Annotation input preparation completed.", log_path)
    log(f"[DONE] Manifest TSV: {manifest_out}", log_path)
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
                    f"n_peaks={normalize_text(r.get('n_peaks'))}",
                    normalize_text(r.get("chrom_style")),
                    normalize_text(r.get("status")),
                ]
            ),
            log_path,
        )

    return 0 if n_err == 0 and n_missing == 0 else 1


if __name__ == "__main__":
    sys.exit(main())