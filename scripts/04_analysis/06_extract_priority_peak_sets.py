#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract priority comparison peak sets for downstream overlap/specific/shared analysis.

Inputs
------
- results/tables/exploratory_priority_comparisons.tsv
- processed/peaks/standardized/standardized_peak_file_summary.tsv
- processed/peaks/annotations/annotation_file_manifest.tsv

Outputs
-------
- processed/peaks/comparison_sets/comparison_peak_set_manifest.tsv
- processed/peaks/comparison_sets/pairs/<comparison_id>/left.standardized.bed
- processed/peaks/comparison_sets/pairs/<comparison_id>/right.standardized.bed
- processed/peaks/comparison_sets/pairs/<comparison_id>/left.annotations.tsv
- processed/peaks/comparison_sets/pairs/<comparison_id>/right.annotations.tsv
- logs/extract_priority_peak_sets.log

Notes
-----
- By default, only sample-scope rows are materialized because aggregated rows do not
  correspond to a single unique file pair.
- Files are symlinked when possible; if symlink fails, they are copied.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


OUTPUT_COLUMNS = [
    "priority_rank",
    "priority_level",
    "comparison_scope",
    "comparison_type",
    "comparison_id",
    "left_target_label",
    "left_biosample_term_name",
    "left_experiment_accession",
    "left_selected_file_accession",
    "right_target_label",
    "right_biosample_term_name",
    "right_experiment_accession",
    "right_selected_file_accession",
    "left_standardized_input_path",
    "right_standardized_input_path",
    "left_annotation_tsv_path",
    "right_annotation_tsv_path",
    "left_materialized_bed_path",
    "right_materialized_bed_path",
    "left_materialized_annotation_path",
    "right_materialized_annotation_path",
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
        return []
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


def build_index(rows: List[Dict[str, str]], key: str) -> Dict[str, Dict[str, str]]:
    idx: Dict[str, Dict[str, str]] = {}
    for row in rows:
        k = normalize_text(row.get(key))
        if k:
            idx[k] = row
    return idx


def materialize_file(src: Path, dst: Path, log_path: Optional[Path] = None) -> str:
    """
    Prefer symlink; fallback to copy.
    Returns one of: symlinked, copied, exists
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists():
        return "exists"

    try:
        os.symlink(src, dst)
        return "symlinked"
    except Exception:
        shutil.copy2(src, dst)
        return "copied"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract materialized peak set inputs for priority comparisons."
    )
    parser.add_argument(
        "--priority-tsv",
        default=None,
        help="Override results/tables/exploratory_priority_comparisons.tsv",
    )
    parser.add_argument(
        "--standardized-summary",
        default=None,
        help="Override processed/peaks/standardized/standardized_peak_file_summary.tsv",
    )
    parser.add_argument(
        "--annotation-file-manifest",
        default=None,
        help="Override processed/peaks/annotations/annotation_file_manifest.tsv",
    )
    parser.add_argument(
        "--priority-levels",
        nargs="*",
        default=["priority_1"],
        help="Priority levels to extract (default: priority_1)",
    )
    parser.add_argument(
        "--include-aggregated",
        action="store_true",
        help="Include aggregated rows if possible (default: off; these are normally skipped).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process first N eligible rows",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    priority_tsv = (
        Path(args.priority_tsv).resolve()
        if args.priority_tsv
        else root / "results" / "tables" / "exploratory_priority_comparisons.tsv"
    )
    standardized_summary = (
        Path(args.standardized_summary).resolve()
        if args.standardized_summary
        else root / "processed" / "peaks" / "standardized" / "standardized_peak_file_summary.tsv"
    )
    annotation_file_manifest = (
        Path(args.annotation_file_manifest).resolve()
        if args.annotation_file_manifest
        else root / "processed" / "peaks" / "annotations" / "annotation_file_manifest.tsv"
    )

    out_dir = root / "processed" / "peaks" / "comparison_sets"
    manifest_out = out_dir / "comparison_peak_set_manifest.tsv"
    log_path = root / "logs" / "extract_priority_peak_sets.log"

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Priority TSV: {priority_tsv}", log_path)
    log(f"[INFO] Standardized summary TSV: {standardized_summary}", log_path)
    log(f"[INFO] Annotation file manifest TSV: {annotation_file_manifest}", log_path)

    priority_rows = read_tsv_rows(priority_tsv)
    standardized_rows = read_tsv_rows(standardized_summary)
    annotation_rows = read_tsv_rows(annotation_file_manifest)

    std_idx = build_index(standardized_rows, "selected_file_accession")
    ann_idx = build_index(annotation_rows, "selected_file_accession")

    allowed_levels = {normalize_text(x) for x in args.priority_levels if normalize_text(x)}

    eligible_rows: List[Dict[str, str]] = []
    for row in priority_rows:
        plevel = normalize_text(row.get("priority_level"))
        scope = normalize_text(row.get("comparison_scope"))

        if plevel not in allowed_levels:
            continue

        if scope == "aggregated" and not args.include_aggregated:
            continue

        if scope != "sample":
            # Aggregated rows normally do not have a unique left/right file pair
            if not args.include_aggregated:
                continue

        eligible_rows.append(row)

    eligible_rows = sorted(
        eligible_rows,
        key=lambda r: (
            normalize_text(r.get("priority_level")),
            normalize_text(r.get("comparison_scope")),
            normalize_text(r.get("priority_rank")),
        ),
    )

    if args.limit is not None:
        eligible_rows = eligible_rows[:args.limit]

    log(f"[INFO] Eligible priority rows: {len(eligible_rows)}", log_path)

    output_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(eligible_rows, start=1):
        priority_rank = normalize_text(row.get("priority_rank"))
        priority_level = normalize_text(row.get("priority_level"))
        comparison_scope = normalize_text(row.get("comparison_scope"))
        comparison_type = normalize_text(row.get("comparison_type"))
        comparison_id = normalize_text(row.get("comparison_id"))

        left_target = normalize_text(row.get("left_target_label"))
        left_biosample = normalize_text(row.get("left_biosample_term_name"))
        left_exp = normalize_text(row.get("left_experiment_accession"))
        left_file = normalize_text(row.get("left_selected_file_accession"))

        right_target = normalize_text(row.get("right_target_label"))
        right_biosample = normalize_text(row.get("right_biosample_term_name"))
        right_exp = normalize_text(row.get("right_experiment_accession"))
        right_file = normalize_text(row.get("right_selected_file_accession"))

        log(
            f"[INFO] [{i}/{len(eligible_rows)}] Extracting {priority_level} | {comparison_scope} | {comparison_id}",
            log_path,
        )

        if comparison_scope != "sample":
            output_rows.append({
                "priority_rank": priority_rank,
                "priority_level": priority_level,
                "comparison_scope": comparison_scope,
                "comparison_type": comparison_type,
                "comparison_id": comparison_id,
                "left_target_label": left_target,
                "left_biosample_term_name": left_biosample,
                "left_experiment_accession": left_exp,
                "left_selected_file_accession": left_file,
                "right_target_label": right_target,
                "right_biosample_term_name": right_biosample,
                "right_experiment_accession": right_exp,
                "right_selected_file_accession": right_file,
                "left_standardized_input_path": "",
                "right_standardized_input_path": "",
                "left_annotation_tsv_path": "",
                "right_annotation_tsv_path": "",
                "left_materialized_bed_path": "",
                "right_materialized_bed_path": "",
                "left_materialized_annotation_path": "",
                "right_materialized_annotation_path": "",
                "status": "skipped_non_sample_scope",
                "note": "only_sample_scope_rows_are_materialized_into_explicit_file_pairs",
            })
            continue

        left_std = std_idx.get(left_file)
        right_std = std_idx.get(right_file)
        left_ann = ann_idx.get(left_file)
        right_ann = ann_idx.get(right_file)

        if left_std is None or right_std is None:
            output_rows.append({
                "priority_rank": priority_rank,
                "priority_level": priority_level,
                "comparison_scope": comparison_scope,
                "comparison_type": comparison_type,
                "comparison_id": comparison_id,
                "left_target_label": left_target,
                "left_biosample_term_name": left_biosample,
                "left_experiment_accession": left_exp,
                "left_selected_file_accession": left_file,
                "right_target_label": right_target,
                "right_biosample_term_name": right_biosample,
                "right_experiment_accession": right_exp,
                "right_selected_file_accession": right_file,
                "left_standardized_input_path": "",
                "right_standardized_input_path": "",
                "left_annotation_tsv_path": "",
                "right_annotation_tsv_path": "",
                "left_materialized_bed_path": "",
                "right_materialized_bed_path": "",
                "left_materialized_annotation_path": "",
                "right_materialized_annotation_path": "",
                "status": "missing_standardized_record",
                "note": "left_or_right_selected_file_not_found_in_standardized_summary",
            })
            continue

        if left_ann is None or right_ann is None:
            output_rows.append({
                "priority_rank": priority_rank,
                "priority_level": priority_level,
                "comparison_scope": comparison_scope,
                "comparison_type": comparison_type,
                "comparison_id": comparison_id,
                "left_target_label": left_target,
                "left_biosample_term_name": left_biosample,
                "left_experiment_accession": left_exp,
                "left_selected_file_accession": left_file,
                "right_target_label": right_target,
                "right_biosample_term_name": right_biosample,
                "right_experiment_accession": right_exp,
                "right_selected_file_accession": right_file,
                "left_standardized_input_path": normalize_text(left_std.get("output_standardized_path")),
                "right_standardized_input_path": normalize_text(right_std.get("output_standardized_path")),
                "left_annotation_tsv_path": "",
                "right_annotation_tsv_path": "",
                "left_materialized_bed_path": "",
                "right_materialized_bed_path": "",
                "left_materialized_annotation_path": "",
                "right_materialized_annotation_path": "",
                "status": "missing_annotation_record",
                "note": "left_or_right_selected_file_not_found_in_annotation_file_manifest",
            })
            continue

        left_std_rel = normalize_text(left_std.get("output_standardized_path"))
        right_std_rel = normalize_text(right_std.get("output_standardized_path"))
        left_ann_rel = normalize_text(left_ann.get("output_annotation_tsv_path"))
        right_ann_rel = normalize_text(right_ann.get("output_annotation_tsv_path"))

        left_std_abs = root / left_std_rel
        right_std_abs = root / right_std_rel
        left_ann_abs = root / left_ann_rel
        right_ann_abs = root / right_ann_rel

        if not left_std_abs.exists() or not right_std_abs.exists():
            output_rows.append({
                "priority_rank": priority_rank,
                "priority_level": priority_level,
                "comparison_scope": comparison_scope,
                "comparison_type": comparison_type,
                "comparison_id": comparison_id,
                "left_target_label": left_target,
                "left_biosample_term_name": left_biosample,
                "left_experiment_accession": left_exp,
                "left_selected_file_accession": left_file,
                "right_target_label": right_target,
                "right_biosample_term_name": right_biosample,
                "right_experiment_accession": right_exp,
                "right_selected_file_accession": right_file,
                "left_standardized_input_path": left_std_rel,
                "right_standardized_input_path": right_std_rel,
                "left_annotation_tsv_path": left_ann_rel,
                "right_annotation_tsv_path": right_ann_rel,
                "left_materialized_bed_path": "",
                "right_materialized_bed_path": "",
                "left_materialized_annotation_path": "",
                "right_materialized_annotation_path": "",
                "status": "missing_standardized_file",
                "note": "left_or_right_standardized_bed_not_found",
            })
            continue

        if not left_ann_abs.exists() or not right_ann_abs.exists():
            output_rows.append({
                "priority_rank": priority_rank,
                "priority_level": priority_level,
                "comparison_scope": comparison_scope,
                "comparison_type": comparison_type,
                "comparison_id": comparison_id,
                "left_target_label": left_target,
                "left_biosample_term_name": left_biosample,
                "left_experiment_accession": left_exp,
                "left_selected_file_accession": left_file,
                "right_target_label": right_target,
                "right_biosample_term_name": right_biosample,
                "right_experiment_accession": right_exp,
                "right_selected_file_accession": right_file,
                "left_standardized_input_path": left_std_rel,
                "right_standardized_input_path": right_std_rel,
                "left_annotation_tsv_path": left_ann_rel,
                "right_annotation_tsv_path": right_ann_rel,
                "left_materialized_bed_path": "",
                "right_materialized_bed_path": "",
                "left_materialized_annotation_path": "",
                "right_materialized_annotation_path": "",
                "status": "missing_annotation_file",
                "note": "left_or_right_annotation_tsv_not_found",
            })
            continue

        pair_dir = out_dir / "pairs" / safe_filename_part(comparison_id)

        left_bed_dst = pair_dir / "left.standardized.bed"
        right_bed_dst = pair_dir / "right.standardized.bed"
        left_ann_dst = pair_dir / "left.annotations.tsv"
        right_ann_dst = pair_dir / "right.annotations.tsv"

        left_bed_mode = materialize_file(left_std_abs, left_bed_dst, log_path)
        right_bed_mode = materialize_file(right_std_abs, right_bed_dst, log_path)
        left_ann_mode = materialize_file(left_ann_abs, left_ann_dst, log_path)
        right_ann_mode = materialize_file(right_ann_abs, right_ann_dst, log_path)

        output_rows.append({
            "priority_rank": priority_rank,
            "priority_level": priority_level,
            "comparison_scope": comparison_scope,
            "comparison_type": comparison_type,
            "comparison_id": comparison_id,
            "left_target_label": left_target,
            "left_biosample_term_name": left_biosample,
            "left_experiment_accession": left_exp,
            "left_selected_file_accession": left_file,
            "right_target_label": right_target,
            "right_biosample_term_name": right_biosample,
            "right_experiment_accession": right_exp,
            "right_selected_file_accession": right_file,
            "left_standardized_input_path": left_std_rel,
            "right_standardized_input_path": right_std_rel,
            "left_annotation_tsv_path": left_ann_rel,
            "right_annotation_tsv_path": right_ann_rel,
            "left_materialized_bed_path": str(left_bed_dst.relative_to(root)),
            "right_materialized_bed_path": str(right_bed_dst.relative_to(root)),
            "left_materialized_annotation_path": str(left_ann_dst.relative_to(root)),
            "right_materialized_annotation_path": str(right_ann_dst.relative_to(root)),
            "status": "ok",
            "note": (
                f"materialized:left_bed={left_bed_mode},right_bed={right_bed_mode},"
                f"left_ann={left_ann_mode},right_ann={right_ann_mode}"
            ),
        })

    output_rows = sorted(
        output_rows,
        key=lambda r: (
            normalize_text(r.get("priority_level")),
            normalize_text(r.get("comparison_scope")),
            normalize_text(r.get("priority_rank")),
            normalize_text(r.get("comparison_id")),
        ),
    )

    write_tsv(manifest_out, output_rows, OUTPUT_COLUMNS)

    n_ok = sum(1 for r in output_rows if normalize_text(r.get("status")) == "ok")
    n_skip = sum(1 for r in output_rows if normalize_text(r.get("status")) == "skipped_non_sample_scope")
    n_bad = len(output_rows) - n_ok - n_skip

    log("", log_path)
    log("[DONE] Priority peak set extraction completed.", log_path)
    log(f"[DONE] Output manifest: {manifest_out}", log_path)
    log(f"[INFO] OK rows       : {n_ok}", log_path)
    log(f"[INFO] Skipped rows  : {n_skip}", log_path)
    log(f"[INFO] Problem rows  : {n_bad}", log_path)

    for row in output_rows:
        if normalize_text(row.get("status")) != "ok":
            continue
        log(
            "  " + " | ".join([
                normalize_text(row.get("priority_level")),
                normalize_text(row.get("comparison_id")),
                normalize_text(row.get("left_materialized_bed_path")),
                normalize_text(row.get("right_materialized_bed_path")),
            ]),
            log_path,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())