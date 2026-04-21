#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build overlap-based peak sets for priority comparison pairs.

Inputs
------
- processed/peaks/comparison_sets/comparison_peak_set_manifest.tsv
- processed/peaks/annotation_inputs/annotation_input_manifest.tsv

Required per-row files
----------------------
From comparison manifest:
- left_materialized_annotation_path
- right_materialized_annotation_path

From annotation input manifest (resolved by selected_file_accession):
- annotation_ready_bed6_path

Outputs
-------
- processed/peaks/comparison_sets/overlaps/overlap_summary.tsv
- processed/peaks/comparison_sets/overlaps/overlap_file_manifest.tsv
- processed/peaks/comparison_sets/overlaps/<comparison_id>/shared_peak_pairs.tsv
- processed/peaks/comparison_sets/overlaps/<comparison_id>/shared_left.annotation_ready.bed
- processed/peaks/comparison_sets/overlaps/<comparison_id>/shared_right.annotation_ready.bed
- processed/peaks/comparison_sets/overlaps/<comparison_id>/left_specific.annotation_ready.bed
- processed/peaks/comparison_sets/overlaps/<comparison_id>/right_specific.annotation_ready.bed
- processed/peaks/comparison_sets/overlaps/<comparison_id>/shared_left.annotations.tsv
- processed/peaks/comparison_sets/overlaps/<comparison_id>/shared_right.annotations.tsv
- processed/peaks/comparison_sets/overlaps/<comparison_id>/left_specific.annotations.tsv
- processed/peaks/comparison_sets/overlaps/<comparison_id>/right_specific.annotations.tsv
- logs/build_peak_overlap_sets.log

Overlap rule
------------
Default: at least 1 bp overlap.
Optional reciprocal overlap threshold can be set.

Notes
-----
- Shared peaks are defined by interval overlap, not exact coordinate equality.
- This script uses annotation-ready BED6 files because they contain stable unique peak IDs
  that are consistent with per-file annotation TSVs.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


BIN_SIZE = 100_000


SUMMARY_COLUMNS = [
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
    "left_total_peaks",
    "right_total_peaks",
    "shared_left_peaks",
    "shared_right_peaks",
    "left_specific_peaks",
    "right_specific_peaks",
    "shared_pair_links",
    "left_shared_fraction",
    "right_shared_fraction",
    "left_specific_fraction",
    "right_specific_fraction",
    "status",
    "note",
]

FILE_MANIFEST_COLUMNS = [
    "comparison_id",
    "left_input_annotation_ready_bed",
    "right_input_annotation_ready_bed",
    "shared_pairs_tsv",
    "shared_left_bed",
    "shared_right_bed",
    "left_specific_bed",
    "right_specific_bed",
    "shared_left_annotations_tsv",
    "shared_right_annotations_tsv",
    "left_specific_annotations_tsv",
    "right_specific_annotations_tsv",
    "status",
    "note",
]

PAIR_COLUMNS = [
    "left_peak_id",
    "left_chrom",
    "left_start",
    "left_end",
    "right_peak_id",
    "right_chrom",
    "right_start",
    "right_end",
    "overlap_bp",
    "left_peak_length",
    "right_peak_length",
    "left_reciprocal_overlap",
    "right_reciprocal_overlap",
]

ANNOTATION_COLUMNS = [
    "peak_id",
    "chrom",
    "start",
    "end",
    "score",
    "strand",
    "target_label",
    "biosample_term_name",
    "experiment_accession",
    "selected_file_accession",
    "primary_region_class",
    "overlap_gene_ids",
    "overlap_gene_names",
    "overlap_transcript_ids",
    "n_primary_region_features",
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


def build_index(rows: List[Dict[str, str]], key: str) -> Dict[str, Dict[str, str]]:
    idx: Dict[str, Dict[str, str]] = {}
    for row in rows:
        k = normalize_text(row.get(key))
        if k:
            idx[k] = row
    return idx


def read_bed6_rows(path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as f:
        for i, raw in enumerate(f, start=1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                fields = line.split()
            if len(fields) < 6:
                continue

            chrom = fields[0]
            start = safe_int(fields[1], None)
            end = safe_int(fields[2], None)
            peak_id = fields[3]
            score = fields[4]
            strand = fields[5]

            if start is None or end is None or end <= start:
                continue

            rows.append({
                "row_index": i,
                "chrom": chrom,
                "start": start,
                "end": end,
                "peak_id": peak_id,
                "score": score,
                "strand": strand,
                "fields": fields,
                "line": line,
            })
    return rows


def read_annotation_index(path: Path) -> Dict[str, Dict[str, str]]:
    rows = read_tsv_rows(path)
    idx: Dict[str, Dict[str, str]] = {}
    for row in rows:
        pid = normalize_text(row.get("peak_id"))
        if pid:
            idx[pid] = row
    return idx


def build_interval_index(rows: List[Dict[str, Any]]) -> Dict[str, Dict[int, List[int]]]:
    idx: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(list))
    for i, row in enumerate(rows):
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        start_bin = start // BIN_SIZE
        end_bin = (end - 1) // BIN_SIZE if end > 0 else start_bin
        for b in range(start_bin, end_bin + 1):
            idx[chrom][b].append(i)
    return idx


def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def find_overlaps(
    left_rows: List[Dict[str, Any]],
    right_rows: List[Dict[str, Any]],
    right_index: Dict[str, Dict[int, List[int]]],
    min_overlap_bp: int,
    min_reciprocal_overlap: float,
) -> Tuple[List[Dict[str, Any]], Set[str], Set[str]]:
    pair_rows: List[Dict[str, Any]] = []
    left_shared_ids: Set[str] = set()
    right_shared_ids: Set[str] = set()

    for left in left_rows:
        chrom = left["chrom"]
        if chrom not in right_index:
            continue

        start = left["start"]
        end = left["end"]
        left_len = end - start
        if left_len <= 0:
            continue

        start_bin = start // BIN_SIZE
        end_bin = (end - 1) // BIN_SIZE if end > 0 else start_bin
        candidate_indices: Set[int] = set()
        for b in range(start_bin, end_bin + 1):
            candidate_indices.update(right_index[chrom].get(b, []))

        for j in candidate_indices:
            right = right_rows[j]
            ov = overlap_bp(start, end, right["start"], right["end"])
            if ov < min_overlap_bp:
                continue

            right_len = right["end"] - right["start"]
            if right_len <= 0:
                continue

            left_recip = ov / left_len
            right_recip = ov / right_len

            if left_recip < min_reciprocal_overlap or right_recip < min_reciprocal_overlap:
                continue

            left_shared_ids.add(left["peak_id"])
            right_shared_ids.add(right["peak_id"])

            pair_rows.append({
                "left_peak_id": left["peak_id"],
                "left_chrom": left["chrom"],
                "left_start": left["start"],
                "left_end": left["end"],
                "right_peak_id": right["peak_id"],
                "right_chrom": right["chrom"],
                "right_start": right["start"],
                "right_end": right["end"],
                "overlap_bp": ov,
                "left_peak_length": left_len,
                "right_peak_length": right_len,
                "left_reciprocal_overlap": round(left_recip, 6),
                "right_reciprocal_overlap": round(right_recip, 6),
            })

    pair_rows = sorted(
        pair_rows,
        key=lambda r: (
            normalize_text(r.get("left_chrom")),
            safe_int(r.get("left_start"), 0),
            safe_int(r.get("left_end"), 0),
            normalize_text(r.get("left_peak_id")),
            normalize_text(r.get("right_peak_id")),
        ),
    )
    return pair_rows, left_shared_ids, right_shared_ids


def subset_bed_rows(rows: List[Dict[str, Any]], keep_ids: Set[str]) -> List[Dict[str, Any]]:
    return [r for r in rows if r["peak_id"] in keep_ids]


def subset_annotation_rows(
    annotation_idx: Dict[str, Dict[str, str]],
    keep_ids: Set[str],
) -> List[Dict[str, str]]:
    rows = [annotation_idx[pid] for pid in keep_ids if pid in annotation_idx]
    rows = sorted(
        rows,
        key=lambda r: (
            normalize_text(r.get("chrom")),
            safe_int(r.get("start"), 0),
            safe_int(r.get("end"), 0),
            normalize_text(r.get("peak_id")),
        ),
    )
    return rows


def write_bed(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for row in rows:
            f.write(row["line"] + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build shared/specific peak overlap sets for priority comparison pairs."
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Override processed/peaks/comparison_sets/comparison_peak_set_manifest.tsv",
    )
    parser.add_argument(
        "--annotation-input-manifest",
        default=None,
        help="Override processed/peaks/annotation_inputs/annotation_input_manifest.tsv",
    )
    parser.add_argument(
        "--priority-levels",
        nargs="*",
        default=None,
        help="Restrict to these priority levels (e.g. priority_1)",
    )
    parser.add_argument(
        "--min-overlap-bp",
        type=int,
        default=1,
        help="Minimum overlap length in bp (default: 1)",
    )
    parser.add_argument(
        "--min-reciprocal-overlap",
        type=float,
        default=0.0,
        help="Minimum reciprocal overlap on both sides (default: 0.0)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process first N eligible comparisons",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    manifest_path = (
        Path(args.manifest).resolve()
        if args.manifest
        else root / "processed" / "peaks" / "comparison_sets" / "comparison_peak_set_manifest.tsv"
    )
    annotation_input_manifest = (
        Path(args.annotation_input_manifest).resolve()
        if args.annotation_input_manifest
        else root / "processed" / "peaks" / "annotation_inputs" / "annotation_input_manifest.tsv"
    )

    out_dir = root / "processed" / "peaks" / "comparison_sets" / "overlaps"
    summary_out = out_dir / "overlap_summary.tsv"
    file_manifest_out = out_dir / "overlap_file_manifest.tsv"
    log_path = root / "logs" / "build_peak_overlap_sets.log"

    rows = read_tsv_rows(manifest_path)
    annotation_input_rows = read_tsv_rows(annotation_input_manifest)
    annotation_input_idx = build_index(annotation_input_rows, "selected_file_accession")

    if args.priority_levels:
        allowed = {normalize_text(x) for x in args.priority_levels if normalize_text(x)}
        rows = [r for r in rows if normalize_text(r.get("priority_level")) in allowed]

    rows = [r for r in rows if normalize_text(r.get("status")) == "ok"]

    rows = sorted(
        rows,
        key=lambda r: (
            normalize_text(r.get("priority_level")),
            normalize_text(r.get("priority_rank")),
            normalize_text(r.get("comparison_id")),
        ),
    )

    if args.limit is not None:
        rows = rows[:args.limit]

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Comparison manifest: {manifest_path}", log_path)
    log(f"[INFO] Annotation input manifest: {annotation_input_manifest}", log_path)
    log(f"[INFO] Eligible comparisons: {len(rows)}", log_path)
    log(f"[INFO] min_overlap_bp: {args.min_overlap_bp}", log_path)
    log(f"[INFO] min_reciprocal_overlap: {args.min_reciprocal_overlap}", log_path)

    summary_rows: List[Dict[str, Any]] = []
    file_manifest_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(rows, start=1):
        comparison_id = normalize_text(row.get("comparison_id"))
        comp_dir = out_dir / safe_filename_part(comparison_id)

        left_file = normalize_text(row.get("left_selected_file_accession"))
        right_file = normalize_text(row.get("right_selected_file_accession"))

        left_ann = root / normalize_text(row.get("left_materialized_annotation_path"))
        right_ann = root / normalize_text(row.get("right_materialized_annotation_path"))

        left_input_row = annotation_input_idx.get(left_file)
        right_input_row = annotation_input_idx.get(right_file)

        log(f"[INFO] [{i}/{len(rows)}] Building overlap sets for {comparison_id}", log_path)

        if left_input_row is None or right_input_row is None:
            summary_rows.append({
                "priority_rank": normalize_text(row.get("priority_rank")),
                "priority_level": normalize_text(row.get("priority_level")),
                "comparison_scope": normalize_text(row.get("comparison_scope")),
                "comparison_type": normalize_text(row.get("comparison_type")),
                "comparison_id": comparison_id,
                "left_target_label": normalize_text(row.get("left_target_label")),
                "left_biosample_term_name": normalize_text(row.get("left_biosample_term_name")),
                "left_experiment_accession": normalize_text(row.get("left_experiment_accession")),
                "left_selected_file_accession": left_file,
                "right_target_label": normalize_text(row.get("right_target_label")),
                "right_biosample_term_name": normalize_text(row.get("right_biosample_term_name")),
                "right_experiment_accession": normalize_text(row.get("right_experiment_accession")),
                "right_selected_file_accession": right_file,
                "left_total_peaks": 0,
                "right_total_peaks": 0,
                "shared_left_peaks": 0,
                "shared_right_peaks": 0,
                "left_specific_peaks": 0,
                "right_specific_peaks": 0,
                "shared_pair_links": 0,
                "left_shared_fraction": 0.0,
                "right_shared_fraction": 0.0,
                "left_specific_fraction": 0.0,
                "right_specific_fraction": 0.0,
                "status": "missing_annotation_ready_bed_record",
                "note": "left_or_right_selected_file_not_found_in_annotation_input_manifest",
            })
            continue

        left_bed6 = root / normalize_text(left_input_row.get("annotation_ready_bed6_path"))
        right_bed6 = root / normalize_text(right_input_row.get("annotation_ready_bed6_path"))

        if not left_bed6.exists() or not right_bed6.exists() or not left_ann.exists() or not right_ann.exists():
            summary_rows.append({
                "priority_rank": normalize_text(row.get("priority_rank")),
                "priority_level": normalize_text(row.get("priority_level")),
                "comparison_scope": normalize_text(row.get("comparison_scope")),
                "comparison_type": normalize_text(row.get("comparison_type")),
                "comparison_id": comparison_id,
                "left_target_label": normalize_text(row.get("left_target_label")),
                "left_biosample_term_name": normalize_text(row.get("left_biosample_term_name")),
                "left_experiment_accession": normalize_text(row.get("left_experiment_accession")),
                "left_selected_file_accession": left_file,
                "right_target_label": normalize_text(row.get("right_target_label")),
                "right_biosample_term_name": normalize_text(row.get("right_biosample_term_name")),
                "right_experiment_accession": normalize_text(row.get("right_experiment_accession")),
                "right_selected_file_accession": right_file,
                "left_total_peaks": 0,
                "right_total_peaks": 0,
                "shared_left_peaks": 0,
                "shared_right_peaks": 0,
                "left_specific_peaks": 0,
                "right_specific_peaks": 0,
                "shared_pair_links": 0,
                "left_shared_fraction": 0.0,
                "right_shared_fraction": 0.0,
                "left_specific_fraction": 0.0,
                "right_specific_fraction": 0.0,
                "status": "missing_input",
                "note": "one_or_more_required_input_files_missing",
            })
            continue

        left_rows = read_bed6_rows(left_bed6)
        right_rows = read_bed6_rows(right_bed6)
        left_ann_idx = read_annotation_index(left_ann)
        right_ann_idx = read_annotation_index(right_ann)

        right_index = build_interval_index(right_rows)

        pair_rows, left_shared_ids, right_shared_ids = find_overlaps(
            left_rows=left_rows,
            right_rows=right_rows,
            right_index=right_index,
            min_overlap_bp=args.min_overlap_bp,
            min_reciprocal_overlap=args.min_reciprocal_overlap,
        )

        left_all_ids = {r["peak_id"] for r in left_rows}
        right_all_ids = {r["peak_id"] for r in right_rows}

        left_specific_ids = left_all_ids - left_shared_ids
        right_specific_ids = right_all_ids - right_shared_ids

        shared_left_rows = subset_bed_rows(left_rows, left_shared_ids)
        shared_right_rows = subset_bed_rows(right_rows, right_shared_ids)
        left_specific_rows = subset_bed_rows(left_rows, left_specific_ids)
        right_specific_rows = subset_bed_rows(right_rows, right_specific_ids)

        shared_left_ann_rows = subset_annotation_rows(left_ann_idx, left_shared_ids)
        shared_right_ann_rows = subset_annotation_rows(right_ann_idx, right_shared_ids)
        left_specific_ann_rows = subset_annotation_rows(left_ann_idx, left_specific_ids)
        right_specific_ann_rows = subset_annotation_rows(right_ann_idx, right_specific_ids)

        shared_pairs_path = comp_dir / "shared_peak_pairs.tsv"
        shared_left_bed_path = comp_dir / "shared_left.annotation_ready.bed"
        shared_right_bed_path = comp_dir / "shared_right.annotation_ready.bed"
        left_specific_bed_path = comp_dir / "left_specific.annotation_ready.bed"
        right_specific_bed_path = comp_dir / "right_specific.annotation_ready.bed"

        shared_left_ann_path = comp_dir / "shared_left.annotations.tsv"
        shared_right_ann_path = comp_dir / "shared_right.annotations.tsv"
        left_specific_ann_path = comp_dir / "left_specific.annotations.tsv"
        right_specific_ann_path = comp_dir / "right_specific.annotations.tsv"

        write_tsv(shared_pairs_path, pair_rows, PAIR_COLUMNS)
        write_bed(shared_left_bed_path, shared_left_rows)
        write_bed(shared_right_bed_path, shared_right_rows)
        write_bed(left_specific_bed_path, left_specific_rows)
        write_bed(right_specific_bed_path, right_specific_rows)

        write_tsv(shared_left_ann_path, shared_left_ann_rows, ANNOTATION_COLUMNS)
        write_tsv(shared_right_ann_path, shared_right_ann_rows, ANNOTATION_COLUMNS)
        write_tsv(left_specific_ann_path, left_specific_ann_rows, ANNOTATION_COLUMNS)
        write_tsv(right_specific_ann_path, right_specific_ann_rows, ANNOTATION_COLUMNS)

        left_total = len(left_rows)
        right_total = len(right_rows)
        shared_left_n = len(left_shared_ids)
        shared_right_n = len(right_shared_ids)
        left_specific_n = len(left_specific_ids)
        right_specific_n = len(right_specific_ids)
        pair_links_n = len(pair_rows)

        summary_rows.append({
            "priority_rank": normalize_text(row.get("priority_rank")),
            "priority_level": normalize_text(row.get("priority_level")),
            "comparison_scope": normalize_text(row.get("comparison_scope")),
            "comparison_type": normalize_text(row.get("comparison_type")),
            "comparison_id": comparison_id,
            "left_target_label": normalize_text(row.get("left_target_label")),
            "left_biosample_term_name": normalize_text(row.get("left_biosample_term_name")),
            "left_experiment_accession": normalize_text(row.get("left_experiment_accession")),
            "left_selected_file_accession": left_file,
            "right_target_label": normalize_text(row.get("right_target_label")),
            "right_biosample_term_name": normalize_text(row.get("right_biosample_term_name")),
            "right_experiment_accession": normalize_text(row.get("right_experiment_accession")),
            "right_selected_file_accession": right_file,
            "left_total_peaks": left_total,
            "right_total_peaks": right_total,
            "shared_left_peaks": shared_left_n,
            "shared_right_peaks": shared_right_n,
            "left_specific_peaks": left_specific_n,
            "right_specific_peaks": right_specific_n,
            "shared_pair_links": pair_links_n,
            "left_shared_fraction": round(shared_left_n / left_total, 6) if left_total > 0 else 0.0,
            "right_shared_fraction": round(shared_right_n / right_total, 6) if right_total > 0 else 0.0,
            "left_specific_fraction": round(left_specific_n / left_total, 6) if left_total > 0 else 0.0,
            "right_specific_fraction": round(right_specific_n / right_total, 6) if right_total > 0 else 0.0,
            "status": "ok",
            "note": f"shared_pairs={pair_links_n}",
        })

        file_manifest_rows.append({
            "comparison_id": comparison_id,
            "left_input_annotation_ready_bed": str(left_bed6.relative_to(root)),
            "right_input_annotation_ready_bed": str(right_bed6.relative_to(root)),
            "shared_pairs_tsv": str(shared_pairs_path.relative_to(root)),
            "shared_left_bed": str(shared_left_bed_path.relative_to(root)),
            "shared_right_bed": str(shared_right_bed_path.relative_to(root)),
            "left_specific_bed": str(left_specific_bed_path.relative_to(root)),
            "right_specific_bed": str(right_specific_bed_path.relative_to(root)),
            "shared_left_annotations_tsv": str(shared_left_ann_path.relative_to(root)),
            "shared_right_annotations_tsv": str(shared_right_ann_path.relative_to(root)),
            "left_specific_annotations_tsv": str(left_specific_ann_path.relative_to(root)),
            "right_specific_annotations_tsv": str(right_specific_ann_path.relative_to(root)),
            "status": "ok",
            "note": f"shared_pairs={pair_links_n}",
        })

    summary_rows = sorted(
        summary_rows,
        key=lambda r: (
            normalize_text(r.get("priority_level")),
            normalize_text(r.get("priority_rank")),
            normalize_text(r.get("comparison_id")),
        ),
    )
    file_manifest_rows = sorted(
        file_manifest_rows,
        key=lambda r: normalize_text(r.get("comparison_id")),
    )

    write_tsv(summary_out, summary_rows, SUMMARY_COLUMNS)
    write_tsv(file_manifest_out, file_manifest_rows, FILE_MANIFEST_COLUMNS)

    n_ok = sum(1 for r in summary_rows if normalize_text(r.get("status")) == "ok")
    n_bad = len(summary_rows) - n_ok

    log("", log_path)
    log("[DONE] Peak overlap set construction completed.", log_path)
    log(f"[DONE] Summary TSV: {summary_out}", log_path)
    log(f"[DONE] File manifest TSV: {file_manifest_out}", log_path)
    log(f"[INFO] OK rows      : {n_ok}", log_path)
    log(f"[INFO] Problem rows : {n_bad}", log_path)

    for row in summary_rows:
        if normalize_text(row.get("status")) != "ok":
            continue
        log(
            "  " + " | ".join([
                normalize_text(row.get("comparison_id")),
                f"left_total={normalize_text(row.get('left_total_peaks'))}",
                f"right_total={normalize_text(row.get('right_total_peaks'))}",
                f"shared_left={normalize_text(row.get('shared_left_peaks'))}",
                f"shared_right={normalize_text(row.get('shared_right_peaks'))}",
                f"left_specific={normalize_text(row.get('left_specific_peaks'))}",
                f"right_specific={normalize_text(row.get('right_specific_peaks'))}",
                f"pair_links={normalize_text(row.get('shared_pair_links'))}",
            ]),
            log_path,
        )

    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())