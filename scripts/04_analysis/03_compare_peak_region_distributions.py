#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare peak region composition across samples, targets, and biosamples.

Inputs
------
- processed/peaks/annotations/region_distribution/sample_region_distribution_wide.tsv
- processed/peaks/annotations/region_distribution/target_biosample_region_distribution.tsv

Outputs
-------
- processed/peaks/annotations/region_distribution/comparisons/sample_pairwise_region_differences_long.tsv
- processed/peaks/annotations/region_distribution/comparisons/sample_pairwise_region_differences_summary.tsv
- processed/peaks/annotations/region_distribution/comparisons/aggregated_pairwise_region_differences_long.tsv
- processed/peaks/annotations/region_distribution/comparisons/aggregated_pairwise_region_differences_summary.tsv
- logs/compare_peak_region_distributions.log

Comparison logic
----------------
1. Sample-level comparisons:
   - same target, different biosample
   - same biosample, different target

2. Aggregated comparisons (target+biosample level):
   - same target, different biosample
   - same biosample, different target

All comparisons are based on region fractions / percentages.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple


REGIONS = [
    "five_prime_utr",
    "three_prime_utr",
    "utr",
    "CDS",
    "exon",
    "intron",
    "gene",
    "intergenic",
]


SAMPLE_LONG_COLUMNS = [
    "comparison_type",
    "left_target_label",
    "left_biosample_term_name",
    "left_experiment_accession",
    "left_selected_file_accession",
    "right_target_label",
    "right_biosample_term_name",
    "right_experiment_accession",
    "right_selected_file_accession",
    "region_class",
    "left_peak_count",
    "right_peak_count",
    "left_fraction",
    "right_fraction",
    "delta_fraction",
    "abs_delta_fraction",
    "left_percent",
    "right_percent",
    "delta_percent",
    "abs_delta_percent",
]

SAMPLE_SUMMARY_COLUMNS = [
    "comparison_type",
    "left_target_label",
    "left_biosample_term_name",
    "left_experiment_accession",
    "left_selected_file_accession",
    "right_target_label",
    "right_biosample_term_name",
    "right_experiment_accession",
    "right_selected_file_accession",
    "left_n_peaks",
    "right_n_peaks",
    "dominant_shift_region",
    "dominant_shift_delta_fraction",
    "dominant_shift_delta_percent",
    "dominant_shift_abs_delta_percent",
    "left_dominant_region",
    "left_dominant_region_percent",
    "right_dominant_region",
    "right_dominant_region_percent",
]

AGG_LONG_COLUMNS = [
    "comparison_type",
    "left_target_label",
    "left_biosample_term_name",
    "right_target_label",
    "right_biosample_term_name",
    "region_class",
    "left_peak_count",
    "right_peak_count",
    "left_fraction",
    "right_fraction",
    "delta_fraction",
    "abs_delta_fraction",
    "left_percent",
    "right_percent",
    "delta_percent",
    "abs_delta_percent",
]

AGG_SUMMARY_COLUMNS = [
    "comparison_type",
    "left_target_label",
    "left_biosample_term_name",
    "right_target_label",
    "right_biosample_term_name",
    "left_n_peaks",
    "right_n_peaks",
    "dominant_shift_region",
    "dominant_shift_delta_fraction",
    "dominant_shift_delta_percent",
    "dominant_shift_abs_delta_percent",
    "left_dominant_region",
    "left_dominant_region_percent",
    "right_dominant_region",
    "right_dominant_region_percent",
]


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, log_path: Path | None = None) -> None:
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


def safe_int(x: Any, default: int = 0) -> int:
    try:
        return int(str(x).strip())
    except Exception:
        return default


def safe_float(x: Any, default: float = 0.0) -> float:
    try:
        return float(str(x).strip())
    except Exception:
        return default


def region_counts(row: Dict[str, Any]) -> Dict[str, int]:
    return {region: safe_int(row.get(f"n_{region}", 0), 0) for region in REGIONS}


def region_fracs(row: Dict[str, Any]) -> Dict[str, float]:
    return {region: safe_float(row.get(f"frac_{region}", 0.0), 0.0) for region in REGIONS}


def dominant_region(row: Dict[str, Any]) -> Tuple[str, float]:
    best_region = "NA"
    best_frac = -1.0
    for region in REGIONS:
        frac = safe_float(row.get(f"frac_{region}", 0.0), 0.0)
        if frac > best_frac:
            best_region = region
            best_frac = frac
    return best_region, best_frac


def sample_pair_comparison_type(left: Dict[str, Any], right: Dict[str, Any]) -> str | None:
    lt = normalize_text(left.get("target_label"))
    rt = normalize_text(right.get("target_label"))
    lb = normalize_text(left.get("biosample_term_name"))
    rb = normalize_text(right.get("biosample_term_name"))

    if lt == rt and lb != rb:
        return "sample_same_target_diff_biosample"
    if lb == rb and lt != rt:
        return "sample_same_biosample_diff_target"
    return None


def aggregated_pair_comparison_type(left: Dict[str, Any], right: Dict[str, Any]) -> str | None:
    lt = normalize_text(left.get("target_label"))
    rt = normalize_text(right.get("target_label"))
    lb = normalize_text(left.get("biosample_term_name"))
    rb = normalize_text(right.get("biosample_term_name"))

    if lt == rt and lb != rb:
        return "aggregated_same_target_diff_biosample"
    if lb == rb and lt != rt:
        return "aggregated_same_biosample_diff_target"
    return None


def compare_two_rows_sample(left: Dict[str, Any], right: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    comparison_type = sample_pair_comparison_type(left, right)
    if comparison_type is None:
        return [], {}

    left_counts = region_counts(left)
    right_counts = region_counts(right)
    left_fr = region_fracs(left)
    right_fr = region_fracs(right)

    long_rows: List[Dict[str, Any]] = []
    max_abs_delta = -1.0
    max_region = "NA"
    max_delta = 0.0

    for region in REGIONS:
        delta_frac = right_fr[region] - left_fr[region]
        abs_delta_frac = abs(delta_frac)

        if abs_delta_frac > max_abs_delta:
            max_abs_delta = abs_delta_frac
            max_region = region
            max_delta = delta_frac

        long_rows.append({
            "comparison_type": comparison_type,
            "left_target_label": normalize_text(left.get("target_label")),
            "left_biosample_term_name": normalize_text(left.get("biosample_term_name")),
            "left_experiment_accession": normalize_text(left.get("experiment_accession")),
            "left_selected_file_accession": normalize_text(left.get("selected_file_accession")),
            "right_target_label": normalize_text(right.get("target_label")),
            "right_biosample_term_name": normalize_text(right.get("biosample_term_name")),
            "right_experiment_accession": normalize_text(right.get("experiment_accession")),
            "right_selected_file_accession": normalize_text(right.get("selected_file_accession")),
            "region_class": region,
            "left_peak_count": left_counts[region],
            "right_peak_count": right_counts[region],
            "left_fraction": round(left_fr[region], 6),
            "right_fraction": round(right_fr[region], 6),
            "delta_fraction": round(delta_frac, 6),
            "abs_delta_fraction": round(abs_delta_frac, 6),
            "left_percent": round(left_fr[region] * 100, 4),
            "right_percent": round(right_fr[region] * 100, 4),
            "delta_percent": round(delta_frac * 100, 4),
            "abs_delta_percent": round(abs_delta_frac * 100, 4),
        })

    left_dom, left_dom_frac = dominant_region(left)
    right_dom, right_dom_frac = dominant_region(right)

    summary = {
        "comparison_type": comparison_type,
        "left_target_label": normalize_text(left.get("target_label")),
        "left_biosample_term_name": normalize_text(left.get("biosample_term_name")),
        "left_experiment_accession": normalize_text(left.get("experiment_accession")),
        "left_selected_file_accession": normalize_text(left.get("selected_file_accession")),
        "right_target_label": normalize_text(right.get("target_label")),
        "right_biosample_term_name": normalize_text(right.get("biosample_term_name")),
        "right_experiment_accession": normalize_text(right.get("experiment_accession")),
        "right_selected_file_accession": normalize_text(right.get("selected_file_accession")),
        "left_n_peaks": safe_int(left.get("n_peaks", 0), 0),
        "right_n_peaks": safe_int(right.get("n_peaks", 0), 0),
        "dominant_shift_region": max_region,
        "dominant_shift_delta_fraction": round(max_delta, 6),
        "dominant_shift_delta_percent": round(max_delta * 100, 4),
        "dominant_shift_abs_delta_percent": round(abs(max_delta) * 100, 4),
        "left_dominant_region": left_dom,
        "left_dominant_region_percent": round(left_dom_frac * 100, 4),
        "right_dominant_region": right_dom,
        "right_dominant_region_percent": round(right_dom_frac * 100, 4),
    }
    return long_rows, summary


def compare_two_rows_aggregated(left: Dict[str, Any], right: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    comparison_type = aggregated_pair_comparison_type(left, right)
    if comparison_type is None:
        return [], {}

    left_counts = region_counts(left)
    right_counts = region_counts(right)
    left_fr = region_fracs(left)
    right_fr = region_fracs(right)

    long_rows: List[Dict[str, Any]] = []
    max_abs_delta = -1.0
    max_region = "NA"
    max_delta = 0.0

    for region in REGIONS:
        delta_frac = right_fr[region] - left_fr[region]
        abs_delta_frac = abs(delta_frac)

        if abs_delta_frac > max_abs_delta:
            max_abs_delta = abs_delta_frac
            max_region = region
            max_delta = delta_frac

        long_rows.append({
            "comparison_type": comparison_type,
            "left_target_label": normalize_text(left.get("target_label")),
            "left_biosample_term_name": normalize_text(left.get("biosample_term_name")),
            "right_target_label": normalize_text(right.get("target_label")),
            "right_biosample_term_name": normalize_text(right.get("biosample_term_name")),
            "region_class": region,
            "left_peak_count": left_counts[region],
            "right_peak_count": right_counts[region],
            "left_fraction": round(left_fr[region], 6),
            "right_fraction": round(right_fr[region], 6),
            "delta_fraction": round(delta_frac, 6),
            "abs_delta_fraction": round(abs_delta_frac, 6),
            "left_percent": round(left_fr[region] * 100, 4),
            "right_percent": round(right_fr[region] * 100, 4),
            "delta_percent": round(delta_frac * 100, 4),
            "abs_delta_percent": round(abs_delta_frac * 100, 4),
        })

    left_dom, left_dom_frac = dominant_region(left)
    right_dom, right_dom_frac = dominant_region(right)

    summary = {
        "comparison_type": comparison_type,
        "left_target_label": normalize_text(left.get("target_label")),
        "left_biosample_term_name": normalize_text(left.get("biosample_term_name")),
        "right_target_label": normalize_text(right.get("target_label")),
        "right_biosample_term_name": normalize_text(right.get("biosample_term_name")),
        "left_n_peaks": safe_int(left.get("n_peaks", 0), 0),
        "right_n_peaks": safe_int(right.get("n_peaks", 0), 0),
        "dominant_shift_region": max_region,
        "dominant_shift_delta_fraction": round(max_delta, 6),
        "dominant_shift_delta_percent": round(max_delta * 100, 4),
        "dominant_shift_abs_delta_percent": round(abs(max_delta) * 100, 4),
        "left_dominant_region": left_dom,
        "left_dominant_region_percent": round(left_dom_frac * 100, 4),
        "right_dominant_region": right_dom,
        "right_dominant_region_percent": round(right_dom_frac * 100, 4),
    }
    return long_rows, summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare peak region distributions across samples, targets, and biosamples."
    )
    parser.add_argument(
        "--sample-wide",
        default=None,
        help="Override sample_region_distribution_wide.tsv",
    )
    parser.add_argument(
        "--target-biosample",
        default=None,
        help="Override target_biosample_region_distribution.tsv",
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

    region_dir = root / "processed" / "peaks" / "annotations" / "region_distribution"
    sample_wide_path = (
        Path(args.sample_wide).resolve()
        if args.sample_wide
        else region_dir / "sample_region_distribution_wide.tsv"
    )
    target_biosample_path = (
        Path(args.target_biosample).resolve()
        if args.target_biosample
        else region_dir / "target_biosample_region_distribution.tsv"
    )

    out_dir = region_dir / "comparisons"
    log_path = root / "logs" / "compare_peak_region_distributions.log"

    sample_long_out = out_dir / "sample_pairwise_region_differences_long.tsv"
    sample_summary_out = out_dir / "sample_pairwise_region_differences_summary.tsv"
    agg_long_out = out_dir / "aggregated_pairwise_region_differences_long.tsv"
    agg_summary_out = out_dir / "aggregated_pairwise_region_differences_summary.tsv"

    sample_rows = read_tsv_rows(sample_wide_path)
    agg_rows = read_tsv_rows(target_biosample_path)

    if args.targets:
        target_filter = {x.strip().lower() for x in args.targets if x.strip()}
        sample_rows = [r for r in sample_rows if normalize_text(r.get("target_label")).lower() in target_filter]
        agg_rows = [r for r in agg_rows if normalize_text(r.get("target_label")).lower() in target_filter]

    if args.biosamples:
        biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()}
        sample_rows = [r for r in sample_rows if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter]
        agg_rows = [r for r in agg_rows if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter]

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Sample wide input: {sample_wide_path}", log_path)
    log(f"[INFO] Aggregated input: {target_biosample_path}", log_path)
    log(f"[INFO] Sample rows: {len(sample_rows)}", log_path)
    log(f"[INFO] Aggregated rows: {len(agg_rows)}", log_path)

    sample_long_rows: List[Dict[str, Any]] = []
    sample_summary_rows: List[Dict[str, Any]] = []

    for left, right in itertools.combinations(sample_rows, 2):
        long_rows, summary = compare_two_rows_sample(left, right)
        if long_rows:
            sample_long_rows.extend(long_rows)
            sample_summary_rows.append(summary)

    agg_long_rows: List[Dict[str, Any]] = []
    agg_summary_rows: List[Dict[str, Any]] = []

    for left, right in itertools.combinations(agg_rows, 2):
        long_rows, summary = compare_two_rows_aggregated(left, right)
        if long_rows:
            agg_long_rows.extend(long_rows)
            agg_summary_rows.append(summary)

    sample_long_rows = sorted(
        sample_long_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_type")),
            normalize_text(r.get("left_target_label")),
            normalize_text(r.get("left_biosample_term_name")),
            normalize_text(r.get("left_experiment_accession")),
            normalize_text(r.get("right_target_label")),
            normalize_text(r.get("right_biosample_term_name")),
            normalize_text(r.get("right_experiment_accession")),
            normalize_text(r.get("region_class")),
        ),
    )

    sample_summary_rows = sorted(
        sample_summary_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_type")),
            -safe_float(r.get("dominant_shift_abs_delta_percent"), 0.0),
            normalize_text(r.get("left_target_label")),
            normalize_text(r.get("left_biosample_term_name")),
            normalize_text(r.get("right_target_label")),
            normalize_text(r.get("right_biosample_term_name")),
        ),
    )

    agg_long_rows = sorted(
        agg_long_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_type")),
            normalize_text(r.get("left_target_label")),
            normalize_text(r.get("left_biosample_term_name")),
            normalize_text(r.get("right_target_label")),
            normalize_text(r.get("right_biosample_term_name")),
            normalize_text(r.get("region_class")),
        ),
    )

    agg_summary_rows = sorted(
        agg_summary_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_type")),
            -safe_float(r.get("dominant_shift_abs_delta_percent"), 0.0),
            normalize_text(r.get("left_target_label")),
            normalize_text(r.get("left_biosample_term_name")),
            normalize_text(r.get("right_target_label")),
            normalize_text(r.get("right_biosample_term_name")),
        ),
    )

    write_tsv(sample_long_out, sample_long_rows, SAMPLE_LONG_COLUMNS)
    write_tsv(sample_summary_out, sample_summary_rows, SAMPLE_SUMMARY_COLUMNS)
    write_tsv(agg_long_out, agg_long_rows, AGG_LONG_COLUMNS)
    write_tsv(agg_summary_out, agg_summary_rows, AGG_SUMMARY_COLUMNS)

    log("", log_path)
    log("[DONE] Peak region distribution comparison completed.", log_path)
    log(f"[DONE] Sample long TSV: {sample_long_out}", log_path)
    log(f"[DONE] Sample summary TSV: {sample_summary_out}", log_path)
    log(f"[DONE] Aggregated long TSV: {agg_long_out}", log_path)
    log(f"[DONE] Aggregated summary TSV: {agg_summary_out}", log_path)
    log("", log_path)

    log("[INFO] Top sample-level shifts:", log_path)
    for row in sample_summary_rows[:10]:
        log(
            "  " + " | ".join([
                normalize_text(row.get("comparison_type")),
                f"{normalize_text(row.get('left_target_label'))}/{normalize_text(row.get('left_biosample_term_name'))}/{normalize_text(row.get('left_experiment_accession'))}",
                f"{normalize_text(row.get('right_target_label'))}/{normalize_text(row.get('right_biosample_term_name'))}/{normalize_text(row.get('right_experiment_accession'))}",
                f"region={normalize_text(row.get('dominant_shift_region'))}",
                f"abs_delta_pct={normalize_text(row.get('dominant_shift_abs_delta_percent'))}",
            ]),
            log_path,
        )

    log("[INFO] Top aggregated shifts:", log_path)
    for row in agg_summary_rows[:10]:
        log(
            "  " + " | ".join([
                normalize_text(row.get("comparison_type")),
                f"{normalize_text(row.get('left_target_label'))}/{normalize_text(row.get('left_biosample_term_name'))}",
                f"{normalize_text(row.get('right_target_label'))}/{normalize_text(row.get('right_biosample_term_name'))}",
                f"region={normalize_text(row.get('dominant_shift_region'))}",
                f"abs_delta_pct={normalize_text(row.get('dominant_shift_abs_delta_percent'))}",
            ]),
            log_path,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())