#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Define exploratory priority comparisons from peak region distribution comparison summaries.

Inputs
------
- processed/peaks/annotations/region_distribution/comparisons/sample_pairwise_region_differences_summary.tsv
- processed/peaks/annotations/region_distribution/comparisons/aggregated_pairwise_region_differences_summary.tsv

Outputs
-------
- results/tables/exploratory_priority_comparisons.tsv
- logs/define_priority_comparisons.log

Goal
----
Turn current comparison results into a ranked priority list for exploratory follow-up.

Priority logic (heuristic)
--------------------------
Base signal:
- dominant_shift_abs_delta_percent

Bonuses:
- same_target_diff_biosample: biologically cleaner for "same RBP across samples"
- dominant shift region in intron/CDS/exon: easier to interpret in current project
- sample-level comparison: more specific
- aggregated-level comparison: broader support

Priority levels
---------------
- priority_1: top N strongest comparisons
- priority_2: next M comparisons
- priority_3: remaining comparisons
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List


OUTPUT_COLUMNS = [
    "priority_rank",
    "priority_level",
    "priority_score",
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
    "dominant_shift_region",
    "dominant_shift_delta_percent",
    "dominant_shift_abs_delta_percent",
    "left_dominant_region",
    "left_dominant_region_percent",
    "right_dominant_region",
    "right_dominant_region_percent",
    "recommended_followup_focus",
    "recommended_peak_subset",
    "followup_note",
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


def safe_float(x: Any, default: float = 0.0) -> float:
    try:
        return float(str(x).strip())
    except Exception:
        return default


def base_priority_score(abs_delta_percent: float) -> float:
    return abs_delta_percent


def scope_bonus(scope: str) -> float:
    if scope == "sample":
        return 1.5
    if scope == "aggregated":
        return 1.0
    return 0.0


def comparison_type_bonus(comparison_type: str) -> float:
    text = normalize_text(comparison_type)
    if "same_target_diff_biosample" in text:
        return 2.5
    if "same_biosample_diff_target" in text:
        return 1.5
    return 0.0


def region_bonus(region: str) -> float:
    region = normalize_text(region)
    if region == "intron":
        return 1.5
    if region == "CDS":
        return 1.5
    if region == "exon":
        return 1.0
    if region in {"five_prime_utr", "three_prime_utr", "utr"}:
        return 0.8
    if region == "gene":
        return 0.5
    return 0.0


def recommended_followup_focus(comparison_type: str, region: str) -> str:
    ctype = normalize_text(comparison_type)
    region = normalize_text(region)

    if "same_target_diff_biosample" in ctype:
        return f"same_RBP_cross_biosample_region_shift:{region}"
    if "same_biosample_diff_target" in ctype:
        return f"same_biosample_cross_RBP_region_shift:{region}"
    return f"region_shift:{region}"


def recommended_peak_subset(comparison_type: str, region: str) -> str:
    ctype = normalize_text(comparison_type)
    region = normalize_text(region)

    if "same_target_diff_biosample" in ctype:
        return f"{region}_specific_and_shared_peaks_between_biosamples"
    if "same_biosample_diff_target" in ctype:
        return f"{region}_specific_and_shared_peaks_between_RBPs"
    return f"{region}_focused_peak_subsets"


def followup_note(comparison_type: str, region: str) -> str:
    ctype = normalize_text(comparison_type)
    region = normalize_text(region)

    if "same_target_diff_biosample" in ctype:
        return (
            f"Prioritize extracting {region} peaks for the same RBP across biosamples; "
            f"compare shared vs biosample-specific sets and map to genes."
        )
    if "same_biosample_diff_target" in ctype:
        return (
            f"Prioritize extracting {region} peaks for different RBPs in the same biosample; "
            f"compare RBP-specific sets and map to genes."
        )
    return (
        f"Prioritize extracting {region} peaks for this comparison and summarize candidate genes."
    )


def build_sample_candidates(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in rows:
        abs_delta = safe_float(row.get("dominant_shift_abs_delta_percent", 0.0), 0.0)
        delta = safe_float(row.get("dominant_shift_delta_percent", 0.0), 0.0)
        region = normalize_text(row.get("dominant_shift_region"))
        ctype = normalize_text(row.get("comparison_type"))

        score = (
            base_priority_score(abs_delta)
            + scope_bonus("sample")
            + comparison_type_bonus(ctype)
            + region_bonus(region)
        )

        left_target = normalize_text(row.get("left_target_label"))
        left_biosample = normalize_text(row.get("left_biosample_term_name"))
        left_exp = normalize_text(row.get("left_experiment_accession"))
        left_file = normalize_text(row.get("left_selected_file_accession"))
        right_target = normalize_text(row.get("right_target_label"))
        right_biosample = normalize_text(row.get("right_biosample_term_name"))
        right_exp = normalize_text(row.get("right_experiment_accession"))
        right_file = normalize_text(row.get("right_selected_file_accession"))

        comparison_id = (
            f"sample__{left_target}__{left_biosample}__{left_exp}__{left_file}"
            f"__VS__{right_target}__{right_biosample}__{right_exp}__{right_file}"
        )

        out.append({
            "comparison_scope": "sample",
            "comparison_type": ctype,
            "comparison_id": comparison_id,
            "left_target_label": left_target,
            "left_biosample_term_name": left_biosample,
            "left_experiment_accession": left_exp,
            "left_selected_file_accession": left_file,
            "right_target_label": right_target,
            "right_biosample_term_name": right_biosample,
            "right_experiment_accession": right_exp,
            "right_selected_file_accession": right_file,
            "dominant_shift_region": region,
            "dominant_shift_delta_percent": round(delta, 4),
            "dominant_shift_abs_delta_percent": round(abs_delta, 4),
            "left_dominant_region": normalize_text(row.get("left_dominant_region")),
            "left_dominant_region_percent": round(
                safe_float(row.get("left_dominant_region_percent", 0.0), 0.0), 4
            ),
            "right_dominant_region": normalize_text(row.get("right_dominant_region")),
            "right_dominant_region_percent": round(
                safe_float(row.get("right_dominant_region_percent", 0.0), 0.0), 4
            ),
            "priority_score": round(score, 4),
            "recommended_followup_focus": recommended_followup_focus(ctype, region),
            "recommended_peak_subset": recommended_peak_subset(ctype, region),
            "followup_note": followup_note(ctype, region),
        })
    return out


def build_aggregated_candidates(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in rows:
        abs_delta = safe_float(row.get("dominant_shift_abs_delta_percent", 0.0), 0.0)
        delta = safe_float(row.get("dominant_shift_delta_percent", 0.0), 0.0)
        region = normalize_text(row.get("dominant_shift_region"))
        ctype = normalize_text(row.get("comparison_type"))

        score = (
            base_priority_score(abs_delta)
            + scope_bonus("aggregated")
            + comparison_type_bonus(ctype)
            + region_bonus(region)
        )

        left_target = normalize_text(row.get("left_target_label"))
        left_biosample = normalize_text(row.get("left_biosample_term_name"))
        right_target = normalize_text(row.get("right_target_label"))
        right_biosample = normalize_text(row.get("right_biosample_term_name"))

        comparison_id = (
            f"aggregated__{left_target}__{left_biosample}"
            f"__VS__{right_target}__{right_biosample}"
        )

        out.append({
            "comparison_scope": "aggregated",
            "comparison_type": ctype,
            "comparison_id": comparison_id,
            "left_target_label": left_target,
            "left_biosample_term_name": left_biosample,
            "left_experiment_accession": "",
            "left_selected_file_accession": "",
            "right_target_label": right_target,
            "right_biosample_term_name": right_biosample,
            "right_experiment_accession": "",
            "right_selected_file_accession": "",
            "dominant_shift_region": region,
            "dominant_shift_delta_percent": round(delta, 4),
            "dominant_shift_abs_delta_percent": round(abs_delta, 4),
            "left_dominant_region": normalize_text(row.get("left_dominant_region")),
            "left_dominant_region_percent": round(
                safe_float(row.get("left_dominant_region_percent", 0.0), 0.0), 4
            ),
            "right_dominant_region": normalize_text(row.get("right_dominant_region")),
            "right_dominant_region_percent": round(
                safe_float(row.get("right_dominant_region_percent", 0.0), 0.0), 4
            ),
            "priority_score": round(score, 4),
            "recommended_followup_focus": recommended_followup_focus(ctype, region),
            "recommended_peak_subset": recommended_peak_subset(ctype, region),
            "followup_note": followup_note(ctype, region),
        })
    return out


def assign_priority_levels(
    rows: List[Dict[str, Any]],
    top_priority1: int,
    top_priority2: int,
) -> List[Dict[str, Any]]:
    rows = sorted(
        rows,
        key=lambda r: (
            -safe_float(r.get("priority_score", 0.0), 0.0),
            -safe_float(r.get("dominant_shift_abs_delta_percent", 0.0), 0.0),
            normalize_text(r.get("comparison_scope")),
            normalize_text(r.get("comparison_id")),
        ),
    )

    out: List[Dict[str, Any]] = []
    for idx, row in enumerate(rows, start=1):
        x = dict(row)
        x["priority_rank"] = idx
        if idx <= top_priority1:
            x["priority_level"] = "priority_1"
        elif idx <= top_priority1 + top_priority2:
            x["priority_level"] = "priority_2"
        else:
            x["priority_level"] = "priority_3"
        out.append(x)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Define exploratory priority comparisons from current comparison summary tables."
    )
    parser.add_argument(
        "--sample-summary",
        default=None,
        help="Override sample_pairwise_region_differences_summary.tsv",
    )
    parser.add_argument(
        "--aggregated-summary",
        default=None,
        help="Override aggregated_pairwise_region_differences_summary.tsv",
    )
    parser.add_argument(
        "--top-priority1",
        type=int,
        default=2,
        help="Top N rows to label as priority_1 (default: 2)",
    )
    parser.add_argument(
        "--top-priority2",
        type=int,
        default=4,
        help="Next M rows to label as priority_2 (default: 4)",
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

    comparison_dir = root / "processed" / "peaks" / "annotations" / "region_distribution" / "comparisons"
    sample_summary_path = (
        Path(args.sample_summary).resolve()
        if args.sample_summary
        else comparison_dir / "sample_pairwise_region_differences_summary.tsv"
    )
    aggregated_summary_path = (
        Path(args.aggregated_summary).resolve()
        if args.aggregated_summary
        else comparison_dir / "aggregated_pairwise_region_differences_summary.tsv"
    )

    out_path = root / "results" / "tables" / "exploratory_priority_comparisons.tsv"
    log_path = root / "logs" / "define_priority_comparisons.log"

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Sample summary input: {sample_summary_path}", log_path)
    log(f"[INFO] Aggregated summary input: {aggregated_summary_path}", log_path)

    sample_rows = read_tsv_rows(sample_summary_path)
    aggregated_rows = read_tsv_rows(aggregated_summary_path)

    if args.targets:
        target_filter = {x.strip().lower() for x in args.targets if x.strip()}
        sample_rows = [
            r for r in sample_rows
            if normalize_text(r.get("left_target_label")).lower() in target_filter
            and normalize_text(r.get("right_target_label")).lower() in target_filter
        ]
        aggregated_rows = [
            r for r in aggregated_rows
            if normalize_text(r.get("left_target_label")).lower() in target_filter
            and normalize_text(r.get("right_target_label")).lower() in target_filter
        ]

    if args.biosamples:
        biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()}
        sample_rows = [
            r for r in sample_rows
            if normalize_text(r.get("left_biosample_term_name")).lower() in biosample_filter
            and normalize_text(r.get("right_biosample_term_name")).lower() in biosample_filter
        ]
        aggregated_rows = [
            r for r in aggregated_rows
            if normalize_text(r.get("left_biosample_term_name")).lower() in biosample_filter
            and normalize_text(r.get("right_biosample_term_name")).lower() in biosample_filter
        ]

    log(f"[INFO] Sample comparison rows: {len(sample_rows)}", log_path)
    log(f"[INFO] Aggregated comparison rows: {len(aggregated_rows)}", log_path)

    candidates = build_sample_candidates(sample_rows) + build_aggregated_candidates(aggregated_rows)
    ranked = assign_priority_levels(
        rows=candidates,
        top_priority1=args.top_priority1,
        top_priority2=args.top_priority2,
    )

    write_tsv(out_path, ranked, OUTPUT_COLUMNS)

    log("", log_path)
    log("[DONE] Priority comparison definition completed.", log_path)
    log(f"[DONE] Output TSV: {out_path}", log_path)
    log("", log_path)

    for row in ranked:
        if normalize_text(row.get("priority_level")) != "priority_1":
            continue
        log(
            "  " + " | ".join([
                normalize_text(row.get("priority_level")),
                normalize_text(row.get("comparison_scope")),
                normalize_text(row.get("comparison_type")),
                f"{normalize_text(row.get('left_target_label'))}/{normalize_text(row.get('left_biosample_term_name'))}/{normalize_text(row.get('left_experiment_accession'))}",
                f"{normalize_text(row.get('right_target_label'))}/{normalize_text(row.get('right_biosample_term_name'))}/{normalize_text(row.get('right_experiment_accession'))}",
                f"region={normalize_text(row.get('dominant_shift_region'))}",
                f"abs_delta_pct={normalize_text(row.get('dominant_shift_abs_delta_percent'))}",
                f"score={normalize_text(row.get('priority_score'))}",
            ]),
            log_path,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())