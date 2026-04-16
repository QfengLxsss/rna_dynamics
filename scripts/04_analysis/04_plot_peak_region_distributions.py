#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot peak region distribution summaries.

Inputs
------
- processed/peaks/annotations/region_distribution/sample_region_distribution_wide.tsv
- processed/peaks/annotations/region_distribution/target_biosample_region_distribution.tsv
- processed/peaks/annotations/region_distribution/comparisons/sample_pairwise_region_differences_summary.tsv
- processed/peaks/annotations/region_distribution/comparisons/aggregated_pairwise_region_differences_summary.tsv

Outputs
-------
- processed/peaks/annotations/region_distribution/plots/sample_region_distribution_stacked.png
- processed/peaks/annotations/region_distribution/plots/target_biosample_region_distribution_stacked.png
- processed/peaks/annotations/region_distribution/plots/sample_region_distribution_heatmap.png
- processed/peaks/annotations/region_distribution/plots/target_biosample_region_distribution_heatmap.png
- processed/peaks/annotations/region_distribution/plots/sample_top_region_shifts.png
- processed/peaks/annotations/region_distribution/plots/aggregated_top_region_shifts.png
- logs/plot_peak_region_distributions.log
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
matplotlib.use("Agg")  # 关键：服务器无图形界面时强制使用非交互后端
import matplotlib.pyplot as plt


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


def normalize_text(x: Any) -> str:
    if x is None:
        return ""
    return str(x).strip()


def safe_float(x: Any, default: float = 0.0) -> float:
    try:
        return float(str(x).strip())
    except Exception:
        return default


def filter_standard_rows(
    rows: List[Dict[str, str]],
    targets: List[str] | None,
    biosamples: List[str] | None,
) -> List[Dict[str, str]]:
    if targets:
        target_filter = {x.strip().lower() for x in targets if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("target_label")).lower() in target_filter]
    if biosamples:
        biosample_filter = {x.strip().lower() for x in biosamples if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter]
    return rows


def filter_comparison_rows(
    rows: List[Dict[str, str]],
    targets: List[str] | None,
    biosamples: List[str] | None,
) -> List[Dict[str, str]]:
    if targets:
        target_filter = {x.strip().lower() for x in targets if x.strip()}
        rows = [
            r for r in rows
            if normalize_text(r.get("left_target_label")).lower() in target_filter
            and normalize_text(r.get("right_target_label")).lower() in target_filter
        ]
    if biosamples:
        biosample_filter = {x.strip().lower() for x in biosamples if x.strip()}
        rows = [
            r for r in rows
            if normalize_text(r.get("left_biosample_term_name")).lower() in biosample_filter
            and normalize_text(r.get("right_biosample_term_name")).lower() in biosample_filter
        ]
    return rows


def sort_sample_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return sorted(
        rows,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
            normalize_text(r.get("selected_file_accession")),
        ),
    )


def sort_target_biosample_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return sorted(
        rows,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
        ),
    )


def sample_label(row: Dict[str, str]) -> str:
    return (
        f"{normalize_text(row.get('target_label'))}\n"
        f"{normalize_text(row.get('biosample_term_name'))}\n"
        f"{normalize_text(row.get('experiment_accession'))}"
    )


def target_biosample_label(row: Dict[str, str]) -> str:
    return (
        f"{normalize_text(row.get('target_label'))}\n"
        f"{normalize_text(row.get('biosample_term_name'))}"
    )


def region_percent_matrix(rows: List[Dict[str, str]]) -> List[List[float]]:
    matrix = []
    for row in rows:
        matrix.append([safe_float(row.get(f"pct_{region}", 0.0), 0.0) for region in REGIONS])
    return matrix


def plot_stacked_bar(
    rows: List[Dict[str, str]],
    label_func,
    title: str,
    out_path: Path,
    log_path: Path,
) -> None:
    if not rows:
        log(f"[WARN] No rows for stacked bar: {out_path.name}", log_path)
        return

    labels = [label_func(r) for r in rows]
    x = list(range(len(rows)))

    fig_width = max(10, min(28, 1.25 * len(rows) + 4))
    fig_height = 7

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    bottoms = [0.0] * len(rows)

    for region in REGIONS:
        vals = [safe_float(r.get(f"pct_{region}", 0.0), 0.0) for r in rows]
        ax.bar(x, vals, bottom=bottoms, label=region)
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    ax.set_title(title)
    ax.set_ylabel("Percent of peaks")
    ax.set_xlabel("Samples")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylim(0, 100)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"[DONE] Wrote plot: {out_path}", log_path)


def plot_heatmap(
    rows: List[Dict[str, str]],
    label_func,
    title: str,
    out_path: Path,
    log_path: Path,
) -> None:
    if not rows:
        log(f"[WARN] No rows for heatmap: {out_path.name}", log_path)
        return

    labels = [label_func(r) for r in rows]
    matrix = region_percent_matrix(rows)

    fig_width = max(10, min(26, 0.8 * len(rows) + 6))
    fig_height = max(5, min(12, 0.6 * len(REGIONS) + 3))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    im = ax.imshow(matrix, aspect="auto")

    ax.set_title(title)
    ax.set_xticks(range(len(REGIONS)))
    ax.set_xticklabels(REGIONS, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Percent of peaks")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"[DONE] Wrote plot: {out_path}", log_path)


def truncate_label(text: str, max_len: int = 70) -> str:
    text = normalize_text(text)
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def plot_top_shifts(
    rows: List[Dict[str, str]],
    label_builder,
    title: str,
    out_path: Path,
    log_path: Path,
    top_n: int = 10,
) -> None:
    if not rows:
        log(f"[WARN] No rows for top shifts plot: {out_path.name}", log_path)
        return

    rows = sorted(
        rows,
        key=lambda r: -safe_float(r.get("dominant_shift_abs_delta_percent", 0.0), 0.0),
    )[:top_n]

    labels = [truncate_label(label_builder(r), 80) for r in rows]
    vals = [safe_float(r.get("dominant_shift_abs_delta_percent", 0.0), 0.0) for r in rows]

    fig_height = max(5, min(12, 0.6 * len(rows) + 2))
    fig, ax = plt.subplots(figsize=(12, fig_height))
    y = list(range(len(rows)))

    ax.barh(y, vals)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Absolute delta percent")
    ax.set_title(title)

    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:.2f}", va="center")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"[DONE] Wrote plot: {out_path}", log_path)


def sample_comparison_label(row: Dict[str, str]) -> str:
    return (
        f"{normalize_text(row.get('comparison_type'))} | "
        f"{normalize_text(row.get('left_target_label'))}/{normalize_text(row.get('left_biosample_term_name'))}/{normalize_text(row.get('left_experiment_accession'))} "
        f"vs "
        f"{normalize_text(row.get('right_target_label'))}/{normalize_text(row.get('right_biosample_term_name'))}/{normalize_text(row.get('right_experiment_accession'))} "
        f"| {normalize_text(row.get('dominant_shift_region'))}"
    )


def aggregated_comparison_label(row: Dict[str, str]) -> str:
    return (
        f"{normalize_text(row.get('comparison_type'))} | "
        f"{normalize_text(row.get('left_target_label'))}/{normalize_text(row.get('left_biosample_term_name'))} "
        f"vs "
        f"{normalize_text(row.get('right_target_label'))}/{normalize_text(row.get('right_biosample_term_name'))} "
        f"| {normalize_text(row.get('dominant_shift_region'))}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot peak region distribution summaries."
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
        "--sample-comparison-summary",
        default=None,
        help="Override sample_pairwise_region_differences_summary.tsv",
    )
    parser.add_argument(
        "--aggregated-comparison-summary",
        default=None,
        help="Override aggregated_pairwise_region_differences_summary.tsv",
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
        "--top-n",
        type=int,
        default=10,
        help="Top N comparisons to plot (default: 10)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    region_dir = root / "processed" / "peaks" / "annotations" / "region_distribution"
    comparison_dir = region_dir / "comparisons"
    plot_dir = region_dir / "plots"
    log_path = root / "logs" / "plot_peak_region_distributions.log"

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
    sample_comparison_summary_path = (
        Path(args.sample_comparison_summary).resolve()
        if args.sample_comparison_summary
        else comparison_dir / "sample_pairwise_region_differences_summary.tsv"
    )
    aggregated_comparison_summary_path = (
        Path(args.aggregated_comparison_summary).resolve()
        if args.aggregated_comparison_summary
        else comparison_dir / "aggregated_pairwise_region_differences_summary.tsv"
    )

    log(f"[INFO] Project root: {root}", log_path)

    sample_rows = read_tsv_rows(sample_wide_path)
    tb_rows = read_tsv_rows(target_biosample_path)
    sample_cmp_rows = read_tsv_rows(sample_comparison_summary_path)
    agg_cmp_rows = read_tsv_rows(aggregated_comparison_summary_path)

    sample_rows = filter_standard_rows(sample_rows, args.targets, args.biosamples)
    tb_rows = filter_standard_rows(tb_rows, args.targets, args.biosamples)
    sample_cmp_rows = filter_comparison_rows(sample_cmp_rows, args.targets, args.biosamples)
    agg_cmp_rows = filter_comparison_rows(agg_cmp_rows, args.targets, args.biosamples)

    sample_rows = sort_sample_rows(sample_rows)
    tb_rows = sort_target_biosample_rows(tb_rows)

    log(f"[INFO] Sample rows: {len(sample_rows)}", log_path)
    log(f"[INFO] Target+biosample rows: {len(tb_rows)}", log_path)
    log(f"[INFO] Sample comparison rows: {len(sample_cmp_rows)}", log_path)
    log(f"[INFO] Aggregated comparison rows: {len(agg_cmp_rows)}", log_path)

    plot_stacked_bar(
        rows=sample_rows,
        label_func=sample_label,
        title="Sample-level peak region distributions",
        out_path=plot_dir / "sample_region_distribution_stacked.png",
        log_path=log_path,
    )

    plot_stacked_bar(
        rows=tb_rows,
        label_func=target_biosample_label,
        title="Target+biosample peak region distributions",
        out_path=plot_dir / "target_biosample_region_distribution_stacked.png",
        log_path=log_path,
    )

    plot_heatmap(
        rows=sample_rows,
        label_func=sample_label,
        title="Sample-level peak region distribution heatmap",
        out_path=plot_dir / "sample_region_distribution_heatmap.png",
        log_path=log_path,
    )

    plot_heatmap(
        rows=tb_rows,
        label_func=target_biosample_label,
        title="Target+biosample peak region distribution heatmap",
        out_path=plot_dir / "target_biosample_region_distribution_heatmap.png",
        log_path=log_path,
    )

    plot_top_shifts(
        rows=sample_cmp_rows,
        label_builder=sample_comparison_label,
        title="Top sample-level region shifts",
        out_path=plot_dir / "sample_top_region_shifts.png",
        log_path=log_path,
        top_n=args.top_n,
    )

    plot_top_shifts(
        rows=agg_cmp_rows,
        label_builder=aggregated_comparison_label,
        title="Top aggregated region shifts",
        out_path=plot_dir / "aggregated_top_region_shifts.png",
        log_path=log_path,
        top_n=args.top_n,
    )

    log("", log_path)
    log("[DONE] Plotting completed.", log_path)
    log(f"[DONE] Plot directory: {plot_dir}", log_path)

    return 0


if __name__ == "__main__":
    main()