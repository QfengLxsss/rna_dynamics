#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarize peak region annotation distributions.

Input
-----
- processed/peaks/annotations/annotation_summary.tsv

Outputs
-------
- processed/peaks/annotations/region_distribution/sample_region_distribution_long.tsv
- processed/peaks/annotations/region_distribution/sample_region_distribution_wide.tsv
- processed/peaks/annotations/region_distribution/target_biosample_region_distribution.tsv
- processed/peaks/annotations/region_distribution/target_region_distribution.tsv
- processed/peaks/annotations/region_distribution/biosample_region_distribution.tsv
- processed/peaks/annotations/region_distribution/region_distribution_report.tsv
- logs/summarize_peak_region_distributions.log
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


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


LONG_COLUMNS = [
    "group_type",
    "target_label",
    "biosample_term_name",
    "experiment_accession",
    "selected_file_accession",
    "region_class",
    "peak_count",
    "peak_fraction",
    "peak_percent",
]

WIDE_COLUMNS_BASE = [
    "target_label",
    "biosample_term_name",
    "experiment_accession",
    "selected_file_accession",
    "n_peaks",
]

REPORT_COLUMNS = [
    "group_type",
    "target_label",
    "biosample_term_name",
    "experiment_accession",
    "selected_file_accession",
    "n_peaks",
    "dominant_region",
    "dominant_region_count",
    "dominant_region_fraction",
    "dominant_region_percent",
    "non_intergenic_fraction",
    "non_intergenic_percent",
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


def make_wide_columns() -> List[str]:
    cols = list(WIDE_COLUMNS_BASE)
    for region in REGIONS:
        cols.append(f"n_{region}")
    for region in REGIONS:
        cols.append(f"frac_{region}")
    for region in REGIONS:
        cols.append(f"pct_{region}")
    cols.append("non_intergenic_fraction")
    cols.append("non_intergenic_percent")
    return cols


def extract_counts(row: Dict[str, str]) -> Dict[str, int]:
    counts = {}
    for region in REGIONS:
        counts[region] = safe_int(row.get(f"n_{region}", 0), 0)
    return counts


def compute_non_intergenic_fraction(counts: Dict[str, int], n_peaks: int) -> float:
    if n_peaks <= 0:
        return 0.0
    non_intergenic = sum(counts[r] for r in REGIONS if r != "intergenic")
    return non_intergenic / n_peaks


def dominant_region(counts: Dict[str, int]) -> Tuple[str, int]:
    best_region = "NA"
    best_count = -1
    for region in REGIONS:
        c = counts.get(region, 0)
        if c > best_count:
            best_region = region
            best_count = c
    return best_region, best_count


def sample_rows_to_long(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in rows:
        target = normalize_text(row.get("target_label"))
        biosample = normalize_text(row.get("biosample_term_name"))
        exp = normalize_text(row.get("experiment_accession"))
        file_acc = normalize_text(row.get("selected_file_accession"))
        n_peaks = safe_int(row.get("n_peaks", 0), 0)
        counts = extract_counts(row)

        for region in REGIONS:
            c = counts[region]
            frac = (c / n_peaks) if n_peaks > 0 else 0.0
            out.append({
                "group_type": "sample",
                "target_label": target,
                "biosample_term_name": biosample,
                "experiment_accession": exp,
                "selected_file_accession": file_acc,
                "region_class": region,
                "peak_count": c,
                "peak_fraction": round(frac, 6),
                "peak_percent": round(frac * 100, 4),
            })
    return out


def sample_rows_to_wide(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in rows:
        target = normalize_text(row.get("target_label"))
        biosample = normalize_text(row.get("biosample_term_name"))
        exp = normalize_text(row.get("experiment_accession"))
        file_acc = normalize_text(row.get("selected_file_accession"))
        n_peaks = safe_int(row.get("n_peaks", 0), 0)
        counts = extract_counts(row)

        x = {
            "target_label": target,
            "biosample_term_name": biosample,
            "experiment_accession": exp,
            "selected_file_accession": file_acc,
            "n_peaks": n_peaks,
        }

        for region in REGIONS:
            frac = (counts[region] / n_peaks) if n_peaks > 0 else 0.0
            x[f"n_{region}"] = counts[region]
            x[f"frac_{region}"] = round(frac, 6)
            x[f"pct_{region}"] = round(frac * 100, 4)

        non_inter_frac = compute_non_intergenic_fraction(counts, n_peaks)
        x["non_intergenic_fraction"] = round(non_inter_frac, 6)
        x["non_intergenic_percent"] = round(non_inter_frac * 100, 4)
        out.append(x)

    out = sorted(
        out,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
            normalize_text(r.get("selected_file_accession")),
        ),
    )
    return out


def aggregate_rows(
    rows: List[Dict[str, str]],
    group_type: str,
) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for row in rows:
        target = normalize_text(row.get("target_label"))
        biosample = normalize_text(row.get("biosample_term_name"))

        if group_type == "target_biosample":
            key = (target, biosample)
            out_target = target
            out_biosample = biosample
        elif group_type == "target":
            key = (target, "")
            out_target = target
            out_biosample = ""
        elif group_type == "biosample":
            key = ("", biosample)
            out_target = ""
            out_biosample = biosample
        else:
            raise ValueError(f"Unsupported group_type: {group_type}")

        if key not in grouped:
            grouped[key] = {
                "target_label": out_target,
                "biosample_term_name": out_biosample,
                "experiment_accession": "",
                "selected_file_accession": "",
                "n_peaks": 0,
            }
            for region in REGIONS:
                grouped[key][f"n_{region}"] = 0

        grouped[key]["n_peaks"] += safe_int(row.get("n_peaks", 0), 0)
        counts = extract_counts(row)
        for region in REGIONS:
            grouped[key][f"n_{region}"] += counts[region]

    out: List[Dict[str, Any]] = []
    for _, g in grouped.items():
        n_peaks = safe_int(g.get("n_peaks", 0), 0)
        for region in REGIONS:
            c = safe_int(g.get(f"n_{region}", 0), 0)
            frac = (c / n_peaks) if n_peaks > 0 else 0.0
            g[f"frac_{region}"] = round(frac, 6)
            g[f"pct_{region}"] = round(frac * 100, 4)

        counts = {region: safe_int(g.get(f"n_{region}", 0), 0) for region in REGIONS}
        non_inter_frac = compute_non_intergenic_fraction(counts, n_peaks)
        g["non_intergenic_fraction"] = round(non_inter_frac, 6)
        g["non_intergenic_percent"] = round(non_inter_frac * 100, 4)
        out.append(g)

    out = sorted(
        out,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
        ),
    )
    return out


def make_report_rows(
    sample_wide_rows: List[Dict[str, Any]],
    tb_rows: List[Dict[str, Any]],
    target_rows: List[Dict[str, Any]],
    biosample_rows: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    report: List[Dict[str, Any]] = []

    def add_rows(rows: List[Dict[str, Any]], group_type: str) -> None:
        for row in rows:
            counts = {region: safe_int(row.get(f"n_{region}", 0), 0) for region in REGIONS}
            n_peaks = safe_int(row.get("n_peaks", 0), 0)
            dom_region, dom_count = dominant_region(counts)
            dom_frac = (dom_count / n_peaks) if n_peaks > 0 else 0.0

            report.append({
                "group_type": group_type,
                "target_label": normalize_text(row.get("target_label")),
                "biosample_term_name": normalize_text(row.get("biosample_term_name")),
                "experiment_accession": normalize_text(row.get("experiment_accession")),
                "selected_file_accession": normalize_text(row.get("selected_file_accession")),
                "n_peaks": n_peaks,
                "dominant_region": dom_region,
                "dominant_region_count": dom_count,
                "dominant_region_fraction": round(dom_frac, 6),
                "dominant_region_percent": round(dom_frac * 100, 4),
                "non_intergenic_fraction": round(
                    safe_float(row.get("non_intergenic_fraction", 0.0), 0.0), 6
                ),
                "non_intergenic_percent": round(
                    safe_float(row.get("non_intergenic_percent", 0.0), 0.0), 4
                ),
            })

    add_rows(sample_wide_rows, "sample")
    add_rows(tb_rows, "target_biosample")
    add_rows(target_rows, "target")
    add_rows(biosample_rows, "biosample")

    report = sorted(
        report,
        key=lambda r: (
            normalize_text(r.get("group_type")),
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
            normalize_text(r.get("selected_file_accession")),
        ),
    )
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize peak region annotation distributions."
    )
    parser.add_argument(
        "--annotation-summary",
        default=None,
        help="Override processed/peaks/annotations/annotation_summary.tsv",
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

    annotation_summary = (
        Path(args.annotation_summary).resolve()
        if args.annotation_summary
        else root / "processed" / "peaks" / "annotations" / "annotation_summary.tsv"
    )

    out_dir = root / "processed" / "peaks" / "annotations" / "region_distribution"
    log_path = root / "logs" / "summarize_peak_region_distributions.log"

    long_out = out_dir / "sample_region_distribution_long.tsv"
    wide_out = out_dir / "sample_region_distribution_wide.tsv"
    tb_out = out_dir / "target_biosample_region_distribution.tsv"
    target_out = out_dir / "target_region_distribution.tsv"
    biosample_out = out_dir / "biosample_region_distribution.tsv"
    report_out = out_dir / "region_distribution_report.tsv"

    rows = read_tsv_rows(annotation_summary)

    if args.targets:
        target_filter = {x.strip().lower() for x in args.targets if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("target_label")).lower() in target_filter]

    if args.biosamples:
        biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter]

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Annotation summary input: {annotation_summary}", log_path)
    log(f"[INFO] Input rows: {len(rows)}", log_path)

    if not rows:
        log("[WARN] No rows available after filtering.", log_path)
        return 0

    sample_long_rows = sample_rows_to_long(rows)
    sample_wide_rows = sample_rows_to_wide(rows)
    tb_rows = aggregate_rows(rows, "target_biosample")
    target_rows = aggregate_rows(rows, "target")
    biosample_rows = aggregate_rows(rows, "biosample")
    report_rows = make_report_rows(sample_wide_rows, tb_rows, target_rows, biosample_rows)

    write_tsv(long_out, sample_long_rows, LONG_COLUMNS)
    write_tsv(wide_out, sample_wide_rows, make_wide_columns())
    write_tsv(tb_out, tb_rows, make_wide_columns())
    write_tsv(target_out, target_rows, make_wide_columns())
    write_tsv(biosample_out, biosample_rows, make_wide_columns())
    write_tsv(report_out, report_rows, REPORT_COLUMNS)

    log("", log_path)
    log("[DONE] Peak region distribution summarization completed.", log_path)
    log(f"[DONE] Long TSV: {long_out}", log_path)
    log(f"[DONE] Sample wide TSV: {wide_out}", log_path)
    log(f"[DONE] Target+biosample TSV: {tb_out}", log_path)
    log(f"[DONE] Target TSV: {target_out}", log_path)
    log(f"[DONE] Biosample TSV: {biosample_out}", log_path)
    log(f"[DONE] Report TSV: {report_out}", log_path)
    log("", log_path)

    for row in report_rows:
        if normalize_text(row.get("group_type")) != "sample":
            continue
        log(
            "  "
            + " | ".join(
                [
                    normalize_text(row.get("target_label")),
                    normalize_text(row.get("biosample_term_name")),
                    normalize_text(row.get("experiment_accession")),
                    normalize_text(row.get("selected_file_accession")),
                    f"dominant={normalize_text(row.get('dominant_region'))}",
                    f"dominant_pct={normalize_text(row.get('dominant_region_percent'))}",
                    f"non_intergenic_pct={normalize_text(row.get('non_intergenic_percent'))}",
                ]
            ),
            log_path,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())