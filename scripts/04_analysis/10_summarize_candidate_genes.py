#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarize candidate genes from per-set gene mapping results.

Inputs
------
- results/tables/gene_mapping/peak_set_gene_mapping_manifest.tsv
- results/tables/gene_mapping/peak_set_gene_mapping_summary.tsv

Per-set required file
---------------------
- output_gene_summary_tsv

Outputs
-------
- results/tables/candidate_genes/candidate_gene_overview.tsv
- results/tables/candidate_genes/top_candidate_genes_by_peak_set.tsv
- results/tables/candidate_genes/candidate_gene_recurrence.tsv
- results/tables/candidate_genes/candidate_gene_set_summary.tsv
- logs/summarize_candidate_genes.log

Purpose
-------
Turn multiple per-set gene summary tables into compact exploratory candidate gene tables.

Main ideas
----------
1. For each peak subset, keep the top-ranked genes.
2. Across all processed subsets, measure recurrence:
   - in how many sets a gene appears
   - in which regions / set types it appears
   - total peak hits accumulated
3. Produce a compact overview table for round-1 exploratory interpretation.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


OVERVIEW_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "n_peak_rows",
    "n_unique_genes",
    "top_gene_name",
    "top_gene_id",
    "top_gene_peak_hits",
    "top5_gene_names",
    "top10_gene_names",
]

TOP_BY_SET_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "gene_rank",
    "gene_id",
    "gene_name",
    "peak_hit_count",
    "unique_peak_count",
    "transcript_count",
    "transcript_ids",
    "peak_hit_fraction_within_set",
]

RECURRENCE_COLUMNS = [
    "gene_id",
    "gene_name",
    "n_sets_present",
    "set_names",
    "region_classes",
    "comparison_ids",
    "total_peak_hits",
    "total_unique_peaks",
    "max_peak_hits_in_one_set",
    "best_set_name",
    "best_region_class",
    "best_comparison_id",
]

SET_SUMMARY_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "n_peak_rows",
    "n_peak_gene_links",
    "n_unique_genes",
    "top_gene_name",
    "top_gene_id",
    "top_gene_peak_hits",
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


def normalize_text(x: Any) -> str:
    if x is None:
        return ""
    return str(x).strip()


def safe_int(x: Any, default: int = 0) -> int:
    try:
        return int(str(x).strip())
    except Exception:
        return default


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


def top_gene_names(rows: List[Dict[str, str]], top_n: int) -> str:
    names = []
    for row in rows[:top_n]:
        gname = normalize_text(row.get("gene_name"))
        gid = normalize_text(row.get("gene_id"))
        if gname:
            names.append(gname)
        elif gid:
            names.append(gid)
    return ";".join(names)


def gene_key(row: Dict[str, str]) -> Tuple[str, str]:
    return (
        normalize_text(row.get("gene_id")),
        normalize_text(row.get("gene_name")),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize candidate genes from per-set gene mapping results."
    )
    parser.add_argument(
        "--gene-mapping-manifest",
        default=None,
        help="Override results/tables/gene_mapping/peak_set_gene_mapping_manifest.tsv",
    )
    parser.add_argument(
        "--gene-mapping-summary",
        default=None,
        help="Override results/tables/gene_mapping/peak_set_gene_mapping_summary.tsv",
    )
    parser.add_argument(
        "--regions",
        nargs="*",
        default=None,
        help="Restrict to these region classes, e.g. intron CDS",
    )
    parser.add_argument(
        "--set-names",
        nargs="*",
        default=None,
        help="Restrict to these set names",
    )
    parser.add_argument(
        "--comparison-ids",
        nargs="*",
        default=None,
        help="Restrict to these comparison IDs",
    )
    parser.add_argument(
        "--top-n-per-set",
        type=int,
        default=20,
        help="Top N genes to keep for each set in top_candidate_genes_by_peak_set.tsv (default: 20)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    gene_mapping_dir = root / "results" / "tables" / "gene_mapping"
    manifest_path = (
        Path(args.gene_mapping_manifest).resolve()
        if args.gene_mapping_manifest
        else gene_mapping_dir / "peak_set_gene_mapping_manifest.tsv"
    )
    summary_path = (
        Path(args.gene_mapping_summary).resolve()
        if args.gene_mapping_summary
        else gene_mapping_dir / "peak_set_gene_mapping_summary.tsv"
    )

    out_dir = root / "results" / "tables" / "candidate_genes"
    overview_out = out_dir / "candidate_gene_overview.tsv"
    top_by_set_out = out_dir / "top_candidate_genes_by_peak_set.tsv"
    recurrence_out = out_dir / "candidate_gene_recurrence.tsv"
    set_summary_out = out_dir / "candidate_gene_set_summary.tsv"
    log_path = root / "logs" / "summarize_candidate_genes.log"

    manifest_rows = read_tsv_rows(manifest_path)
    summary_rows = read_tsv_rows(summary_path)

    summary_idx: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    for row in summary_rows:
        key = (
            normalize_text(row.get("comparison_id")),
            normalize_text(row.get("set_name")),
            normalize_text(row.get("region_class")),
        )
        summary_idx[key] = row

    if args.regions:
        allowed_regions = {normalize_text(x) for x in args.regions if normalize_text(x)}
        manifest_rows = [r for r in manifest_rows if normalize_text(r.get("region_class")) in allowed_regions]

    if args.set_names:
        allowed_sets = {normalize_text(x) for x in args.set_names if normalize_text(x)}
        manifest_rows = [r for r in manifest_rows if normalize_text(r.get("set_name")) in allowed_sets]

    if args.comparison_ids:
        allowed_cmp = {normalize_text(x) for x in args.comparison_ids if normalize_text(x)}
        manifest_rows = [r for r in manifest_rows if normalize_text(r.get("comparison_id")) in allowed_cmp]

    manifest_rows = [r for r in manifest_rows if normalize_text(r.get("status")) == "ok"]

    manifest_rows = sorted(
        manifest_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_id")),
            normalize_text(r.get("set_name")),
            normalize_text(r.get("region_class")),
        ),
    )

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Gene mapping manifest TSV: {manifest_path}", log_path)
    log(f"[INFO] Gene mapping summary TSV: {summary_path}", log_path)
    log(f"[INFO] Eligible peak sets: {len(manifest_rows)}", log_path)

    overview_rows: List[Dict[str, Any]] = []
    top_by_set_rows: List[Dict[str, Any]] = []
    set_summary_rows: List[Dict[str, Any]] = []
    recurrence_grouped: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for i, row in enumerate(manifest_rows, start=1):
        comparison_id = normalize_text(row.get("comparison_id"))
        set_name = normalize_text(row.get("set_name"))
        region_class = normalize_text(row.get("region_class"))
        gene_summary_rel = normalize_text(row.get("output_gene_summary_tsv"))
        gene_summary_path = root / gene_summary_rel

        log(
            f"[INFO] [{i}/{len(manifest_rows)}] Summarizing candidates for {comparison_id} | {set_name} | {region_class}",
            log_path,
        )

        if not gene_summary_path.exists():
            continue

        gene_rows = read_tsv_rows(gene_summary_path)
        key = (comparison_id, set_name, region_class)
        summary_row = summary_idx.get(key, {})

        n_peak_rows = safe_int(summary_row.get("n_peak_rows", 0), 0)
        n_unique_genes = len(gene_rows)

        top_gene_name = ""
        top_gene_id = ""
        top_gene_peak_hits = 0
        if gene_rows:
            top_gene_name = normalize_text(gene_rows[0].get("gene_name"))
            top_gene_id = normalize_text(gene_rows[0].get("gene_id"))
            top_gene_peak_hits = safe_int(gene_rows[0].get("peak_hit_count", 0), 0)

        overview_rows.append({
            "comparison_id": comparison_id,
            "set_name": set_name,
            "region_class": region_class,
            "n_peak_rows": n_peak_rows,
            "n_unique_genes": n_unique_genes,
            "top_gene_name": top_gene_name,
            "top_gene_id": top_gene_id,
            "top_gene_peak_hits": top_gene_peak_hits,
            "top5_gene_names": top_gene_names(gene_rows, 5),
            "top10_gene_names": top_gene_names(gene_rows, 10),
        })

        set_summary_rows.append({
            "comparison_id": comparison_id,
            "set_name": set_name,
            "region_class": region_class,
            "n_peak_rows": safe_int(summary_row.get("n_peak_rows", 0), 0),
            "n_peak_gene_links": safe_int(summary_row.get("n_peak_gene_links", 0), 0),
            "n_unique_genes": safe_int(summary_row.get("n_unique_genes", 0), 0),
            "top_gene_name": top_gene_name,
            "top_gene_id": top_gene_id,
            "top_gene_peak_hits": top_gene_peak_hits,
            "status": normalize_text(summary_row.get("status")) or "ok",
            "note": normalize_text(summary_row.get("note")) or "candidate_gene_summary_completed",
        })

        for rank, gene_row in enumerate(gene_rows[: args.top_n_per_set], start=1):
            peak_hit_count = safe_int(gene_row.get("peak_hit_count", 0), 0)
            peak_hit_fraction = (peak_hit_count / n_peak_rows) if n_peak_rows > 0 else 0.0

            top_by_set_rows.append({
                "comparison_id": comparison_id,
                "set_name": set_name,
                "region_class": region_class,
                "gene_rank": rank,
                "gene_id": normalize_text(gene_row.get("gene_id")),
                "gene_name": normalize_text(gene_row.get("gene_name")),
                "peak_hit_count": peak_hit_count,
                "unique_peak_count": safe_int(gene_row.get("unique_peak_count", 0), 0),
                "transcript_count": safe_int(gene_row.get("transcript_count", 0), 0),
                "transcript_ids": normalize_text(gene_row.get("transcript_ids")),
                "peak_hit_fraction_within_set": round(peak_hit_fraction, 6),
            })

        for gene_row in gene_rows:
            gkey = gene_key(gene_row)
            if not normalize_text(gkey[0]) and not normalize_text(gkey[1]):
                continue

            if gkey not in recurrence_grouped:
                recurrence_grouped[gkey] = {
                    "gene_id": gkey[0],
                    "gene_name": gkey[1],
                    "set_keys": set(),
                    "set_names": set(),
                    "region_classes": set(),
                    "comparison_ids": set(),
                    "total_peak_hits": 0,
                    "total_unique_peaks": 0,
                    "max_peak_hits_in_one_set": 0,
                    "best_set_name": "",
                    "best_region_class": "",
                    "best_comparison_id": "",
                }

            rec = recurrence_grouped[gkey]
            set_key = (comparison_id, set_name, region_class)
            rec["set_keys"].add(set_key)
            rec["set_names"].add(set_name)
            rec["region_classes"].add(region_class)
            rec["comparison_ids"].add(comparison_id)

            peak_hits = safe_int(gene_row.get("peak_hit_count", 0), 0)
            unique_peaks = safe_int(gene_row.get("unique_peak_count", 0), 0)

            rec["total_peak_hits"] += peak_hits
            rec["total_unique_peaks"] += unique_peaks

            if peak_hits > rec["max_peak_hits_in_one_set"]:
                rec["max_peak_hits_in_one_set"] = peak_hits
                rec["best_set_name"] = set_name
                rec["best_region_class"] = region_class
                rec["best_comparison_id"] = comparison_id

    recurrence_rows: List[Dict[str, Any]] = []
    for _, rec in recurrence_grouped.items():
        recurrence_rows.append({
            "gene_id": rec["gene_id"],
            "gene_name": rec["gene_name"],
            "n_sets_present": len(rec["set_keys"]),
            "set_names": ";".join(sorted(rec["set_names"])),
            "region_classes": ";".join(sorted(rec["region_classes"])),
            "comparison_ids": ";".join(sorted(rec["comparison_ids"])),
            "total_peak_hits": rec["total_peak_hits"],
            "total_unique_peaks": rec["total_unique_peaks"],
            "max_peak_hits_in_one_set": rec["max_peak_hits_in_one_set"],
            "best_set_name": rec["best_set_name"],
            "best_region_class": rec["best_region_class"],
            "best_comparison_id": rec["best_comparison_id"],
        })

    overview_rows = sorted(
        overview_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_id")),
            normalize_text(r.get("set_name")),
            normalize_text(r.get("region_class")),
        ),
    )

    top_by_set_rows = sorted(
        top_by_set_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_id")),
            normalize_text(r.get("set_name")),
            normalize_text(r.get("region_class")),
            safe_int(r.get("gene_rank"), 0),
        ),
    )

    recurrence_rows = sorted(
        recurrence_rows,
        key=lambda r: (
            -safe_int(r.get("n_sets_present"), 0),
            -safe_int(r.get("total_peak_hits"), 0),
            normalize_text(r.get("gene_name")),
            normalize_text(r.get("gene_id")),
        ),
    )

    set_summary_rows = sorted(
        set_summary_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_id")),
            normalize_text(r.get("set_name")),
            normalize_text(r.get("region_class")),
        ),
    )

    write_tsv(overview_out, overview_rows, OVERVIEW_COLUMNS)
    write_tsv(top_by_set_out, top_by_set_rows, TOP_BY_SET_COLUMNS)
    write_tsv(recurrence_out, recurrence_rows, RECURRENCE_COLUMNS)
    write_tsv(set_summary_out, set_summary_rows, SET_SUMMARY_COLUMNS)

    log("", log_path)
    log("[DONE] Candidate gene summarization completed.", log_path)
    log(f"[DONE] Overview TSV: {overview_out}", log_path)
    log(f"[DONE] Top-by-set TSV: {top_by_set_out}", log_path)
    log(f"[DONE] Recurrence TSV: {recurrence_out}", log_path)
    log(f"[DONE] Set summary TSV: {set_summary_out}", log_path)

    top_overview = sorted(
        overview_rows,
        key=lambda r: (
            -safe_int(r.get("top_gene_peak_hits"), 0),
            normalize_text(r.get("comparison_id")),
        ),
    )[:10]

    log("[INFO] Top candidate-gene subsets:", log_path)
    for r in top_overview:
        log(
            "  " + " | ".join([
                normalize_text(r.get("comparison_id")),
                normalize_text(r.get("set_name")),
                normalize_text(r.get("region_class")),
                f"top_gene={normalize_text(r.get('top_gene_name')) or normalize_text(r.get('top_gene_id'))}",
                f"top_hits={normalize_text(r.get('top_gene_peak_hits'))}",
            ]),
            log_path,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())