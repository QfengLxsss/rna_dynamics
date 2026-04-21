#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare enrichment-ready functional gene lists from candidate gene tables.

Inputs
------
- results/tables/candidate_genes/top_candidate_genes_by_peak_set.tsv
- results/tables/candidate_genes/candidate_gene_overview.tsv
- results/tables/candidate_genes/candidate_gene_recurrence.tsv
- results/tables/candidate_genes/candidate_gene_set_summary.tsv

Outputs
-------
- results/tables/functional_inputs/functional_gene_list_manifest.tsv
- results/tables/functional_inputs/functional_gene_list_summary.tsv
- results/tables/functional_inputs/per_set/<comparison_id>/<set_name>.<region>.top_genes.tsv
- results/tables/functional_inputs/per_set/<comparison_id>/<set_name>.<region>.gene_symbols.txt
- results/tables/functional_inputs/per_set/<comparison_id>/<set_name>.<region>.gene_ids.txt
- results/tables/functional_inputs/combined/selected_sets_union.gene_symbols.txt
- results/tables/functional_inputs/combined/selected_sets_union.gene_ids.txt
- results/tables/functional_inputs/combined/recurrent_genes_minN.tsv
- logs/prepare_functional_gene_lists.log

Purpose
-------
Convert exploratory candidate gene tables into lightweight, enrichment-ready gene lists.

Key behaviors
-------------
1. For each selected peak set:
   - keep top N genes by rank from top_candidate_genes_by_peak_set.tsv
   - export TSV + gene_symbols.txt + gene_ids.txt
2. Across all selected peak sets:
   - export union gene symbols / IDs
   - export recurrent genes appearing in at least N sets
3. Preserve enough metadata to trace each list back to:
   - comparison_id
   - set_name
   - region_class
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


PER_SET_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "gene_rank",
    "gene_id",
    "gene_name",
    "gene_label",
    "peak_hit_count",
    "unique_peak_count",
    "transcript_count",
    "transcript_ids",
    "peak_hit_fraction_within_set",
]

MANIFEST_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "n_exported_genes",
    "top_gene_name",
    "top_gene_id",
    "top_gene_peak_hits",
    "top_genes_tsv",
    "gene_symbols_txt",
    "gene_ids_txt",
    "status",
    "note",
]

SUMMARY_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "n_exported_genes",
    "top_gene_name",
    "top_gene_id",
    "top_gene_peak_hits",
    "top5_gene_labels",
    "status",
    "note",
]

RECURRENT_COLUMNS = [
    "gene_id",
    "gene_name",
    "gene_label",
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


def safe_float(x: Any, default: float = 0.0) -> float:
    try:
        return float(str(x).strip())
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


def write_txt_list(path: Path, items: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for item in items:
            f.write(item + "\n")


def strip_ensembl_version(gene_id: str) -> str:
    gene_id = normalize_text(gene_id)
    if "." in gene_id:
        return gene_id.split(".", 1)[0]
    return gene_id


def gene_label(gene_name: str, gene_id: str) -> str:
    gname = normalize_text(gene_name)
    gid = strip_ensembl_version(gene_id)
    return gname if gname else gid


def build_key(row: Dict[str, str]) -> Tuple[str, str, str]:
    return (
        normalize_text(row.get("comparison_id")),
        normalize_text(row.get("set_name")),
        normalize_text(row.get("region_class")),
    )


def top5_labels(rows: List[Dict[str, Any]]) -> str:
    labels = []
    for row in rows[:5]:
        label = normalize_text(row.get("gene_label"))
        if label:
            labels.append(label)
    return ";".join(labels)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare functional gene lists from candidate gene tables."
    )
    parser.add_argument(
        "--top-by-set",
        default=None,
        help="Override results/tables/candidate_genes/top_candidate_genes_by_peak_set.tsv",
    )
    parser.add_argument(
        "--candidate-overview",
        default=None,
        help="Override results/tables/candidate_genes/candidate_gene_overview.tsv",
    )
    parser.add_argument(
        "--candidate-recurrence",
        default=None,
        help="Override results/tables/candidate_genes/candidate_gene_recurrence.tsv",
    )
    parser.add_argument(
        "--candidate-set-summary",
        default=None,
        help="Override results/tables/candidate_genes/candidate_gene_set_summary.tsv",
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
        help="Restrict to these set names, e.g. left_specific right_specific",
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
        default=100,
        help="How many genes to export for each selected peak set (default: 100)",
    )
    parser.add_argument(
        "--min-peak-hits",
        type=int,
        default=1,
        help="Minimum peak_hit_count required to keep a gene in per-set exports (default: 1)",
    )
    parser.add_argument(
        "--min-recurrence",
        type=int,
        default=2,
        help="Minimum number of sets a gene must appear in to be exported as recurrent (default: 2)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    cand_dir = root / "results" / "tables" / "candidate_genes"
    top_by_set_path = (
        Path(args.top_by_set).resolve()
        if args.top_by_set
        else cand_dir / "top_candidate_genes_by_peak_set.tsv"
    )
    candidate_overview_path = (
        Path(args.candidate_overview).resolve()
        if args.candidate_overview
        else cand_dir / "candidate_gene_overview.tsv"
    )
    candidate_recurrence_path = (
        Path(args.candidate_recurrence).resolve()
        if args.candidate_recurrence
        else cand_dir / "candidate_gene_recurrence.tsv"
    )
    candidate_set_summary_path = (
        Path(args.candidate_set_summary).resolve()
        if args.candidate_set_summary
        else cand_dir / "candidate_gene_set_summary.tsv"
    )

    out_dir = root / "results" / "tables" / "functional_inputs"
    per_set_dir = out_dir / "per_set"
    combined_dir = out_dir / "combined"
    manifest_out = out_dir / "functional_gene_list_manifest.tsv"
    summary_out = out_dir / "functional_gene_list_summary.tsv"
    recurrent_out = combined_dir / f"recurrent_genes_min{args.min_recurrence}.tsv"
    union_symbols_out = combined_dir / "selected_sets_union.gene_symbols.txt"
    union_ids_out = combined_dir / "selected_sets_union.gene_ids.txt"
    log_path = root / "logs" / "prepare_functional_gene_lists.log"

    top_rows = read_tsv_rows(top_by_set_path)
    overview_rows = read_tsv_rows(candidate_overview_path)
    recurrence_rows = read_tsv_rows(candidate_recurrence_path)
    set_summary_rows_input = read_tsv_rows(candidate_set_summary_path)

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Top-by-set TSV: {top_by_set_path}", log_path)
    log(f"[INFO] Candidate overview TSV: {candidate_overview_path}", log_path)
    log(f"[INFO] Candidate recurrence TSV: {candidate_recurrence_path}", log_path)
    log(f"[INFO] Candidate set summary TSV: {candidate_set_summary_path}", log_path)

    if args.regions:
        allowed_regions = {normalize_text(x) for x in args.regions if normalize_text(x)}
        top_rows = [r for r in top_rows if normalize_text(r.get("region_class")) in allowed_regions]
        overview_rows = [r for r in overview_rows if normalize_text(r.get("region_class")) in allowed_regions]
        recurrence_rows = [
            r for r in recurrence_rows
            if any(region in allowed_regions for region in normalize_text(r.get("region_classes")).split(";") if region)
        ]
        set_summary_rows_input = [r for r in set_summary_rows_input if normalize_text(r.get("region_class")) in allowed_regions]

    if args.set_names:
        allowed_sets = {normalize_text(x) for x in args.set_names if normalize_text(x)}
        top_rows = [r for r in top_rows if normalize_text(r.get("set_name")) in allowed_sets]
        overview_rows = [r for r in overview_rows if normalize_text(r.get("set_name")) in allowed_sets]
        recurrence_rows = [
            r for r in recurrence_rows
            if any(s in allowed_sets for s in normalize_text(r.get("set_names")).split(";") if s)
        ]
        set_summary_rows_input = [r for r in set_summary_rows_input if normalize_text(r.get("set_name")) in allowed_sets]

    if args.comparison_ids:
        allowed_cmp = {normalize_text(x) for x in args.comparison_ids if normalize_text(x)}
        top_rows = [r for r in top_rows if normalize_text(r.get("comparison_id")) in allowed_cmp]
        overview_rows = [r for r in overview_rows if normalize_text(r.get("comparison_id")) in allowed_cmp]
        recurrence_rows = [
            r for r in recurrence_rows
            if any(c in allowed_cmp for c in normalize_text(r.get("comparison_ids")).split(";") if c)
        ]
        set_summary_rows_input = [r for r in set_summary_rows_input if normalize_text(r.get("comparison_id")) in allowed_cmp]

    # group top rows by set
    grouped_top: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in top_rows:
        row2 = dict(row)
        row2["gene_label"] = gene_label(row2.get("gene_name"), row2.get("gene_id"))
        grouped_top[build_key(row2)].append(row2)

    for key in grouped_top:
        grouped_top[key] = sorted(
            grouped_top[key],
            key=lambda r: safe_int(r.get("gene_rank"), 999999),
        )

    overview_idx = {build_key(r): r for r in overview_rows}
    set_summary_idx = {build_key(r): r for r in set_summary_rows_input}

    manifest_rows: List[Dict[str, Any]] = []
    summary_rows: List[Dict[str, Any]] = []

    union_gene_symbols: Set[str] = set()
    union_gene_ids: Set[str] = set()

    eligible_keys = sorted(grouped_top.keys(), key=lambda x: (x[0], x[1], x[2]))
    log(f"[INFO] Eligible peak sets: {len(eligible_keys)}", log_path)

    for i, key in enumerate(eligible_keys, start=1):
        comparison_id, set_name, region_class = key
        rows = grouped_top[key]

        rows = [
            r for r in rows
            if safe_int(r.get("peak_hit_count"), 0) >= args.min_peak_hits
        ][: args.top_n_per_set]

        rows_out: List[Dict[str, Any]] = []
        gene_symbols: List[str] = []
        gene_ids: List[str] = []

        for row in rows:
            label = normalize_text(row.get("gene_label"))
            gid = strip_ensembl_version(row.get("gene_id"))
            if label:
                gene_symbols.append(label)
                union_gene_symbols.add(label)
            if gid:
                gene_ids.append(gid)
                union_gene_ids.add(gid)

            rows_out.append({
                "comparison_id": comparison_id,
                "set_name": set_name,
                "region_class": region_class,
                "gene_rank": safe_int(row.get("gene_rank"), 0),
                "gene_id": normalize_text(row.get("gene_id")),
                "gene_name": normalize_text(row.get("gene_name")),
                "gene_label": label,
                "peak_hit_count": safe_int(row.get("peak_hit_count"), 0),
                "unique_peak_count": safe_int(row.get("unique_peak_count"), 0),
                "transcript_count": safe_int(row.get("transcript_count"), 0),
                "transcript_ids": normalize_text(row.get("transcript_ids")),
                "peak_hit_fraction_within_set": safe_float(row.get("peak_hit_fraction_within_set"), 0.0),
            })

        comparison_subdir = per_set_dir / safe_filename_part(comparison_id)
        top_genes_tsv = comparison_subdir / f"{set_name}.{region_class}.top_genes.tsv"
        gene_symbols_txt = comparison_subdir / f"{set_name}.{region_class}.gene_symbols.txt"
        gene_ids_txt = comparison_subdir / f"{set_name}.{region_class}.gene_ids.txt"

        write_tsv(top_genes_tsv, rows_out, PER_SET_COLUMNS)
        write_txt_list(gene_symbols_txt, gene_symbols)
        write_txt_list(gene_ids_txt, gene_ids)

        ov = overview_idx.get(key, {})
        top_gene_name = normalize_text(ov.get("top_gene_name"))
        top_gene_id = normalize_text(ov.get("top_gene_id"))
        top_gene_peak_hits = safe_int(ov.get("top_gene_peak_hits"), 0)

        manifest_rows.append({
            "comparison_id": comparison_id,
            "set_name": set_name,
            "region_class": region_class,
            "n_exported_genes": len(rows_out),
            "top_gene_name": top_gene_name,
            "top_gene_id": top_gene_id,
            "top_gene_peak_hits": top_gene_peak_hits,
            "top_genes_tsv": str(top_genes_tsv.relative_to(root)),
            "gene_symbols_txt": str(gene_symbols_txt.relative_to(root)),
            "gene_ids_txt": str(gene_ids_txt.relative_to(root)),
            "status": "ok",
            "note": f"exported_top_n={args.top_n_per_set};min_peak_hits={args.min_peak_hits}",
        })

        summary_rows.append({
            "comparison_id": comparison_id,
            "set_name": set_name,
            "region_class": region_class,
            "n_exported_genes": len(rows_out),
            "top_gene_name": top_gene_name,
            "top_gene_id": top_gene_id,
            "top_gene_peak_hits": top_gene_peak_hits,
            "top5_gene_labels": top5_labels(rows_out),
            "status": "ok",
            "note": f"exported_top_n={args.top_n_per_set};min_peak_hits={args.min_peak_hits}",
        })

        log(
            f"[INFO] [{i}/{len(eligible_keys)}] Exported {comparison_id} | {set_name} | {region_class} | n={len(rows_out)}",
            log_path,
        )

    # recurrent genes table
    filtered_recurrence_rows: List[Dict[str, Any]] = []
    for row in recurrence_rows:
        n_sets_present = safe_int(row.get("n_sets_present"), 0)
        if n_sets_present < args.min_recurrence:
            continue

        gname = normalize_text(row.get("gene_name"))
        gid = strip_ensembl_version(row.get("gene_id"))
        label = gene_label(gname, gid)

        filtered_recurrence_rows.append({
            "gene_id": normalize_text(row.get("gene_id")),
            "gene_name": gname,
            "gene_label": label,
            "n_sets_present": n_sets_present,
            "set_names": normalize_text(row.get("set_names")),
            "region_classes": normalize_text(row.get("region_classes")),
            "comparison_ids": normalize_text(row.get("comparison_ids")),
            "total_peak_hits": safe_int(row.get("total_peak_hits"), 0),
            "total_unique_peaks": safe_int(row.get("total_unique_peaks"), 0),
            "max_peak_hits_in_one_set": safe_int(row.get("max_peak_hits_in_one_set"), 0),
            "best_set_name": normalize_text(row.get("best_set_name")),
            "best_region_class": normalize_text(row.get("best_region_class")),
            "best_comparison_id": normalize_text(row.get("best_comparison_id")),
        })

    filtered_recurrence_rows = sorted(
        filtered_recurrence_rows,
        key=lambda r: (
            -safe_int(r.get("n_sets_present"), 0),
            -safe_int(r.get("total_peak_hits"), 0),
            normalize_text(r.get("gene_label")),
        ),
    )

    write_tsv(recurrent_out, filtered_recurrence_rows, RECURRENT_COLUMNS)
    write_txt_list(union_symbols_out, sorted(union_gene_symbols))
    write_txt_list(union_ids_out, sorted(union_gene_ids))

    manifest_rows = sorted(
        manifest_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_id")),
            normalize_text(r.get("set_name")),
            normalize_text(r.get("region_class")),
        ),
    )
    summary_rows = sorted(
        summary_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_id")),
            normalize_text(r.get("set_name")),
            normalize_text(r.get("region_class")),
        ),
    )

    write_tsv(manifest_out, manifest_rows, MANIFEST_COLUMNS)
    write_tsv(summary_out, summary_rows, SUMMARY_COLUMNS)

    log("", log_path)
    log("[DONE] Functional gene list preparation completed.", log_path)
    log(f"[DONE] Manifest TSV: {manifest_out}", log_path)
    log(f"[DONE] Summary TSV: {summary_out}", log_path)
    log(f"[DONE] Union gene symbols: {union_symbols_out}", log_path)
    log(f"[DONE] Union gene IDs: {union_ids_out}", log_path)
    log(f"[DONE] Recurrent genes TSV: {recurrent_out}", log_path)

    top_rows_for_log = sorted(
        summary_rows,
        key=lambda r: (
            -safe_int(r.get("top_gene_peak_hits"), 0),
            normalize_text(r.get("comparison_id")),
        ),
    )[:10]

    log("[INFO] Top exported functional sets:", log_path)
    for row in top_rows_for_log:
        log(
            "  " + " | ".join([
                normalize_text(row.get("comparison_id")),
                normalize_text(row.get("set_name")),
                normalize_text(row.get("region_class")),
                f"top_gene={normalize_text(row.get('top_gene_name')) or normalize_text(row.get('top_gene_id'))}",
                f"top_hits={normalize_text(row.get('top_gene_peak_hits'))}",
                f"top5={normalize_text(row.get('top5_gene_labels'))}",
            ]),
            log_path,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())