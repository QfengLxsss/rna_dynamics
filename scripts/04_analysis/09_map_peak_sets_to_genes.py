#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Map region-specific peak subsets to gene-level tables.

Inputs
------
- processed/peaks/comparison_sets/by_region/region_specific_file_manifest.tsv
- processed/peaks/comparison_sets/by_region/region_specific_peak_summary.tsv

Per-subset required file
------------------------
- output_annotation_tsv_path

Outputs
-------
- results/tables/gene_mapping/peak_set_gene_mapping_manifest.tsv
- results/tables/gene_mapping/peak_set_gene_mapping_summary.tsv
- results/tables/gene_mapping/per_set/<comparison_id>/<set_name>.<region>.peak_to_gene.tsv
- results/tables/gene_mapping/per_set/<comparison_id>/<set_name>.<region>.gene_summary.tsv
- logs/map_peak_sets_to_genes.log

Notes
-----
- A peak may map to multiple genes.
- This script is explicitly gene-centric.
- transcript_ids are kept as auxiliary information and are NOT used to expand row count.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


PEAK_TO_GENE_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "peak_id",
    "chrom",
    "start",
    "end",
    "primary_region_class",
    "gene_id",
    "gene_name",
    "transcript_ids",
]

GENE_SUMMARY_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "gene_id",
    "gene_name",
    "peak_hit_count",
    "unique_peak_count",
    "transcript_count",
    "transcript_ids",
]

MANIFEST_COLUMNS = [
    "comparison_id",
    "set_name",
    "region_class",
    "input_annotation_tsv_path",
    "output_peak_to_gene_tsv",
    "output_gene_summary_tsv",
    "n_peak_rows",
    "n_peak_gene_links",
    "n_unique_genes",
    "status",
    "note",
]

SUMMARY_COLUMNS = [
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


def split_semicolon_field(text: str) -> List[str]:
    x = normalize_text(text)
    if not x:
        return []
    parts = [p.strip() for p in x.split(";")]
    return [p for p in parts if p]


def build_gene_pairs(gene_ids: List[str], gene_names: List[str]) -> List[Tuple[str, str]]:
    """
    Build (gene_id, gene_name) pairs conservatively.

    We only expand by gene dimension, not transcript dimension.
    """
    pairs: List[Tuple[str, str]] = []

    if gene_ids and gene_names:
        n = min(len(gene_ids), len(gene_names))
        for i in range(n):
            gid = normalize_text(gene_ids[i])
            gname = normalize_text(gene_names[i])
            if gid or gname:
                pairs.append((gid, gname))

        # if one side is longer, keep remaining info instead of losing it
        if len(gene_ids) > n:
            for gid in gene_ids[n:]:
                gid = normalize_text(gid)
                if gid:
                    pairs.append((gid, ""))
        if len(gene_names) > n:
            for gname in gene_names[n:]:
                gname = normalize_text(gname)
                if gname:
                    pairs.append(("", gname))

    elif gene_ids:
        for gid in gene_ids:
            gid = normalize_text(gid)
            if gid:
                pairs.append((gid, ""))

    elif gene_names:
        for gname in gene_names:
            gname = normalize_text(gname)
            if gname:
                pairs.append(("", gname))

    # de-duplicate while preserving order
    seen = set()
    uniq_pairs = []
    for pair in pairs:
        if pair not in seen:
            seen.add(pair)
            uniq_pairs.append(pair)
    return uniq_pairs


def parse_peak_annotation_to_gene_rows(
    rows: List[Dict[str, str]],
    comparison_id: str,
    set_name: str,
    region_class: str,
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []

    for row in rows:
        peak_id = normalize_text(row.get("peak_id"))
        chrom = normalize_text(row.get("chrom"))
        start = normalize_text(row.get("start"))
        end = normalize_text(row.get("end"))
        primary_region = normalize_text(row.get("primary_region_class"))

        gene_ids = split_semicolon_field(row.get("overlap_gene_ids"))
        gene_names = split_semicolon_field(row.get("overlap_gene_names"))
        transcript_ids = split_semicolon_field(row.get("overlap_transcript_ids"))
        transcript_ids_str = ";".join(sorted(set([t for t in transcript_ids if normalize_text(t)])))

        gene_pairs = build_gene_pairs(gene_ids, gene_names)

        # 如果一个 peak 没有 gene 信息，就跳过 gene mapping
        if not gene_pairs:
            continue

        for gene_id, gene_name in gene_pairs:
            out.append({
                "comparison_id": comparison_id,
                "set_name": set_name,
                "region_class": region_class,
                "peak_id": peak_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "primary_region_class": primary_region,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "transcript_ids": transcript_ids_str,
            })

    return out


def summarize_genes(peak_to_gene_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for row in peak_to_gene_rows:
        gene_id = normalize_text(row.get("gene_id"))
        gene_name = normalize_text(row.get("gene_name"))

        # 忽略完全空白的“gene”
        if not gene_id and not gene_name:
            continue

        key = (gene_id, gene_name)
        if key not in grouped:
            grouped[key] = {
                "comparison_id": normalize_text(row.get("comparison_id")),
                "set_name": normalize_text(row.get("set_name")),
                "region_class": normalize_text(row.get("region_class")),
                "gene_id": gene_id,
                "gene_name": gene_name,
                "peak_ids": set(),
                "transcript_ids": set(),
                "peak_hit_count": 0,
            }

        grouped[key]["peak_hit_count"] += 1

        peak_id = normalize_text(row.get("peak_id"))
        transcript_ids = split_semicolon_field(row.get("transcript_ids"))

        if peak_id:
            grouped[key]["peak_ids"].add(peak_id)
        for tid in transcript_ids:
            if tid:
                grouped[key]["transcript_ids"].add(tid)

    out: List[Dict[str, Any]] = []
    for _, g in grouped.items():
        transcript_ids = sorted(g["transcript_ids"])
        out.append({
            "comparison_id": g["comparison_id"],
            "set_name": g["set_name"],
            "region_class": g["region_class"],
            "gene_id": g["gene_id"],
            "gene_name": g["gene_name"],
            "peak_hit_count": g["peak_hit_count"],
            "unique_peak_count": len(g["peak_ids"]),
            "transcript_count": len(transcript_ids),
            "transcript_ids": ";".join(transcript_ids),
        })

    out = sorted(
        out,
        key=lambda r: (
            -safe_int(r.get("peak_hit_count"), 0),
            -safe_int(r.get("unique_peak_count"), 0),
            normalize_text(r.get("gene_name")),
            normalize_text(r.get("gene_id")),
        ),
    )
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map region-specific peak subsets to gene-level tables."
    )
    parser.add_argument(
        "--region-file-manifest",
        default=None,
        help="Override processed/peaks/comparison_sets/by_region/region_specific_file_manifest.tsv",
    )
    parser.add_argument(
        "--region-summary",
        default=None,
        help="Override processed/peaks/comparison_sets/by_region/region_specific_peak_summary.tsv",
    )
    parser.add_argument(
        "--regions",
        nargs="*",
        default=None,
        help="Restrict to these regions, e.g. intron CDS",
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
        "--limit",
        type=int,
        default=None,
        help="Only process first N eligible rows",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    by_region_dir = root / "processed" / "peaks" / "comparison_sets" / "by_region"
    region_file_manifest = (
        Path(args.region_file_manifest).resolve()
        if args.region_file_manifest
        else by_region_dir / "region_specific_file_manifest.tsv"
    )
    region_summary = (
        Path(args.region_summary).resolve()
        if args.region_summary
        else by_region_dir / "region_specific_peak_summary.tsv"
    )

    out_dir = root / "results" / "tables" / "gene_mapping"
    per_set_dir = out_dir / "per_set"
    manifest_out = out_dir / "peak_set_gene_mapping_manifest.tsv"
    summary_out = out_dir / "peak_set_gene_mapping_summary.tsv"
    log_path = root / "logs" / "map_peak_sets_to_genes.log"

    file_rows = read_tsv_rows(region_file_manifest)
    _summary_rows_input = read_tsv_rows(region_summary)

    if args.regions:
        allowed_regions = {normalize_text(x) for x in args.regions if normalize_text(x)}
        file_rows = [r for r in file_rows if normalize_text(r.get("region_class")) in allowed_regions]

    if args.set_names:
        allowed_sets = {normalize_text(x) for x in args.set_names if normalize_text(x)}
        file_rows = [r for r in file_rows if normalize_text(r.get("set_name")) in allowed_sets]

    if args.comparison_ids:
        allowed_cmp = {normalize_text(x) for x in args.comparison_ids if normalize_text(x)}
        file_rows = [r for r in file_rows if normalize_text(r.get("comparison_id")) in allowed_cmp]

    file_rows = [r for r in file_rows if normalize_text(r.get("status")) == "ok"]

    file_rows = sorted(
        file_rows,
        key=lambda r: (
            normalize_text(r.get("comparison_id")),
            normalize_text(r.get("set_name")),
            normalize_text(r.get("region_class")),
        ),
    )

    if args.limit is not None:
        file_rows = file_rows[:args.limit]

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Region file manifest TSV: {region_file_manifest}", log_path)
    log(f"[INFO] Region summary TSV: {region_summary}", log_path)
    log(f"[INFO] Eligible region subsets: {len(file_rows)}", log_path)

    manifest_rows: List[Dict[str, Any]] = []
    summary_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(file_rows, start=1):
        comparison_id = normalize_text(row.get("comparison_id"))
        set_name = normalize_text(row.get("set_name"))
        region_class = normalize_text(row.get("region_class"))
        input_ann_rel = normalize_text(row.get("output_annotation_tsv_path"))
        input_ann = root / input_ann_rel

        log(
            f"[INFO] [{i}/{len(file_rows)}] Mapping genes for {comparison_id} | {set_name} | {region_class}",
            log_path,
        )

        if not input_ann.exists():
            manifest_rows.append({
                "comparison_id": comparison_id,
                "set_name": set_name,
                "region_class": region_class,
                "input_annotation_tsv_path": input_ann_rel,
                "output_peak_to_gene_tsv": "",
                "output_gene_summary_tsv": "",
                "n_peak_rows": 0,
                "n_peak_gene_links": 0,
                "n_unique_genes": 0,
                "status": "missing_input",
                "note": "input_annotation_tsv_missing",
            })
            continue

        ann_rows = read_tsv_rows(input_ann)
        peak_to_gene_rows = parse_peak_annotation_to_gene_rows(
            rows=ann_rows,
            comparison_id=comparison_id,
            set_name=set_name,
            region_class=region_class,
        )
        gene_summary_rows = summarize_genes(peak_to_gene_rows)

        comparison_subdir = per_set_dir / safe_filename_part(comparison_id)
        peak_to_gene_out = comparison_subdir / f"{set_name}.{region_class}.peak_to_gene.tsv"
        gene_summary_out = comparison_subdir / f"{set_name}.{region_class}.gene_summary.tsv"

        write_tsv(peak_to_gene_out, peak_to_gene_rows, PEAK_TO_GENE_COLUMNS)
        write_tsv(gene_summary_out, gene_summary_rows, GENE_SUMMARY_COLUMNS)

        n_peak_rows = len(ann_rows)
        n_peak_gene_links = len(peak_to_gene_rows)
        n_unique_genes = len(gene_summary_rows)

        top_gene_name = ""
        top_gene_id = ""
        top_gene_peak_hits = 0
        if gene_summary_rows:
            top_gene_name = normalize_text(gene_summary_rows[0].get("gene_name"))
            top_gene_id = normalize_text(gene_summary_rows[0].get("gene_id"))
            top_gene_peak_hits = safe_int(gene_summary_rows[0].get("peak_hit_count"), 0)

        manifest_rows.append({
            "comparison_id": comparison_id,
            "set_name": set_name,
            "region_class": region_class,
            "input_annotation_tsv_path": input_ann_rel,
            "output_peak_to_gene_tsv": str(peak_to_gene_out.relative_to(root)),
            "output_gene_summary_tsv": str(gene_summary_out.relative_to(root)),
            "n_peak_rows": n_peak_rows,
            "n_peak_gene_links": n_peak_gene_links,
            "n_unique_genes": n_unique_genes,
            "status": "ok",
            "note": "gene_mapping_completed",
        })

        summary_rows.append({
            "comparison_id": comparison_id,
            "set_name": set_name,
            "region_class": region_class,
            "n_peak_rows": n_peak_rows,
            "n_peak_gene_links": n_peak_gene_links,
            "n_unique_genes": n_unique_genes,
            "top_gene_name": top_gene_name,
            "top_gene_id": top_gene_id,
            "top_gene_peak_hits": top_gene_peak_hits,
            "status": "ok",
            "note": "gene_mapping_completed",
        })

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

    n_ok = sum(1 for r in manifest_rows if normalize_text(r.get("status")) == "ok")
    n_bad = len(manifest_rows) - n_ok

    log("", log_path)
    log("[DONE] Peak-set to gene mapping completed.", log_path)
    log(f"[DONE] Manifest TSV: {manifest_out}", log_path)
    log(f"[DONE] Summary TSV: {summary_out}", log_path)
    log(f"[INFO] OK rows      : {n_ok}", log_path)
    log(f"[INFO] Problem rows : {n_bad}", log_path)

    best_rows = sorted(
        [r for r in summary_rows if normalize_text(r.get("status")) == "ok"],
        key=lambda r: (
            -safe_int(r.get("n_unique_genes"), 0),
            normalize_text(r.get("comparison_id")),
        ),
    )[:12]

    log("[INFO] Top mapped subsets:", log_path)
    for r in best_rows:
        log(
            "  " + " | ".join([
                normalize_text(r.get("comparison_id")),
                normalize_text(r.get("set_name")),
                normalize_text(r.get("region_class")),
                f"genes={normalize_text(r.get('n_unique_genes'))}",
                f"top_gene={normalize_text(r.get('top_gene_name'))}",
                f"top_hits={normalize_text(r.get('top_gene_peak_hits'))}",
            ]),
            log_path,
        )

    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())