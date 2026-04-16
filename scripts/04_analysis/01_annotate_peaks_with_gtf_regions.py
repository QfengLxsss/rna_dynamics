#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Annotate prepared peak BED6 files with GTF-derived genomic regions.

Required input
--------------
- processed/peaks/annotation_inputs/annotation_input_manifest.tsv
- processed/peaks/annotation_inputs/bed6/<...>/*.annotation_ready.bed
- --gtf /raw/annotations/gencode.v49.basic.annotation.gtf.gz

Outputs
-------
- processed/peaks/annotations/per_file/<target>/<biosample>/<experiment>/<file>.peak_region_annotations.tsv
- processed/peaks/annotations/annotation_summary.tsv
- processed/peaks/annotations/annotation_file_manifest.tsv
- logs/annotate_peaks_with_gtf_regions.log

Primary region priority
-----------------------
five_prime_utr
three_prime_utr
utr
CDS
exon
intron
gene
intergenic

Notes
-----
- GTF is converted to BED-like 0-based half-open intervals internally.
- Introns are derived from exon coordinates grouped by transcript_id.
- If chromosome naming between peaks and GTF is incompatible, the script exits by default.
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


REGION_PRIORITY = [
    "five_prime_utr",
    "three_prime_utr",
    "utr",
    "CDS",
    "exon",
    "intron",
    "gene",
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

FILE_MANIFEST_COLUMNS = [
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "selected_file_accession",
    "annotation_input_bed6_path",
    "output_annotation_tsv_path",
    "n_peaks",
    "gtf_path",
    "gtf_chrom_style",
    "peak_chrom_style",
    "chrom_compatible",
    "status",
    "note",
]

SUMMARY_COLUMNS = [
    "target_label",
    "biosample_term_name",
    "experiment_accession",
    "selected_file_accession",
    "n_peaks",
    "n_five_prime_utr",
    "n_three_prime_utr",
    "n_utr",
    "n_CDS",
    "n_exon",
    "n_intron",
    "n_gene",
    "n_intergenic",
    "annotation_rate_non_intergenic",
    "status",
    "note",
]

BIN_SIZE = 100_000


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


def parse_gtf_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, value = part.split(" ", 1)
        value = value.strip().strip('"')
        attrs[key] = value
    return attrs


def canonical_gtf_feature(feature: str) -> Optional[str]:
    x = feature.strip()
    lower = x.lower()
    if lower == "gene":
        return "gene"
    if lower == "exon":
        return "exon"
    if lower == "cds":
        return "CDS"
    if lower in {"five_prime_utr", "five_prime_utr_region", "5utr", "five_prime_utr"}:
        return "five_prime_utr"
    if lower in {"three_prime_utr", "three_prime_utr_region", "3utr", "three_prime_utr"}:
        return "three_prime_utr"
    if lower == "utr":
        return "utr"
    if lower == "transcript":
        return "transcript"
    return None


def add_interval_to_index(
    region_records: Dict[str, Dict[str, List[Tuple[int, int, str, str, str, str]]]],
    region_bins: Dict[str, Dict[str, Dict[int, List[int]]]],
    region_class: str,
    chrom: str,
    record: Tuple[int, int, str, str, str, str],
) -> None:
    if chrom not in region_records[region_class]:
        region_records[region_class][chrom] = []
        region_bins[region_class][chrom] = defaultdict(list)

    idx = len(region_records[region_class][chrom])
    region_records[region_class][chrom].append(record)

    start, end = record[0], record[1]
    start_bin = start // BIN_SIZE
    end_bin = (end - 1) // BIN_SIZE if end > 0 else start_bin
    for b in range(start_bin, end_bin + 1):
        region_bins[region_class][chrom][b].append(idx)


def inspect_first_peak_chrom_style(bed_path: Path, max_lines: int = 1000) -> Tuple[str, str]:
    chroms: List[str] = []
    with open_maybe_gzip(bed_path) as f:
        for i, raw in enumerate(f, start=1):
            if i > max_lines:
                break
            line = raw.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            chroms.append(fields[0])
            if len(chroms) >= 200:
                break
    return infer_chrom_style(chroms)


def make_output_annotation_path(
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
        / "annotations"
        / "per_file"
        / safe_filename_part(target_label)
        / safe_filename_part(biosample_term_name)
        / safe_filename_part(experiment_accession)
        / f"{safe_filename_part(file_accession)}.peak_region_annotations.tsv"
    )


def parse_gtf_and_build_indices(
    gtf_path: Path,
    log_path: Optional[Path] = None,
) -> Tuple[
    Dict[str, Dict[str, List[Tuple[int, int, str, str, str, str]]]],
    Dict[str, Dict[str, Dict[int, List[int]]]],
    Tuple[str, str],
]:
    log(f"[INFO] Parsing GTF: {gtf_path}", log_path)

    region_records: Dict[str, Dict[str, List[Tuple[int, int, str, str, str, str]]]] = defaultdict(dict)
    region_bins: Dict[str, Dict[str, Dict[int, List[int]]]] = defaultdict(dict)

    # transcript_id -> transcript info
    transcripts: Dict[str, Dict[str, Any]] = {}

    # fallback gene spans if gene features are absent
    gene_spans: Dict[Tuple[str, str, str], List[int]] = {}

    gtf_chroms: List[str] = []

    n_lines = 0
    n_features_kept = 0

    with open_maybe_gzip(gtf_path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            n_lines += 1

            chrom, source, feature, start_s, end_s, score, strand, frame, attrs_s = fields
            feature_class = canonical_gtf_feature(feature)
            if feature_class is None:
                continue

            start_1 = safe_int(start_s, None)
            end_1 = safe_int(end_s, None)
            if start_1 is None or end_1 is None:
                continue
            if end_1 < start_1:
                continue

            # GTF is 1-based inclusive; convert to 0-based half-open BED-like
            start_0 = start_1 - 1
            end_0 = end_1

            attrs = parse_gtf_attributes(attrs_s)
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id)
            transcript_id = attrs.get("transcript_id", "")

            gtf_chroms.append(chrom)
            record = (start_0, end_0, gene_id, gene_name, transcript_id, strand)

            if feature_class in {"gene", "exon", "CDS", "five_prime_utr", "three_prime_utr", "utr"}:
                add_interval_to_index(region_records, region_bins, feature_class, chrom, record)
                n_features_kept += 1

            if transcript_id:
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = {
                        "chrom": chrom,
                        "strand": strand,
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "start": start_0,
                        "end": end_0,
                        "exons": [],
                    }
                else:
                    transcripts[transcript_id]["start"] = min(transcripts[transcript_id]["start"], start_0)
                    transcripts[transcript_id]["end"] = max(transcripts[transcript_id]["end"], end_0)

                if feature_class == "exon":
                    transcripts[transcript_id]["exons"].append((start_0, end_0))

            if gene_id:
                gkey = (chrom, gene_id, gene_name)
                if gkey not in gene_spans:
                    gene_spans[gkey] = [start_0, end_0]
                else:
                    gene_spans[gkey][0] = min(gene_spans[gkey][0], start_0)
                    gene_spans[gkey][1] = max(gene_spans[gkey][1], end_0)

    # If gene features are absent for some genes, add gene spans from aggregated transcript/exon data
    existing_gene_keys = set()
    for chrom, recs in region_records["gene"].items():
        for start_0, end_0, gene_id, gene_name, transcript_id, strand in recs:
            existing_gene_keys.add((chrom, gene_id, gene_name))

    for (chrom, gene_id, gene_name), (start_0, end_0) in gene_spans.items():
        key = (chrom, gene_id, gene_name)
        if key in existing_gene_keys:
            continue
        record = (start_0, end_0, gene_id, gene_name, "", ".")
        add_interval_to_index(region_records, region_bins, "gene", chrom, record)
        n_features_kept += 1

    # Derive introns from transcript exon structures
    intron_count = 0
    for transcript_id, tinfo in transcripts.items():
        chrom = tinfo["chrom"]
        strand = tinfo["strand"]
        gene_id = tinfo["gene_id"]
        gene_name = tinfo["gene_name"]
        exons = sorted(set(tinfo["exons"]))
        if len(exons) < 2:
            continue
        for (s1, e1), (s2, e2) in zip(exons[:-1], exons[1:]):
            intron_start = e1
            intron_end = s2
            if intron_end <= intron_start:
                continue
            record = (intron_start, intron_end, gene_id, gene_name, transcript_id, strand)
            add_interval_to_index(region_records, region_bins, "intron", chrom, record)
            intron_count += 1

    gtf_style = infer_chrom_style(gtf_chroms)
    log(f"[INFO] GTF lines parsed: {n_lines}", log_path)
    log(f"[INFO] GTF feature intervals kept: {n_features_kept}", log_path)
    log(f"[INFO] Derived introns: {intron_count}", log_path)
    log(f"[INFO] GTF chrom style: {gtf_style[0]} | has_chr_prefix={gtf_style[1]}", log_path)

    return region_records, region_bins, gtf_style


def find_overlaps_for_region(
    chrom: str,
    start: int,
    end: int,
    region_class: str,
    region_records: Dict[str, Dict[str, List[Tuple[int, int, str, str, str, str]]]],
    region_bins: Dict[str, Dict[str, Dict[int, List[int]]]],
) -> List[Tuple[int, int, str, str, str, str]]:
    if chrom not in region_records.get(region_class, {}):
        return []

    recs = region_records[region_class][chrom]
    bins = region_bins[region_class][chrom]

    start_bin = start // BIN_SIZE
    end_bin = (end - 1) // BIN_SIZE if end > 0 else start_bin

    candidate_indices = set()
    for b in range(start_bin, end_bin + 1):
        if b in bins:
            candidate_indices.update(bins[b])

    overlaps: List[Tuple[int, int, str, str, str, str]] = []
    for idx in candidate_indices:
        feat_start, feat_end, gene_id, gene_name, transcript_id, strand = recs[idx]
        if feat_start < end and feat_end > start:
            overlaps.append((feat_start, feat_end, gene_id, gene_name, transcript_id, strand))
    return overlaps


def annotate_one_peak_file(
    input_bed6: Path,
    output_tsv: Path,
    region_records: Dict[str, Dict[str, List[Tuple[int, int, str, str, str, str]]]],
    region_bins: Dict[str, Dict[str, Dict[int, List[int]]]],
    target_label: str,
    biosample_term_name: str,
    experiment_accession: str,
    selected_file_accession: str,
) -> Dict[str, Any]:
    counts = Counter()
    n_peaks = 0

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with input_bed6.open("r", encoding="utf-8") as f, output_tsv.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=ANNOTATION_COLUMNS, delimiter="\t")
        writer.writeheader()

        for raw in f:
            line = raw.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 6:
                continue

            chrom, start_s, end_s, peak_id, score, strand = fields[:6]
            start = safe_int(start_s, None)
            end = safe_int(end_s, None)
            if start is None or end is None or end <= start:
                continue

            n_peaks += 1
            primary_region = "intergenic"
            overlap_gene_ids = ""
            overlap_gene_names = ""
            overlap_transcript_ids = ""
            n_primary_region_features = 0

            for region_class in REGION_PRIORITY:
                overlaps = find_overlaps_for_region(
                    chrom=chrom,
                    start=start,
                    end=end,
                    region_class=region_class,
                    region_records=region_records,
                    region_bins=region_bins,
                )
                if overlaps:
                    primary_region = region_class
                    gene_ids = sorted({x[2] for x in overlaps if normalize_text(x[2])})
                    gene_names = sorted({x[3] for x in overlaps if normalize_text(x[3])})
                    transcript_ids = sorted({x[4] for x in overlaps if normalize_text(x[4])})

                    overlap_gene_ids = ";".join(gene_ids)
                    overlap_gene_names = ";".join(gene_names)
                    overlap_transcript_ids = ";".join(transcript_ids)
                    n_primary_region_features = len(overlaps)
                    break

            counts[primary_region] += 1

            writer.writerow({
                "peak_id": peak_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "score": score,
                "strand": strand,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "experiment_accession": experiment_accession,
                "selected_file_accession": selected_file_accession,
                "primary_region_class": primary_region,
                "overlap_gene_ids": overlap_gene_ids,
                "overlap_gene_names": overlap_gene_names,
                "overlap_transcript_ids": overlap_transcript_ids,
                "n_primary_region_features": n_primary_region_features,
            })

    return {
        "n_peaks": n_peaks,
        "n_five_prime_utr": counts["five_prime_utr"],
        "n_three_prime_utr": counts["three_prime_utr"],
        "n_utr": counts["utr"],
        "n_CDS": counts["CDS"],
        "n_exon": counts["exon"],
        "n_intron": counts["intron"],
        "n_gene": counts["gene"],
        "n_intergenic": counts["intergenic"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate prepared peak BED6 files with GTF-derived genomic regions."
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Override processed/peaks/annotation_inputs/annotation_input_manifest.tsv",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="Path to GTF or GTF.GZ annotation file",
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
        help="Only process first N files",
    )
    parser.add_argument(
        "--allow-chrom-mismatch",
        action="store_true",
        help="Allow running even if peak/GTF chromosome naming styles mismatch",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    manifest_tsv = (
        Path(args.manifest).resolve()
        if args.manifest
        else root / "processed" / "peaks" / "annotation_inputs" / "annotation_input_manifest.tsv"
    )
    gtf_path = Path(args.gtf).resolve()
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    out_dir = root / "processed" / "peaks" / "annotations"
    file_manifest_out = out_dir / "annotation_file_manifest.tsv"
    summary_out = out_dir / "annotation_summary.tsv"
    log_path = root / "logs" / "annotate_peaks_with_gtf_regions.log"

    rows = read_tsv_rows(manifest_tsv)

    if args.targets:
        target_filter = {x.strip().lower() for x in args.targets if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("target_label")).lower() in target_filter]

    if args.biosamples:
        biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()}
        rows = [r for r in rows if normalize_text(r.get("biosample_term_name")).lower() in biosample_filter]

    if args.limit is not None:
        rows = rows[:args.limit]

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] Annotation input manifest: {manifest_tsv}", log_path)
    log(f"[INFO] GTF path: {gtf_path}", log_path)
    log(f"[INFO] Files to annotate: {len(rows)}", log_path)

    if not rows:
        log("[WARN] No files to annotate after filtering.", log_path)
        return 0

    region_records, region_bins, gtf_style = parse_gtf_and_build_indices(gtf_path, log_path)

    file_manifest_rows: List[Dict[str, Any]] = []
    summary_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(rows, start=1):
        experiment_accession = normalize_text(row.get("experiment_accession"))
        target_label = normalize_text(row.get("target_label"))
        biosample_term_name = normalize_text(row.get("biosample_term_name"))
        selected_file_accession = normalize_text(row.get("selected_file_accession"))
        input_bed_rel = normalize_text(row.get("annotation_ready_bed6_path"))
        input_bed = root / input_bed_rel

        peak_style = normalize_text(row.get("chrom_style"))
        peak_has_chr = normalize_text(row.get("has_chr_prefix"))
        gtf_has_chr = gtf_style[1]
        chrom_compatible = "yes" if peak_has_chr == gtf_has_chr else "no"

        log(
            f"[INFO] [{i}/{len(rows)}] Annotating {selected_file_accession} | "
            f"{target_label} | {biosample_term_name} | {experiment_accession}",
            log_path,
        )

        output_tsv = make_output_annotation_path(
            root=root,
            target_label=target_label,
            biosample_term_name=biosample_term_name,
            experiment_accession=experiment_accession,
            file_accession=selected_file_accession,
        )

        if not input_bed.exists():
            file_manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": selected_file_accession,
                "annotation_input_bed6_path": input_bed_rel,
                "output_annotation_tsv_path": str(output_tsv.relative_to(root)),
                "n_peaks": 0,
                "gtf_path": str(gtf_path),
                "gtf_chrom_style": gtf_style[0],
                "peak_chrom_style": peak_style,
                "chrom_compatible": chrom_compatible,
                "status": "missing_input_bed6",
                "note": "annotation_ready_bed6_not_found",
            })
            continue

        if chrom_compatible != "yes" and not args.allow_chrom_mismatch:
            file_manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": selected_file_accession,
                "annotation_input_bed6_path": input_bed_rel,
                "output_annotation_tsv_path": str(output_tsv.relative_to(root)),
                "n_peaks": 0,
                "gtf_path": str(gtf_path),
                "gtf_chrom_style": gtf_style[0],
                "peak_chrom_style": peak_style,
                "chrom_compatible": chrom_compatible,
                "status": "chrom_mismatch",
                "note": "peak_and_gtf_chrom_naming_incompatible_use_matching_gtf_or_allow_chrom_mismatch",
            })
            log(
                f"[ERROR] Chrom naming mismatch for {selected_file_accession}: "
                f"peak={peak_style}, gtf={gtf_style[0]}",
                log_path,
            )
            continue

        try:
            stats = annotate_one_peak_file(
                input_bed6=input_bed,
                output_tsv=output_tsv,
                region_records=region_records,
                region_bins=region_bins,
                target_label=target_label,
                biosample_term_name=biosample_term_name,
                experiment_accession=experiment_accession,
                selected_file_accession=selected_file_accession,
            )

            n_peaks = stats["n_peaks"]
            non_intergenic = (
                stats["n_five_prime_utr"]
                + stats["n_three_prime_utr"]
                + stats["n_utr"]
                + stats["n_CDS"]
                + stats["n_exon"]
                + stats["n_intron"]
                + stats["n_gene"]
            )
            annotation_rate = round(non_intergenic / n_peaks, 6) if n_peaks > 0 else 0.0

            file_manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": selected_file_accession,
                "annotation_input_bed6_path": input_bed_rel,
                "output_annotation_tsv_path": str(output_tsv.relative_to(root)),
                "n_peaks": n_peaks,
                "gtf_path": str(gtf_path),
                "gtf_chrom_style": gtf_style[0],
                "peak_chrom_style": peak_style,
                "chrom_compatible": chrom_compatible,
                "status": "ok",
                "note": "peak_region_annotation_completed",
            })

            summary_rows.append({
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "experiment_accession": experiment_accession,
                "selected_file_accession": selected_file_accession,
                "n_peaks": n_peaks,
                "n_five_prime_utr": stats["n_five_prime_utr"],
                "n_three_prime_utr": stats["n_three_prime_utr"],
                "n_utr": stats["n_utr"],
                "n_CDS": stats["n_CDS"],
                "n_exon": stats["n_exon"],
                "n_intron": stats["n_intron"],
                "n_gene": stats["n_gene"],
                "n_intergenic": stats["n_intergenic"],
                "annotation_rate_non_intergenic": annotation_rate,
                "status": "ok",
                "note": "peak_region_annotation_completed",
            })

        except Exception as e:
            file_manifest_rows.append({
                "experiment_accession": experiment_accession,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "selected_file_accession": selected_file_accession,
                "annotation_input_bed6_path": input_bed_rel,
                "output_annotation_tsv_path": str(output_tsv.relative_to(root)),
                "n_peaks": 0,
                "gtf_path": str(gtf_path),
                "gtf_chrom_style": gtf_style[0],
                "peak_chrom_style": peak_style,
                "chrom_compatible": chrom_compatible,
                "status": "error",
                "note": str(e),
            })
            log(f"[ERROR] Failed to annotate {selected_file_accession}: {e}", log_path)

    file_manifest_rows = sorted(
        file_manifest_rows,
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

    write_tsv(file_manifest_out, file_manifest_rows, FILE_MANIFEST_COLUMNS)
    write_tsv(summary_out, summary_rows, SUMMARY_COLUMNS)

    n_ok = sum(1 for r in file_manifest_rows if normalize_text(r.get("status")) == "ok")
    n_err = sum(1 for r in file_manifest_rows if normalize_text(r.get("status")) == "error")
    n_mismatch = sum(1 for r in file_manifest_rows if normalize_text(r.get("status")) == "chrom_mismatch")
    n_missing = sum(1 for r in file_manifest_rows if normalize_text(r.get("status")) == "missing_input_bed6")

    log("", log_path)
    log("[DONE] Peak region annotation completed.", log_path)
    log(f"[DONE] File manifest: {file_manifest_out}", log_path)
    log(f"[DONE] Summary TSV: {summary_out}", log_path)
    log(f"[INFO] OK files         : {n_ok}", log_path)
    log(f"[INFO] Error files      : {n_err}", log_path)
    log(f"[INFO] Chrom mismatches : {n_mismatch}", log_path)
    log(f"[INFO] Missing inputs   : {n_missing}", log_path)

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
                    f"CDS={normalize_text(r.get('n_CDS'))}",
                    f"intron={normalize_text(r.get('n_intron'))}",
                    f"intergenic={normalize_text(r.get('n_intergenic'))}",
                    f"annot_rate={normalize_text(r.get('annotation_rate_non_intergenic'))}",
                ]
            ),
            log_path,
        )

    return 0 if n_err == 0 and n_mismatch == 0 and n_missing == 0 else 1


if __name__ == "__main__":
    sys.exit(main())