#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build download manifest from ENCODE experiment/file metadata.

Inputs:
- metadata/encode/experiments.tsv
- metadata/encode/files.tsv

Output:
- metadata/manifest/download_manifest.tsv

Default conservative filters:
- experiments.qc_keep == yes
- files.output_type == peaks
- files.file_format == bed
- files.assembly == GRCh38

Recommended usage:
    cd ~/wangshuo/rna_dynamics
    python scripts/01_metadata/02_build_download_manifest.py

Optional examples:
    python scripts/01_metadata/02_build_download_manifest.py --dry-run
    python scripts/01_metadata/02_build_download_manifest.py --include-manual-check
    python scripts/01_metadata/02_build_download_manifest.py --assembly hg19
    python scripts/01_metadata/02_build_download_manifest.py --output-types peaks "conservative IDR thresholded peaks"
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set


DEFAULT_ASSEMBLY = "GRCh38"
DEFAULT_OUTPUT_TYPES = ["peaks"]
DEFAULT_FILE_FORMATS = ["bed"]
DEFAULT_QC_KEEP = ["yes"]


MANIFEST_COLUMNS = [
    "file_accession",
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "file_format",
    "output_type",
    "assembly",
    "download_url",
    "local_path",
    "download_flag",
]


def log(msg: str) -> None:
    print(msg, flush=True)


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


def normalize_lower(x: Any) -> str:
    return normalize_text(x).lower()


def safe_filename_part(text: str) -> str:
    """
    Turn arbitrary text into a filesystem-safe path component.
    """
    s = normalize_text(text)
    if not s:
        return "NA"
    allowed = []
    for ch in s:
        if ch.isalnum() or ch in {"-", "_", "."}:
            allowed.append(ch)
        elif ch in {" ", "/", "\\", ":", ";", ","}:
            allowed.append("_")
        else:
            allowed.append("_")
    out = "".join(allowed)
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_") or "NA"


def split_semicolon_values(text: str) -> List[str]:
    s = normalize_text(text)
    if not s:
        return []
    return [x.strip() for x in s.split(";") if x.strip()]


def choose_primary_biosample(biosample_text: str) -> str:
    """
    experiments.tsv may contain one or multiple biosamples joined by ';'
    For current seed-scale work, keep the first one as primary display value.
    """
    vals = split_semicolon_values(biosample_text)
    if not vals:
        return ""
    return vals[0]


def build_experiment_index(experiment_rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    idx: Dict[str, Dict[str, str]] = {}
    for row in experiment_rows:
        exp_acc = normalize_text(row.get("experiment_accession"))
        if exp_acc:
            idx[exp_acc] = row
    return idx


def make_local_path(
    target_label: str,
    biosample_term_name: str,
    experiment_accession: str,
    file_accession: str,
    output_type: str,
    file_format: str,
    assembly: str,
) -> str:
    target = safe_filename_part(target_label)
    biosample = safe_filename_part(biosample_term_name)
    exp = safe_filename_part(experiment_accession)
    file_acc = safe_filename_part(file_accession)
    output = safe_filename_part(output_type)
    fmt = safe_filename_part(file_format)
    asm = safe_filename_part(assembly)

    return f"raw/encode/peaks/{target}/{biosample}/{exp}/{file_acc}.{output}.{asm}.{fmt}"


def row_passes_qc(
    exp_row: Dict[str, str],
    accepted_qc_keep: Set[str],
) -> bool:
    qc_keep = normalize_lower(exp_row.get("qc_keep"))
    return qc_keep in accepted_qc_keep


def row_passes_output_type(
    file_row: Dict[str, str],
    accepted_output_types: Set[str],
) -> bool:
    output_type = normalize_lower(file_row.get("output_type"))
    return output_type in accepted_output_types


def row_passes_file_format(
    file_row: Dict[str, str],
    accepted_file_formats: Set[str],
) -> bool:
    file_format = normalize_lower(file_row.get("file_format"))
    return file_format in accepted_file_formats


def row_passes_assembly(
    file_row: Dict[str, str],
    accepted_assemblies: Set[str],
) -> bool:
    assembly = normalize_lower(file_row.get("assembly"))
    return assembly in accepted_assemblies


def deduplicate_rows(rows: List[Dict[str, Any]], key_fields: List[str]) -> List[Dict[str, Any]]:
    seen = set()
    out = []
    for row in rows:
        key = tuple(normalize_text(row.get(k)) for k in key_fields)
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build ENCODE download manifest for peak files."
    )
    parser.add_argument(
        "--assembly",
        nargs="*",
        default=[DEFAULT_ASSEMBLY],
        help=f"Accepted assemblies (default: {DEFAULT_ASSEMBLY})",
    )
    parser.add_argument(
        "--output-types",
        nargs="*",
        default=DEFAULT_OUTPUT_TYPES,
        help=f"Accepted output types (default: {DEFAULT_OUTPUT_TYPES})",
    )
    parser.add_argument(
        "--file-formats",
        nargs="*",
        default=DEFAULT_FILE_FORMATS,
        help=f"Accepted file formats (default: {DEFAULT_FILE_FORMATS})",
    )
    parser.add_argument(
        "--include-manual-check",
        action="store_true",
        help="Also include experiments where qc_keep == manual_check",
    )
    parser.add_argument(
        "--include-no",
        action="store_true",
        help="Also include experiments where qc_keep == no",
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
        "--dry-run",
        action="store_true",
        help="Print summary only; do not write output file",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    exp_path = root / "metadata" / "encode" / "experiments.tsv"
    file_path = root / "metadata" / "encode" / "files.tsv"
    out_path = root / "metadata" / "manifest" / "download_manifest.tsv"

    experiment_rows = read_tsv_rows(exp_path)
    file_rows = read_tsv_rows(file_path)

    exp_index = build_experiment_index(experiment_rows)

    accepted_qc_keep = set(DEFAULT_QC_KEEP)
    if args.include_manual_check:
        accepted_qc_keep.add("manual_check")
    if args.include_no:
        accepted_qc_keep.add("no")

    accepted_output_types = {normalize_lower(x) for x in args.output_types if normalize_text(x)}
    accepted_file_formats = {normalize_lower(x) for x in args.file_formats if normalize_text(x)}
    accepted_assemblies = {normalize_lower(x) for x in args.assembly if normalize_text(x)}

    target_filter = {normalize_lower(x) for x in args.targets} if args.targets else None
    biosample_filter = {normalize_lower(x) for x in args.biosamples} if args.biosamples else None

    log(f"[INFO] Project root: {root}")
    log(f"[INFO] Input experiments: {len(experiment_rows)}")
    log(f"[INFO] Input files: {len(file_rows)}")
    log(f"[INFO] Accepted qc_keep: {sorted(accepted_qc_keep)}")
    log(f"[INFO] Accepted output_type: {sorted(accepted_output_types)}")
    log(f"[INFO] Accepted file_format: {sorted(accepted_file_formats)}")
    log(f"[INFO] Accepted assembly: {sorted(accepted_assemblies)}")
    if target_filter:
        log(f"[INFO] Target filter: {sorted(target_filter)}")
    if biosample_filter:
        log(f"[INFO] Biosample filter: {sorted(biosample_filter)}")

    total_seen = 0
    missing_experiment = 0
    filtered_qc = 0
    filtered_output = 0
    filtered_format = 0
    filtered_assembly = 0
    filtered_target = 0
    filtered_biosample = 0

    manifest_rows: List[Dict[str, Any]] = []

    for file_row in file_rows:
        total_seen += 1
        exp_acc = normalize_text(file_row.get("experiment_accession"))
        if not exp_acc or exp_acc not in exp_index:
            missing_experiment += 1
            continue

        exp_row = exp_index[exp_acc]

        if not row_passes_qc(exp_row, accepted_qc_keep):
            filtered_qc += 1
            continue

        if not row_passes_output_type(file_row, accepted_output_types):
            filtered_output += 1
            continue

        if not row_passes_file_format(file_row, accepted_file_formats):
            filtered_format += 1
            continue

        if not row_passes_assembly(file_row, accepted_assemblies):
            filtered_assembly += 1
            continue

        target_label = normalize_text(exp_row.get("target_label"))
        biosample_term_name = choose_primary_biosample(exp_row.get("biosample_term_name", ""))

        if target_filter and normalize_lower(target_label) not in target_filter:
            filtered_target += 1
            continue

        if biosample_filter and normalize_lower(biosample_term_name) not in biosample_filter:
            filtered_biosample += 1
            continue

        file_accession = normalize_text(file_row.get("file_accession"))
        output_type = normalize_text(file_row.get("output_type"))
        file_format = normalize_text(file_row.get("file_format"))
        assembly = normalize_text(file_row.get("assembly"))
        download_url = normalize_text(file_row.get("download_url"))

        manifest_rows.append(
            {
                "file_accession": file_accession,
                "experiment_accession": exp_acc,
                "target_label": target_label,
                "biosample_term_name": biosample_term_name,
                "file_format": file_format,
                "output_type": output_type,
                "assembly": assembly,
                "download_url": download_url,
                "local_path": make_local_path(
                    target_label=target_label,
                    biosample_term_name=biosample_term_name,
                    experiment_accession=exp_acc,
                    file_accession=file_accession,
                    output_type=output_type,
                    file_format=file_format,
                    assembly=assembly,
                ),
                "download_flag": "yes",
            }
        )

    manifest_rows = deduplicate_rows(manifest_rows, ["file_accession"])
    manifest_rows = sorted(
        manifest_rows,
        key=lambda x: (
            normalize_text(x.get("target_label")),
            normalize_text(x.get("biosample_term_name")),
            normalize_text(x.get("experiment_accession")),
            normalize_text(x.get("file_accession")),
        ),
    )

    kept_experiments = sorted({row["experiment_accession"] for row in manifest_rows})
    kept_targets = sorted({row["target_label"] for row in manifest_rows})
    kept_biosamples = sorted({row["biosample_term_name"] for row in manifest_rows})

    log("")
    log("[INFO] Filtering summary:")
    log(f"  Total file rows scanned      : {total_seen}")
    log(f"  Missing experiment rows      : {missing_experiment}")
    log(f"  Filtered by qc_keep          : {filtered_qc}")
    log(f"  Filtered by output_type      : {filtered_output}")
    log(f"  Filtered by file_format      : {filtered_format}")
    log(f"  Filtered by assembly         : {filtered_assembly}")
    log(f"  Filtered by target           : {filtered_target}")
    log(f"  Filtered by biosample        : {filtered_biosample}")
    log(f"  Manifest rows kept           : {len(manifest_rows)}")
    log(f"  Unique experiments kept      : {len(kept_experiments)}")
    log(f"  Unique targets kept          : {len(kept_targets)} -> {', '.join(kept_targets) if kept_targets else 'NA'}")
    log(f"  Unique biosamples kept       : {len(kept_biosamples)} -> {', '.join(kept_biosamples) if kept_biosamples else 'NA'}")

    if manifest_rows:
        log("")
        log("[INFO] First 5 manifest rows preview:")
        for row in manifest_rows[:5]:
            log(
                "  "
                + " | ".join(
                    [
                        row["file_accession"],
                        row["experiment_accession"],
                        row["target_label"],
                        row["biosample_term_name"],
                        row["output_type"],
                        row["assembly"],
                    ]
                )
            )

    if args.dry_run:
        log("")
        log("[DONE] Dry run only. Output file was not written.")
        return 0

    write_tsv(out_path, manifest_rows, MANIFEST_COLUMNS)

    log("")
    log(f"[DONE] Manifest written: {out_path}")
    log("[INFO] Recommended next step:")
    log("  python scripts/02_download/download_encode_data.py")

    return 0


if __name__ == "__main__":
    sys.exit(main())