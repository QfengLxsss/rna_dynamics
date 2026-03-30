#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import sys


DIRS = [
    # docs
    "docs",
    # metadata
    "metadata",
    "metadata/rbp",
    "metadata/encode",
    "metadata/biosample",
    "metadata/manifest",
    # raw
    "raw",
    "raw/encode",
    "raw/encode/peaks",
    "raw/encode/bam",
    "raw/encode/fastq",
    "raw/encode/bigwig",
    "raw/external",
    # processed
    "processed",
    "processed/peaks",
    "processed/peaks/merged",
    "processed/peaks/filtered",
    "processed/peaks/annotated",
    "processed/matrices",
    "processed/features",
    "processed/features/motif",
    "processed/features/region_annotation",
    # results
    "results",
    "results/figures",
    "results/tables",
    "results/enrichment",
    "results/comparisons",
    # scripts
    "scripts",
    "scripts/01_metadata",
    "scripts/02_download",
    "scripts/03_processing",
    "scripts/04_analysis",
    "scripts/utils",
    # notebooks
    "notebooks",
    "notebooks/exploratory",
    "notebooks/figures",
    # logs
    "logs",
]

FILES = {
    # docs
    "docs/project_scope.md": """# Project Scope

## Project
RBP binding site collection and analysis across multiple biosamples/tissues, starting from ENCODE eCLIP data.

## Initial strategy
1. Start with a small seed set of RBP targets
2. Focus on high-quality ENCODE eCLIP experiments
3. Build metadata tables first
4. Download only necessary files
5. Standardize peaks before downstream analysis
""",
    "docs/data_standard.md": """# Data Standard

## Core metadata tables
- metadata/rbp/rbp_seed_list.tsv
- metadata/encode/experiments.tsv
- metadata/encode/files.tsv
- metadata/biosample/tissue_mapping.tsv
- metadata/manifest/download_manifest.tsv

## Core file formats
- Peaks: BED/BED-like
- Signal tracks: bigWig
- Alignments: BAM
""",
    "docs/encode_notes.md": """# ENCODE Notes

Write ENCODE-specific filtering rules, QC criteria, and API field notes here.
""",
    "docs/rbp_background.md": """# RBP Background

Write background notes for seed RBPs here.
""",
    # metadata
    "metadata/rbp/rbp_seed_list.tsv": "rbp_symbol\tpriority\treason\n",
    "metadata/rbp/rbp_expanded_list.tsv": "rbp_symbol\tpriority\tsource\tnote\n",
    "metadata/encode/experiments.tsv": (
        "experiment_accession\ttarget_label\tassay_title\tbiosample_term_name\t"
        "biosample_type\torganism\tlab\tstatus\taudit_warning\tqc_keep\tqc_reason\n"
    ),
    "metadata/encode/files.tsv": (
        "file_accession\texperiment_accession\tfile_format\toutput_type\tassembly\t"
        "biological_replicates\ttechnical_replicates\tdownload_url\n"
    ),
    "metadata/encode/audit.tsv": (
        "experiment_accession\tfile_accession\taudit_level\taudit_category\tmessage\n"
    ),
    "metadata/biosample/biosample_list.tsv": (
        "biosample_term_name\tbiosample_type\torganism\tsource\n"
    ),
    "metadata/biosample/tissue_mapping.tsv": (
        "biosample_term_name\tcanonical_tissue\tsystem\tmapping_note\n"
    ),
    "metadata/manifest/download_manifest.tsv": (
        "file_accession\texperiment_accession\ttarget_label\tbiosample_term_name\t"
        "file_format\toutput_type\tdownload_url\tlocal_path\tdownload_flag\n"
    ),
    # results placeholders
    "results/tables/.gitkeep": "",
    "results/figures/.gitkeep": "",
    "results/enrichment/.gitkeep": "",
    "results/comparisons/.gitkeep": "",
    # notebooks placeholders
    "notebooks/exploratory/.gitkeep": "",
    "notebooks/figures/.gitkeep": "",
    # raw / processed placeholders
    "raw/encode/peaks/.gitkeep": "",
    "raw/encode/bam/.gitkeep": "",
    "raw/encode/fastq/.gitkeep": "",
    "raw/encode/bigwig/.gitkeep": "",
    "raw/external/.gitkeep": "",
    "processed/peaks/merged/.gitkeep": "",
    "processed/peaks/filtered/.gitkeep": "",
    "processed/peaks/annotated/.gitkeep": "",
    "processed/matrices/.gitkeep": "",
    "processed/features/motif/.gitkeep": "",
    "processed/features/region_annotation/.gitkeep": "",
    "logs/.gitkeep": "",
}

ROOT_SENTINELS = ["README.md", "docs", "metadata", "raw", "processed", "scripts"]


def find_project_root() -> Path:
    cwd = Path.cwd().resolve()
    if all((cwd / s).exists() for s in ROOT_SENTINELS[:1]):
        return cwd

    current = cwd
    for _ in range(6):
        if any((current / s).exists() for s in ROOT_SENTINELS):
            return current
        if current.parent == current:
            break
        current = current.parent

    return cwd


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def ensure_file(path: Path, content: str) -> str:
    if path.exists():
        return "exists"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return "created"


def main() -> int:
    root = find_project_root()

    print(f"[INFO] Project root: {root}")

    created_dirs = 0
    for rel in DIRS:
        p = root / rel
        if not p.exists():
            p.mkdir(parents=True, exist_ok=True)
            created_dirs += 1

    created_files = 0
    existing_files = 0
    for rel, content in FILES.items():
        status = ensure_file(root / rel, content)
        if status == "created":
            created_files += 1
        else:
            existing_files += 1

    print(f"[DONE] Directories ensured: {len(DIRS)}")
    print(f"[DONE] Newly created directories: {created_dirs}")
    print(f"[DONE] Newly created files: {created_files}")
    print(f"[DONE] Existing files kept unchanged: {existing_files}")

    print("\n[INFO] Recommended next files to edit first:")
    print("  1. metadata/rbp/rbp_seed_list.tsv")
    print("  2. metadata/biosample/tissue_mapping.tsv")
    print("  3. docs/encode_notes.md")
    print("  4. scripts/01_metadata/01_collect_encode_eclip_metadata.py")

    return 0


if __name__ == "__main__":
    sys.exit(main())