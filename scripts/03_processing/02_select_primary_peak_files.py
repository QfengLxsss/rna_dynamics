#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Select one primary peak file per experiment from downloaded ENCODE peak files.

Upgraded selection strategy:
1. inventory status == ok
2. prefer larger peak_count (most important)
3. prefer output_type keywords that look like primary peak sets:
   idr / optimal / conservative / pseudoreplicated / overlap / pooled
4. prefer larger file size
5. tie-break by file_accession

Inputs:
- metadata/encode/files.tsv
- metadata/manifest/download_manifest.tsv
- processed/peaks/inventory/peak_file_inventory.tsv

Outputs:
- processed/peaks/selected/primary_peak_files.tsv
- processed/peaks/selected/peak_file_ranking.tsv
- logs/select_primary_peak_files.log
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


PRIMARY_COLUMNS = [
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "selected_file_accession",
    "selected_output_type",
    "selected_file_format",
    "selected_assembly",
    "selected_local_path",
    "selected_file_size_bytes",
    "selected_peak_count",
    "selection_score",
    "selection_reason",
]

RANKING_COLUMNS = [
    "experiment_accession",
    "target_label",
    "biosample_term_name",
    "file_accession",
    "output_type",
    "file_format",
    "assembly",
    "local_path",
    "file_size_bytes",
    "peak_count",
    "inventory_status",
    "ranking_score",
    "ranking_reason",
    "is_selected",
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


def safe_int(x: Any, default: int = 0) -> int:
    try:
        return int(str(x).strip())
    except Exception:
        return default


def build_index(rows: List[Dict[str, str]], key: str) -> Dict[str, Dict[str, str]]:
    idx: Dict[str, Dict[str, str]] = {}
    for row in rows:
        k = normalize_text(row.get(key))
        if k:
            idx[k] = row
    return idx


def open_maybe_gzip(path: Path):
    with path.open("rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def count_peak_lines(path: Path) -> int:
    """
    Count non-empty, non-comment, non-header lines in a BED-like file.
    """
    n = 0
    with open_maybe_gzip(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            lower = line.lower()
            if line.startswith("#") or lower.startswith("track") or lower.startswith("browser"):
                continue
            n += 1
    return n


def keyword_score(output_type: str) -> Tuple[int, List[str]]:
    text = normalize_lower(output_type)
    score = 0
    reasons: List[str] = []

    strong_keywords = {
        "idr": 80,
        "optimal": 70,
        "conservative": 65,
        "pseudoreplicated": 60,
        "pseudoreplicate": 60,
        "overlap": 55,
        "pooled": 50,
    }

    moderate_keywords = {
        "replicated": 35,
        "filtered": 20,
        "thresholded": 18,
        "peaks": 10,
    }

    penalty_keywords = {
        "replicate 1": -12,
        "replicate 2": -12,
        "biological replicate 1": -12,
        "biological replicate 2": -12,
        "pr1": -8,
        "pr2": -8,
    }

    for kw, val in strong_keywords.items():
        if kw in text:
            score += val
            reasons.append(f"+{val}:{kw}")

    for kw, val in moderate_keywords.items():
        if kw in text:
            score += val
            reasons.append(f"+{val}:{kw}")

    for kw, val in penalty_keywords.items():
        if kw in text:
            score += val
            reasons.append(f"{val}:{kw}")

    return score, reasons


def inventory_score(inventory_status: str) -> Tuple[int, List[str]]:
    status = normalize_lower(inventory_status)
    reasons: List[str] = []

    if status == "ok":
        return 100, ["+100:inventory_ok"]
    if status == "manual_check":
        return 20, ["+20:inventory_manual_check"]
    if status == "inconsistent_columns":
        return -50, ["-50:inventory_inconsistent_columns"]
    if status == "invalid_bed_like":
        return -80, ["-80:inventory_invalid_bed_like"]
    if status == "empty":
        return -100, ["-100:inventory_empty"]
    if status == "empty_or_no_data":
        return -100, ["-100:inventory_empty_or_no_data"]
    if status == "missing":
        return -120, ["-120:inventory_missing"]
    if status == "read_error":
        return -120, ["-120:inventory_read_error"]

    return 0, ["+0:inventory_unknown"]


def peak_count_score(peak_count: int) -> Tuple[int, List[str]]:
    """
    Peak count is the main upgraded signal.
    Uses coarse bins to avoid overwhelming all other criteria.
    """
    if peak_count >= 30000:
        return 60, [f"+60:peak_count_ge_30000({peak_count})"]
    if peak_count >= 20000:
        return 50, [f"+50:peak_count_ge_20000({peak_count})"]
    if peak_count >= 10000:
        return 40, [f"+40:peak_count_ge_10000({peak_count})"]
    if peak_count >= 5000:
        return 30, [f"+30:peak_count_ge_5000({peak_count})"]
    if peak_count >= 1000:
        return 20, [f"+20:peak_count_ge_1000({peak_count})"]
    if peak_count > 0:
        return 10, [f"+10:peak_count_gt_0({peak_count})"]
    return -30, [f"-30:peak_count_zero({peak_count})"]


def size_score(file_size_bytes: int) -> Tuple[int, List[str]]:
    reasons: List[str] = []
    if file_size_bytes >= 2_000_000:
        return 15, ["+15:size_ge_2MB"]
    if file_size_bytes >= 1_000_000:
        return 10, ["+10:size_ge_1MB"]
    if file_size_bytes >= 300_000:
        return 5, ["+5:size_ge_300KB"]
    if file_size_bytes > 0:
        return 1, ["+1:size_gt_0"]
    return -20, ["-20:size_zero"]


def format_score(file_format: str, assembly: str) -> Tuple[int, List[str]]:
    score = 0
    reasons: List[str] = []

    if normalize_lower(file_format) == "bed":
        score += 5
        reasons.append("+5:file_format_bed")

    if normalize_lower(assembly) == "grch38":
        score += 5
        reasons.append("+5:assembly_grch38")

    return score, reasons


def combined_score(candidate: Dict[str, Any]) -> Tuple[int, str]:
    total = 0
    reasons: List[str] = []

    inv_s, inv_r = inventory_score(candidate.get("inventory_status", ""))
    total += inv_s
    reasons.extend(inv_r)

    pc_s, pc_r = peak_count_score(safe_int(candidate.get("peak_count", 0)))
    total += pc_s
    reasons.extend(pc_r)

    kw_s, kw_r = keyword_score(candidate.get("output_type", ""))
    total += kw_s
    reasons.extend(kw_r)

    sz_s, sz_r = size_score(safe_int(candidate.get("file_size_bytes", 0)))
    total += sz_s
    reasons.extend(sz_r)

    fmt_s, fmt_r = format_score(
        candidate.get("file_format", ""),
        candidate.get("assembly", ""),
    )
    total += fmt_s
    reasons.extend(fmt_r)

    return total, ",".join(reasons)


def choose_best_candidate(rows: List[Dict[str, Any]]) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    scored_rows: List[Dict[str, Any]] = []

    for row in rows:
        score, reason = combined_score(row)
        x = dict(row)
        x["ranking_score"] = score
        x["ranking_reason"] = reason
        scored_rows.append(x)

    scored_rows = sorted(
        scored_rows,
        key=lambda r: (
            -safe_int(r.get("ranking_score"), -10**9),
            -safe_int(r.get("peak_count"), 0),
            -safe_int(r.get("file_size_bytes"), 0),
            normalize_text(r.get("file_accession")),
        ),
    )

    best = scored_rows[0]
    return best, scored_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select one primary peak file per experiment."
    )
    parser.add_argument(
        "--files-tsv",
        default=None,
        help="Override metadata/encode/files.tsv",
    )
    parser.add_argument(
        "--manifest-tsv",
        default=None,
        help="Override metadata/manifest/download_manifest.tsv",
    )
    parser.add_argument(
        "--inventory-tsv",
        default=None,
        help="Override processed/peaks/inventory/peak_file_inventory.tsv",
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

    files_tsv = (
        Path(args.files_tsv).resolve()
        if args.files_tsv
        else root / "metadata" / "encode" / "files.tsv"
    )
    manifest_tsv = (
        Path(args.manifest_tsv).resolve()
        if args.manifest_tsv
        else root / "metadata" / "manifest" / "download_manifest.tsv"
    )
    inventory_tsv = (
        Path(args.inventory_tsv).resolve()
        if args.inventory_tsv
        else root / "processed" / "peaks" / "inventory" / "peak_file_inventory.tsv"
    )

    out_dir = root / "processed" / "peaks" / "selected"
    primary_out = out_dir / "primary_peak_files.tsv"
    ranking_out = out_dir / "peak_file_ranking.tsv"
    log_path = root / "logs" / "select_primary_peak_files.log"

    log(f"[INFO] Project root: {root}", log_path)
    log(f"[INFO] files.tsv: {files_tsv}", log_path)
    log(f"[INFO] manifest.tsv: {manifest_tsv}", log_path)
    log(f"[INFO] inventory.tsv: {inventory_tsv}", log_path)

    file_rows = read_tsv_rows(files_tsv)
    manifest_rows = read_tsv_rows(manifest_tsv)
    inventory_rows = read_tsv_rows(inventory_tsv)

    file_idx = build_index(file_rows, "file_accession")
    inventory_idx = build_index(inventory_rows, "file_accession")

    target_filter = {x.strip().lower() for x in args.targets if x.strip()} if args.targets else None
    biosample_filter = {x.strip().lower() for x in args.biosamples if x.strip()} if args.biosamples else None

    candidates_by_experiment: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    for m in manifest_rows:
        file_accession = normalize_text(m.get("file_accession"))
        exp_acc = normalize_text(m.get("experiment_accession"))
        target_label = normalize_text(m.get("target_label"))
        biosample = normalize_text(m.get("biosample_term_name"))

        if not file_accession or not exp_acc:
            continue

        if target_filter and target_label.lower() not in target_filter:
            continue
        if biosample_filter and biosample.lower() not in biosample_filter:
            continue

        f = file_idx.get(file_accession, {})
        inv = inventory_idx.get(file_accession, {})
        local_path = normalize_text(m.get("local_path"))
        abs_path = root / local_path

        peak_count = 0
        if abs_path.exists() and normalize_lower(inv.get("status")) == "ok":
            try:
                peak_count = count_peak_lines(abs_path)
            except Exception as e:
                log(f"[WARN] Failed to count peaks for {file_accession}: {e}", log_path)
                peak_count = 0

        candidate = {
            "experiment_accession": exp_acc,
            "target_label": target_label,
            "biosample_term_name": biosample,
            "file_accession": file_accession,
            "output_type": normalize_text(f.get("output_type", m.get("output_type", ""))),
            "file_format": normalize_text(f.get("file_format", m.get("file_format", ""))),
            "assembly": normalize_text(f.get("assembly", m.get("assembly", ""))),
            "local_path": local_path,
            "file_size_bytes": safe_int(inv.get("file_size_bytes", 0)),
            "peak_count": peak_count,
            "inventory_status": normalize_text(inv.get("status", "")),
        }

        candidates_by_experiment[exp_acc].append(candidate)

    if not candidates_by_experiment:
        log("[WARN] No candidate peak files found after filtering.", log_path)
        return 0

    primary_rows: List[Dict[str, Any]] = []
    ranking_rows: List[Dict[str, Any]] = []

    for exp_acc in sorted(candidates_by_experiment.keys()):
        candidates = candidates_by_experiment[exp_acc]
        if not candidates:
            continue

        best, ranked = choose_best_candidate(candidates)

        for r in ranked:
            ranking_rows.append(
                {
                    "experiment_accession": r["experiment_accession"],
                    "target_label": r["target_label"],
                    "biosample_term_name": r["biosample_term_name"],
                    "file_accession": r["file_accession"],
                    "output_type": r["output_type"],
                    "file_format": r["file_format"],
                    "assembly": r["assembly"],
                    "local_path": r["local_path"],
                    "file_size_bytes": r["file_size_bytes"],
                    "peak_count": r["peak_count"],
                    "inventory_status": r["inventory_status"],
                    "ranking_score": r["ranking_score"],
                    "ranking_reason": r["ranking_reason"],
                    "is_selected": "yes" if r["file_accession"] == best["file_accession"] else "no",
                }
            )

        primary_rows.append(
            {
                "experiment_accession": best["experiment_accession"],
                "target_label": best["target_label"],
                "biosample_term_name": best["biosample_term_name"],
                "selected_file_accession": best["file_accession"],
                "selected_output_type": best["output_type"],
                "selected_file_format": best["file_format"],
                "selected_assembly": best["assembly"],
                "selected_local_path": best["local_path"],
                "selected_file_size_bytes": best["file_size_bytes"],
                "selected_peak_count": best["peak_count"],
                "selection_score": best["ranking_score"],
                "selection_reason": best["ranking_reason"],
            }
        )

    primary_rows = sorted(
        primary_rows,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
        ),
    )

    ranking_rows = sorted(
        ranking_rows,
        key=lambda r: (
            normalize_text(r.get("target_label")),
            normalize_text(r.get("biosample_term_name")),
            normalize_text(r.get("experiment_accession")),
            -safe_int(r.get("ranking_score"), -10**9),
            -safe_int(r.get("peak_count"), 0),
            normalize_text(r.get("file_accession")),
        ),
    )

    write_tsv(primary_out, primary_rows, PRIMARY_COLUMNS)
    write_tsv(ranking_out, ranking_rows, RANKING_COLUMNS)

    log("", log_path)
    log("[DONE] Primary peak selection completed.", log_path)
    log(f"[DONE] Primary output: {primary_out}", log_path)
    log(f"[DONE] Ranking output: {ranking_out}", log_path)
    log("", log_path)

    for row in primary_rows:
        log(
            "  "
            + " | ".join(
                [
                    normalize_text(row.get("experiment_accession")),
                    normalize_text(row.get("target_label")),
                    normalize_text(row.get("biosample_term_name")),
                    normalize_text(row.get("selected_file_accession")),
                    normalize_text(row.get("selected_output_type")),
                    f"peak_count={normalize_text(row.get('selected_peak_count'))}",
                    f"score={normalize_text(row.get('selection_score'))}",
                ]
            ),
            log_path,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())