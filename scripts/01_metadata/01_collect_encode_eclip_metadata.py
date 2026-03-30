#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Collect ENCODE eCLIP metadata for seed RBPs and optional biosamples.

Outputs:
- metadata/encode/experiments.tsv
- metadata/encode/files.tsv
- metadata/encode/audit.tsv

Recommended usage:
    cd ~/wangshuo/rna_dynamics
    python scripts/01_metadata/01_collect_encode_eclip_metadata.py

Optional usage:
    python scripts/01_metadata/01_collect_encode_eclip_metadata.py \
        --rbps TARDBP HNRNPK \
        --biosamples HepG2 K562

Notes:
- Default filters: assay_title=eCLIP, organism=Homo sapiens, status=released
- The script reads:
    metadata/rbp/rbp_seed_list.tsv
    metadata/biosample/biosample_list.tsv   (optional)
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import requests


BASE_URL = "https://www.encodeproject.org"
HEADERS = {
    "accept": "application/json",
    "user-agent": "rna_dynamics_encode_metadata_collector/1.0",
}

DEFAULT_ASSAY = "eCLIP"
DEFAULT_ORGANISM = "Homo sapiens"
DEFAULT_STATUS = "released"
REQUEST_TIMEOUT = 90
REQUEST_RETRY = 3
REQUEST_SLEEP = 1.5


EXPERIMENT_COLUMNS = [
    "experiment_accession",
    "target_label",
    "assay_title",
    "biosample_term_name",
    "biosample_type",
    "organism",
    "lab",
    "status",
    "audit_warning",
    "qc_keep",
    "qc_reason",
]

FILE_COLUMNS = [
    "file_accession",
    "experiment_accession",
    "file_format",
    "output_type",
    "assembly",
    "biological_replicates",
    "technical_replicates",
    "download_url",
]

AUDIT_COLUMNS = [
    "experiment_accession",
    "file_accession",
    "audit_level",
    "audit_category",
    "message",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def safe_join(items: Iterable[Any], sep: str = ";") -> str:
    vals = []
    for x in items:
        if x is None:
            continue
        s = str(x).strip()
        if s:
            vals.append(s)
    return sep.join(vals)


def ensure_project_root() -> Path:
    """
    Find the project root based on current file location:
    rna_dynamics/scripts/01_metadata/01_collect_encode_eclip_metadata.py
    """
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    return root


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


def load_seed_rbps(path: Path) -> List[str]:
    rows = read_tsv_rows(path)
    rbps: List[str] = []
    for row in rows:
        rbp = (row.get("rbp_symbol") or "").strip()
        if rbp:
            rbps.append(rbp)
    return sorted(set(rbps))


def load_biosamples(path: Path) -> List[str]:
    rows = read_tsv_rows(path)
    biosamples: List[str] = []
    for row in rows:
        bs = (row.get("biosample_term_name") or "").strip()
        if bs:
            biosamples.append(bs)
    return sorted(set(biosamples))


def request_json(url: str, params: Dict[str, Any]) -> Dict[str, Any]:
    last_err: Optional[Exception] = None
    for attempt in range(1, REQUEST_RETRY + 1):
        try:
            resp = requests.get(
                url,
                params=params,
                headers=HEADERS,
                timeout=REQUEST_TIMEOUT,
            )
            resp.raise_for_status()
            return resp.json()
        except Exception as e:
            last_err = e
            log(f"[WARN] Request failed (attempt {attempt}/{REQUEST_RETRY}): {e}")
            if attempt < REQUEST_RETRY:
                time.sleep(REQUEST_SLEEP * attempt)
    raise RuntimeError(f"Failed request after {REQUEST_RETRY} attempts: {url} | params={params}") from last_err


def encode_search_experiments(
    rbp: str,
    biosample: Optional[str],
    assay_title: str = DEFAULT_ASSAY,
    organism: str = DEFAULT_ORGANISM,
    status: str = DEFAULT_STATUS,
    limit: str = "all",
) -> List[Dict[str, Any]]:
    """
    Search ENCODE experiment objects.
    """
    url = f"{BASE_URL}/search/"
    params: Dict[str, Any] = {
        "type": "Experiment",
        "assay_title": assay_title,
        "status": status,
        "target.label": rbp,
        "replicates.library.biosample.donor.organism.scientific_name": organism,
        "format": "json",
        "limit": limit,
        "frame": "embedded",
    }
    if biosample:
        params["biosample_ontology.term_name"] = biosample

    data = request_json(url, params)
    return data.get("@graph", [])


def extract_target_label(exp: Dict[str, Any]) -> str:
    target = exp.get("target")
    if isinstance(target, dict):
        return (target.get("label") or "").strip()
    return ""


def extract_lab_title(exp: Dict[str, Any]) -> str:
    lab = exp.get("lab")
    if isinstance(lab, dict):
        return (lab.get("title") or "").strip()
    return ""


def extract_audit_warning(exp: Dict[str, Any]) -> str:
    """
    Summarize audit levels from ENCODE object.
    """
    audits = exp.get("audit", {}) or {}
    found: List[str] = []

    if isinstance(audits, dict):
        for level, items in audits.items():
            if not items:
                continue
            if isinstance(items, list) and len(items) > 0:
                found.append(level)

    # Fallback: if audit is not embedded as dict, inspect possible aliases
    if not found and exp.get("audit"):
        found.append("present")

    return safe_join(sorted(set(found)))


def infer_qc_keep_and_reason(exp: Dict[str, Any]) -> Tuple[str, str]:
    """
    Conservative initial QC label:
    - keep=yes if released and no obvious severe audit keyword
    - keep=manual_check when warning/error/not_compliant exists
    """
    status = (exp.get("status") or "").strip().lower()
    audits = exp.get("audit", {}) or {}

    if status != "released":
        return "no", "status_not_released"

    flagged_levels = set()
    if isinstance(audits, dict):
        for level, items in audits.items():
            if items:
                flagged_levels.add(level.lower())

    severe_terms = {"error", "not compliant", "not_compliant", "warning"}
    if any(level in severe_terms for level in flagged_levels):
        return "manual_check", f"audit_flag:{safe_join(sorted(flagged_levels), sep=',')}"

    possible_controls = exp.get("possible_controls", [])
    if not possible_controls:
        return "manual_check", "missing_or_unembedded_control_info"

    return "yes", "released_no_obvious_severe_audit"


def collect_replicate_biosample_info(exp: Dict[str, Any]) -> Tuple[str, str, str]:
    """
    Return:
    - biosample_term_name: semicolon-joined
    - biosample_type: semicolon-joined classification(s)
    - organism: semicolon-joined scientific names
    """
    biosample_names: Set[str] = set()
    biosample_types: Set[str] = set()
    organisms: Set[str] = set()

    for rep in exp.get("replicates", []) or []:
        try:
            biosample = rep["library"]["biosample"]
            bo = biosample.get("biosample_ontology", {}) or {}
            term_name = (bo.get("term_name") or "").strip()
            classification = (bo.get("classification") or "").strip()
            donor = biosample.get("donor", {}) or {}
            organism = (donor.get("organism", {}).get("scientific_name") or "").strip()

            if term_name:
                biosample_names.add(term_name)
            if classification:
                biosample_types.add(classification)
            if organism:
                organisms.add(organism)
        except Exception:
            continue

    return (
        safe_join(sorted(biosample_names)),
        safe_join(sorted(biosample_types)),
        safe_join(sorted(organisms)),
    )


def build_experiment_row(exp: Dict[str, Any]) -> Dict[str, Any]:
    biosample_term_name, biosample_type, organism = collect_replicate_biosample_info(exp)
    qc_keep, qc_reason = infer_qc_keep_and_reason(exp)

    return {
        "experiment_accession": (exp.get("accession") or "").strip(),
        "target_label": extract_target_label(exp),
        "assay_title": (exp.get("assay_title") or "").strip(),
        "biosample_term_name": biosample_term_name,
        "biosample_type": biosample_type,
        "organism": organism,
        "lab": extract_lab_title(exp),
        "status": (exp.get("status") or "").strip(),
        "audit_warning": extract_audit_warning(exp),
        "qc_keep": qc_keep,
        "qc_reason": qc_reason,
    }


def normalize_replicate_list(values: Any) -> str:
    if not values:
        return ""
    if isinstance(values, list):
        return safe_join(values)
    return str(values)


def make_download_url(href: str) -> str:
    href = (href or "").strip()
    if not href:
        return ""
    if href.startswith("http://") or href.startswith("https://"):
        return href
    return f"{BASE_URL}{href}"


def build_file_rows(exp: Dict[str, Any]) -> List[Dict[str, Any]]:
    exp_acc = (exp.get("accession") or "").strip()
    rows: List[Dict[str, Any]] = []

    for f in exp.get("files", []) or []:
        row = {
            "file_accession": (f.get("accession") or "").strip(),
            "experiment_accession": exp_acc,
            "file_format": (f.get("file_format") or "").strip(),
            "output_type": (f.get("output_type") or "").strip(),
            "assembly": (f.get("assembly") or "").strip(),
            "biological_replicates": normalize_replicate_list(f.get("biological_replicates")),
            "technical_replicates": normalize_replicate_list(f.get("technical_replicates")),
            "download_url": make_download_url(f.get("href") or ""),
        }
        rows.append(row)

    return rows


def build_audit_rows(exp: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Collect audit rows from experiment and files.
    """
    rows: List[Dict[str, Any]] = []
    exp_acc = (exp.get("accession") or "").strip()

    # Experiment-level audit
    exp_audit = exp.get("audit", {}) or {}
    if isinstance(exp_audit, dict):
        for level, items in exp_audit.items():
            if not items:
                continue
            for item in items:
                category = ""
                message = ""
                if isinstance(item, dict):
                    category = (item.get("category") or "").strip()
                    message = (item.get("detail") or item.get("message") or "").strip()
                rows.append({
                    "experiment_accession": exp_acc,
                    "file_accession": "",
                    "audit_level": str(level),
                    "audit_category": category,
                    "message": message,
                })

    # File-level audit
    for f in exp.get("files", []) or []:
        file_acc = (f.get("accession") or "").strip()
        file_audit = f.get("audit", {}) or {}
        if isinstance(file_audit, dict):
            for level, items in file_audit.items():
                if not items:
                    continue
                for item in items:
                    category = ""
                    message = ""
                    if isinstance(item, dict):
                        category = (item.get("category") or "").strip()
                        message = (item.get("detail") or item.get("message") or "").strip()
                    rows.append({
                        "experiment_accession": exp_acc,
                        "file_accession": file_acc,
                        "audit_level": str(level),
                        "audit_category": category,
                        "message": message,
                    })

    return rows


def deduplicate_rows(rows: List[Dict[str, Any]], key_fields: List[str]) -> List[Dict[str, Any]]:
    seen = set()
    out = []
    for row in rows:
        key = tuple(row.get(k, "") for k in key_fields)
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect ENCODE eCLIP metadata for seed RBPs and optional biosamples."
    )
    parser.add_argument(
        "--rbps",
        nargs="*",
        default=None,
        help="Override RBP list from metadata/rbp/rbp_seed_list.tsv",
    )
    parser.add_argument(
        "--biosamples",
        nargs="*",
        default=None,
        help="Optional biosample filters. If omitted, use metadata/biosample/biosample_list.tsv if present; otherwise no biosample restriction.",
    )
    parser.add_argument(
        "--assay-title",
        default=DEFAULT_ASSAY,
        help=f"Assay title to search (default: {DEFAULT_ASSAY})",
    )
    parser.add_argument(
        "--organism",
        default=DEFAULT_ORGANISM,
        help=f"Organism scientific name (default: {DEFAULT_ORGANISM})",
    )
    parser.add_argument(
        "--status",
        default=DEFAULT_STATUS,
        help=f"Experiment status (default: {DEFAULT_STATUS})",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.4,
        help="Sleep seconds between queries to avoid hammering the API (default: 0.4)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = ensure_project_root()

    rbp_seed_path = root / "metadata" / "rbp" / "rbp_seed_list.tsv"
    biosample_list_path = root / "metadata" / "biosample" / "biosample_list.tsv"

    out_exp_path = root / "metadata" / "encode" / "experiments.tsv"
    out_file_path = root / "metadata" / "encode" / "files.tsv"
    out_audit_path = root / "metadata" / "encode" / "audit.tsv"

    if args.rbps:
        rbps = sorted(set([x.strip() for x in args.rbps if x.strip()]))
    else:
        rbps = load_seed_rbps(rbp_seed_path)

    if not rbps:
        log("[ERROR] No RBP seeds found.")
        log(f"[ERROR] Please fill: {rbp_seed_path}")
        log("[ERROR] Example:")
        log("rbp_symbol\tpriority\treason")
        log("TARDBP\t1\tclassic RBP example")
        log("HNRNPK\t1\tclassic hnRNP family")
        return 1

    if args.biosamples is not None:
        biosamples = sorted(set([x.strip() for x in args.biosamples if x.strip()]))
    else:
        biosamples = load_biosamples(biosample_list_path)

    use_biosample_filter = len(biosamples) > 0

    log(f"[INFO] Project root: {root}")
    log(f"[INFO] RBP count: {len(rbps)} -> {', '.join(rbps)}")
    if use_biosample_filter:
        log(f"[INFO] Biosample filter count: {len(biosamples)} -> {', '.join(biosamples)}")
    else:
        log("[INFO] No biosample filter provided; searching all matching biosamples.")
    log(f"[INFO] Search filters: assay={args.assay_title}, organism={args.organism}, status={args.status}")

    experiment_rows: List[Dict[str, Any]] = []
    file_rows: List[Dict[str, Any]] = []
    audit_rows: List[Dict[str, Any]] = []

    total_queries = 0

    for rbp in rbps:
        query_biosamples = biosamples if use_biosample_filter else [None]

        for biosample in query_biosamples:
            total_queries += 1
            label = biosample if biosample else "ALL_BIOSAMPLES"
            log(f"[INFO] Querying ENCODE: RBP={rbp} | biosample={label}")

            try:
                experiments = encode_search_experiments(
                    rbp=rbp,
                    biosample=biosample,
                    assay_title=args.assay_title,
                    organism=args.organism,
                    status=args.status,
                )
            except Exception as e:
                log(f"[ERROR] Query failed for RBP={rbp}, biosample={label}: {e}")
                continue

            log(f"[INFO] Retrieved {len(experiments)} experiments for RBP={rbp}, biosample={label}")

            for exp in experiments:
                experiment_rows.append(build_experiment_row(exp))
                file_rows.extend(build_file_rows(exp))
                audit_rows.extend(build_audit_rows(exp))

            time.sleep(args.sleep)

    experiment_rows = deduplicate_rows(experiment_rows, ["experiment_accession"])
    file_rows = deduplicate_rows(file_rows, ["file_accession"])
    audit_rows = deduplicate_rows(
        audit_rows,
        ["experiment_accession", "file_accession", "audit_level", "audit_category", "message"],
    )

    # Sort for stable outputs
    experiment_rows = sorted(
        experiment_rows,
        key=lambda x: (
            x.get("target_label", ""),
            x.get("biosample_term_name", ""),
            x.get("experiment_accession", ""),
        ),
    )
    file_rows = sorted(
        file_rows,
        key=lambda x: (
            x.get("experiment_accession", ""),
            x.get("output_type", ""),
            x.get("file_accession", ""),
        ),
    )
    audit_rows = sorted(
        audit_rows,
        key=lambda x: (
            x.get("experiment_accession", ""),
            x.get("file_accession", ""),
            x.get("audit_level", ""),
            x.get("audit_category", ""),
        ),
    )

    write_tsv(out_exp_path, experiment_rows, EXPERIMENT_COLUMNS)
    write_tsv(out_file_path, file_rows, FILE_COLUMNS)
    write_tsv(out_audit_path, audit_rows, AUDIT_COLUMNS)

    log("")
    log("[DONE] ENCODE metadata collection completed.")
    log(f"[DONE] Queries sent: {total_queries}")
    log(f"[DONE] Experiments written: {len(experiment_rows)} -> {out_exp_path}")
    log(f"[DONE] Files written: {len(file_rows)} -> {out_file_path}")
    log(f"[DONE] Audit rows written: {len(audit_rows)} -> {out_audit_path}")

    # Quick next-step summary
    peak_like = [r for r in file_rows if "peak" in (r.get("output_type", "").lower()) or r.get("file_format", "").lower() in {"bed", "narrowpeak", "broadpeak"}]
    bigwig_like = [r for r in file_rows if r.get("file_format", "").lower() == "bigwig"]
    bam_like = [r for r in file_rows if r.get("file_format", "").lower() == "bam"]

    log("")
    log("[INFO] Quick file summary:")
    log(f"  Peak-like files : {len(peak_like)}")
    log(f"  bigWig files    : {len(bigwig_like)}")
    log(f"  BAM files       : {len(bam_like)}")

    log("")
    log("[INFO] Recommended next steps:")
    log("  1. Inspect metadata/encode/experiments.tsv")
    log("  2. Inspect metadata/encode/files.tsv")
    log("  3. Build metadata/manifest/download_manifest.tsv from peak-like files first")
    log("  4. Then write scripts/02_download/download_encode_data.py")

    return 0


if __name__ == "__main__":
    sys.exit(main())