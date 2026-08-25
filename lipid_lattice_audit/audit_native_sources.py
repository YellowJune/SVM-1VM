#!/usr/bin/env python3
"""Audit preregistered native coarse-resolution lipid MS/MS availability."""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import pathlib
import re
import statistics
import time
import urllib.error
import urllib.request

BASE = "https://external.gnps2.org/gnpslibrary"
LIBRARIES = (
    "PNNL-LIPIDS-POSITIVE",
    "PNNL-LIPIDS-NEGATIVE",
    "HCE-CELL-LYSATE-LIPIDS",
    "GNPS-D2-AMINO-LIPID-LIBRARY",
    "GNPS-N-ACYL-LIPIDS-MASSQL",
    "GNPS-LIPID-MAPS-STANDARDS-SPECTRA-DB",
)
PNNL = {"PNNL-LIPIDS-POSITIVE", "PNNL-LIPIDS-NEGATIVE"}
ARITY = {
    "PC": 2,
    "PE": 2,
    "PG": 2,
    "PI": 2,
    "PS": 2,
    "PA": 2,
    "DG": 2,
    "TG": 3,
}
CHAIN_RE = re.compile(r"(?<![A-Za-z0-9])(?:[dmt]|[OP]-)?(\d{1,3}):(\d{1,2})(?!\d)")
CLASS_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9-]*)")
ADDUCT_RE = re.compile(r"\[(?:M|m)[^\]]+\][+-]?")
EXCLUDED_CHEMISTRY_RE = re.compile(r"(?:[OP]-\d|\b[dmt]\d{1,3}:\d)", re.I)


def fetch(url: str, retries: int = 4) -> bytes:
    request = urllib.request.Request(
        url, headers={"User-Agent": "native-lipid-resolution-audit/0.1"}
    )
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=240) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {url}: {last}")


def field(record: dict, *names: str):
    lowered = {str(key).lower(): value for key, value in record.items()}
    for name in names:
        value = lowered.get(name.lower())
        if value is not None and str(value).strip():
            return value
    return ""


def parse_label(name: str):
    cleaned = ADDUCT_RE.sub("", name).strip()
    match = CLASS_RE.search(cleaned)
    if not match:
        return None
    lipid_class = match.group(1).upper()
    if lipid_class not in ARITY:
        return None
    if EXCLUDED_CHEMISTRY_RE.search(cleaned):
        return {"status": "excluded_chemistry", "class": lipid_class, "name": name}
    chains = [(int(carbon), int(double)) for carbon, double in CHAIN_RE.findall(cleaned)]
    arity = ARITY[lipid_class]
    if len(chains) == arity:
        species = (lipid_class, tuple(sorted(chains)))
        return {
            "status": "fine",
            "class": lipid_class,
            "species": species,
            "sum_key": (
                lipid_class,
                sum(carbon for carbon, _ in chains),
                sum(double for _, double in chains),
            ),
            "name": name,
        }
    if len(chains) == 1 and arity > 1:
        carbon, double = chains[0]
        return {
            "status": "coarse",
            "class": lipid_class,
            "sum_key": (lipid_class, carbon, double),
            "name": name,
        }
    return {"status": "unparsed", "class": lipid_class, "name": name}


def stable_json(value):
    try:
        if isinstance(value, str):
            value = json.loads(value)
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    except Exception:
        return str(value).strip()


def spectrum_signature(record: dict) -> str:
    splash = field(record, "splash", "SPLASH")
    if splash:
        return "splash:" + str(splash).strip().lower()
    peaks = field(
        record,
        "peaks_json",
        "peaks",
        "spectrum",
        "mz_intensity",
        "MZ_INTENSITY",
    )
    precursor = field(
        record,
        "Precursor_MZ",
        "precursor_mz",
        "precursor",
        "EXACTMASS",
    )
    adduct = field(record, "Adduct", "adduct", "Ion_Source", "ion_source")
    ion_mode = field(record, "Ion_Mode", "ion_mode", "ionmode")
    if peaks:
        material = {
            "peaks": stable_json(peaks),
            "precursor": str(precursor),
            "adduct": str(adduct),
            "ion_mode": str(ion_mode),
        }
    else:
        material = {
            "spectrum_id": str(
                field(
                    record,
                    "spectrum_id",
                    "SpectrumID",
                    "CCMSLIBRARY",
                    "library_membership",
                )
            ),
            "precursor": str(precursor),
            "adduct": str(adduct),
            "ion_mode": str(ion_mode),
            "instrument": str(field(record, "Instrument", "instrument")),
        }
    return hashlib.sha256(
        json.dumps(material, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def entropy(values):
    counts = collections.Counter(values)
    total = sum(counts.values())
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/lipid_native_audit")
    args = parser.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifests = []
    records_by_library = {}
    for library in LIBRARIES:
        url = f"{BASE}/{library}.json"
        try:
            raw = fetch(url)
            payload = json.loads(raw)
            rows = payload if isinstance(payload, list) else payload.get(
                "spectra", payload.get("data", [])
            )
            if not isinstance(rows, list):
                raise TypeError(f"unexpected payload for {library}: {type(rows)!r}")
            rows = [row for row in rows if isinstance(row, dict)]
            records_by_library[library] = rows
            manifests.append(
                {
                    "library": library,
                    "url": url,
                    "bytes": len(raw),
                    "sha256": hashlib.sha256(raw).hexdigest(),
                    "records": len(rows),
                    "error": None,
                }
            )
        except Exception as exc:
            records_by_library[library] = []
            manifests.append(
                {
                    "library": library,
                    "url": url,
                    "bytes": 0,
                    "sha256": None,
                    "records": 0,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )

    pnnl_species = collections.defaultdict(set)
    pnnl_fine_counts = collections.Counter()
    for library in LIBRARIES:
        if library not in PNNL:
            continue
        for record in records_by_library[library]:
            name = str(field(record, "Compound_Name", "compound_name", "name", "NAME"))
            parsed = parse_label(name)
            if parsed and parsed["status"] == "fine":
                pnnl_species[parsed["sum_key"]].add(parsed["species"])
                pnnl_fine_counts[parsed["sum_key"]] += 1
    ambiguous_pnnl = {
        key: species for key, species in pnnl_species.items() if len(species) >= 2
    }

    status_by_library = {}
    coarse_raw_by_library = collections.Counter()
    native_by_signature = {}
    unparsed_examples = collections.defaultdict(list)
    for library in LIBRARIES:
        statuses = collections.Counter()
        for record in records_by_library[library]:
            name = str(field(record, "Compound_Name", "compound_name", "name", "NAME"))
            parsed = parse_label(name)
            if parsed is None:
                statuses["not_target_class"] += 1
                continue
            statuses[parsed["status"]] += 1
            if parsed["status"] == "coarse":
                coarse_raw_by_library[library] += 1
                signature = spectrum_signature(record)
                current = native_by_signature.get(signature)
                item = {
                    "library": library,
                    "sum_key": parsed["sum_key"],
                    "name": name,
                    "signature": signature,
                }
                if current is None:
                    native_by_signature[signature] = item
                elif current["sum_key"] != item["sum_key"]:
                    current.setdefault("label_conflicts", []).append(
                        {
                            "library": library,
                            "sum_key": item["sum_key"],
                            "name": name,
                        }
                    )
            elif parsed["status"] == "unparsed" and len(unparsed_examples[library]) < 30:
                unparsed_examples[library].append(name)
        status_by_library[library] = dict(statuses)

    native_unique = list(native_by_signature.values())
    mapped = [
        item for item in native_unique if item["sum_key"] in ambiguous_pnnl
    ]
    mapped_by_library = collections.Counter(item["library"] for item in mapped)
    mapped_by_class = collections.Counter(item["sum_key"][0] for item in mapped)
    mapped_keys = {item["sum_key"] for item in mapped}
    conflicts = sum(len(item.get("label_conflicts", [])) for item in native_unique)
    candidate_counts = [len(ambiguous_pnnl[item["sum_key"]]) for item in mapped]
    max_source_fraction = (
        max(mapped_by_library.values(), default=0) / len(mapped) if mapped else 1.0
    )

    metrics = {
        "downloaded_records": sum(len(rows) for rows in records_by_library.values()),
        "records_by_library": {
            library: len(records_by_library[library]) for library in LIBRARIES
        },
        "status_by_library": status_by_library,
        "raw_native_coarse_by_library": dict(coarse_raw_by_library),
        "unique_native_coarse_spectra": len(native_unique),
        "mapped_ambiguous_native_coarse_spectra": len(mapped),
        "mapped_fraction_of_unique_native_coarse": (
            len(mapped) / len(native_unique) if native_unique else 0.0
        ),
        "mapped_ambiguous_sum_keys": len(mapped_keys),
        "mapped_target_classes": len(mapped_by_class),
        "mapped_class_counts": dict(mapped_by_class),
        "mapped_contributing_libraries": len(mapped_by_library),
        "mapped_library_counts": dict(mapped_by_library),
        "largest_library_fraction": max_source_fraction,
        "candidate_count_median": (
            statistics.median(candidate_counts) if candidate_counts else 0
        ),
        "candidate_count_mean": (
            statistics.mean(candidate_counts) if candidate_counts else 0
        ),
        "cross_label_conflicts_after_spectral_dedup": conflicts,
        "pnnl_ambiguous_fine_keys": len(ambiguous_pnnl),
        "unparsed_examples": dict(unparsed_examples),
    }
    gates = {
        "native_coarse_spectra_ge_5000": len(native_unique) >= 5000,
        "mapped_ambiguous_sum_keys_ge_150": len(mapped_keys) >= 150,
        "mapped_target_classes_ge_4": len(mapped_by_class) >= 4,
        "mapped_fraction_ge_0_80": (
            len(mapped) / len(native_unique) >= 0.80 if native_unique else False
        ),
        "mapped_libraries_ge_2": len(mapped_by_library) >= 2,
        "largest_library_fraction_le_0_90": max_source_fraction <= 0.90,
        "all_sources_hashed": (
            len(manifests) == len(LIBRARIES)
            and all(
                item["error"] is None
                and item["bytes"] > 0
                and item["sha256"] is not None
                and len(item["sha256"]) == 64
                for item in manifests
            )
        ),
    }
    result = {
        "audit_version": "0.1.0",
        "preregistration_commit": "adbaa6bf3a8969e4cfac1bdd7da457463dcf7c2b",
        "sources": manifests,
        "metrics": metrics,
        "gates": gates,
        "numeric_gate_pass": all(gates.values()),
        "decision": "KEEP_FOR_METHOD_PROOF" if all(gates.values()) else "REJECT",
        "literature_gate": "manual_pending",
    }
    (outdir / "native_audit.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    (outdir / "source_manifest.json").write_text(
        json.dumps(manifests, indent=2), encoding="utf-8"
    )

    lines = [
        "# Native coarse-resolution lipid source audit",
        "",
        f"- Sources fetched and hashed: {len(manifests)}",
        f"- Records downloaded: {metrics['downloaded_records']:,}",
        f"- Unique native coarse spectra: {len(native_unique):,}",
        f"- Mapped ambiguous native coarse spectra: {len(mapped):,}",
        f"- Mapped fraction: {metrics['mapped_fraction_of_unique_native_coarse']:.4%}",
        f"- Mapped ambiguous sum keys: {len(mapped_keys):,}",
        f"- Mapped target classes: {len(mapped_by_class):,}",
        f"- Contributing libraries: {len(mapped_by_library):,}",
        f"- Largest-library fraction: {max_source_fraction:.4%}",
        f"- Cross-label conflicts after spectral deduplication: {conflicts:,}",
        "",
        "## Frozen gates",
        "",
    ]
    lines.extend(
        f"- [{'x' if passed else ' '}] {name}" for name, passed in gates.items()
    )
    lines.extend(
        [
            "",
            f"Decision: **{result['decision']}**",
            "",
            "Passing this audit would only authorize a method proof; it would not establish novelty.",
        ]
    )
    (outdir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
