#!/usr/bin/env python3
"""Audit GNPS PNNL lipid libraries for mixed-resolution learning feasibility."""

from __future__ import annotations

import argparse
import collections
import gzip
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
PREFERRED = ("PNNL-LIPIDS-POSITIVE", "PNNL-LIPIDS-NEGATIVE")
ARITY = {
    "PC": 2, "PE": 2, "PG": 2, "PI": 2, "PS": 2, "PA": 2,
    "DG": 2, "TG": 3,
}
CHAIN_RE = re.compile(r"(?<![A-Za-z0-9])(?:[dmt]|[OP]-)?(\d{1,3}):(\d{1,2})(?!\d)")
CLASS_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9-]*)")
ADDUCT_RE = re.compile(r"\[(?:M|m)[^\]]+\][+-]?")
EXCLUDED_CHEMISTRY_RE = re.compile(r"(?:[OP]-\d|\b[dmt]\d{1,3}:\d)", re.I)


def fetch(url: str, retries: int = 4) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "resolution-lattice-audit/0.1"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=180) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to download {url}: {last}")


def field(record: dict, *names: str) -> str:
    lowered = {str(k).lower(): v for k, v in record.items()}
    for name in names:
        value = lowered.get(name.lower())
        if value is not None:
            return str(value).strip()
    return ""


def discover_libraries() -> list[str]:
    try:
        obj = json.loads(fetch(f"{BASE}.json"))
    except Exception:
        return list(PREFERRED)
    names = []
    rows = obj if isinstance(obj, list) else obj.get("libraries", obj.get("data", []))
    if isinstance(rows, list):
        for row in rows:
            if isinstance(row, str):
                name = row
            elif isinstance(row, dict):
                name = field(row, "libraryname", "library_name", "name")
            else:
                continue
            if "PNNL-LIPIDS" in name.upper():
                names.append(name)
    ordered = [name for name in PREFERRED if name in names]
    ordered.extend(sorted(set(names) - set(ordered)))
    return ordered or list(PREFERRED)


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
    chains = [(int(c), int(u)) for c, u in CHAIN_RE.findall(cleaned)]
    needed = ARITY[lipid_class]
    if len(chains) == needed:
        species = (lipid_class, tuple(sorted(chains)))
        sum_key = (lipid_class, sum(c for c, _ in chains), sum(u for _, u in chains))
        return {
            "status": "fine",
            "class": lipid_class,
            "species": species,
            "sum_key": sum_key,
            "name": name,
        }
    if len(chains) == 1 and needed > 1:
        c, u = chains[0]
        return {
            "status": "coarse",
            "class": lipid_class,
            "sum_key": (lipid_class, c, u),
            "name": name,
        }
    return {"status": "unparsed", "class": lipid_class, "name": name}


def entropy(counts):
    total = sum(counts)
    if not total:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counts if n)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/lipid_lattice_audit")
    args = parser.parse_args()
    outdir = pathlib.Path(args.outdir)
    rawdir = outdir / "raw"
    rawdir.mkdir(parents=True, exist_ok=True)

    source_manifest = []
    records = []
    for library in discover_libraries():
        url = f"{BASE}/{library}.json"
        raw = fetch(url)
        source_manifest.append({
            "library": library,
            "url": url,
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        })
        with gzip.open(rawdir / f"{library}.json.gz", "wb", compresslevel=9) as handle:
            handle.write(raw)
        payload = json.loads(raw)
        rows = payload if isinstance(payload, list) else payload.get("spectra", payload.get("data", []))
        if not isinstance(rows, list):
            raise TypeError(f"unexpected payload for {library}: {type(rows)!r}")
        for row in rows:
            if isinstance(row, dict):
                row["_audit_library"] = library
                records.append(row)

    status_counts = collections.Counter()
    class_counts = collections.Counter()
    fine_spectrum_counts = collections.Counter()
    fine_label_examples = {}
    coarse_counts = collections.Counter()
    unparsed_examples = []
    instrument_counts = collections.Counter()
    ion_counts = collections.Counter()

    for record in records:
        name = field(record, "Compound_Name", "compound_name", "name", "NAME")
        parsed = parse_label(name)
        if parsed is None:
            status_counts["not_target_class"] += 1
            continue
        status = parsed["status"]
        status_counts[status] += 1
        class_counts[parsed["class"]] += 1
        instrument_counts[field(record, "Instrument", "instrument") or "missing"] += 1
        ion_counts[field(record, "Ion_Mode", "ion_mode", "ionmode") or "missing"] += 1
        if status == "fine":
            key = (parsed["sum_key"], parsed["species"])
            fine_spectrum_counts[key] += 1
            fine_label_examples.setdefault(parsed["species"], name)
        elif status == "coarse":
            coarse_counts[parsed["sum_key"]] += 1
        elif status == "unparsed" and len(unparsed_examples) < 100:
            unparsed_examples.append(name)

    by_sum = collections.defaultdict(collections.Counter)
    for (sum_key, species), count in fine_spectrum_counts.items():
        by_sum[sum_key][species] += count

    fine_spectra = sum(fine_spectrum_counts.values())
    unique_species = len({species for _, species in fine_spectrum_counts})
    ambiguous = {key: counts for key, counts in by_sum.items() if len(counts) >= 2}
    ambiguous_spectra = sum(sum(counts.values()) for counts in ambiguous.values())
    candidate_counts = [len(counts) for counts in ambiguous.values()]
    weighted_entropy = (
        sum(sum(counts.values()) * entropy(counts.values()) for counts in by_sum.values())
        / fine_spectra if fine_spectra else 0.0
    )
    matched_coarse_spectra = sum(n for key, n in coarse_counts.items() if key in by_sum)
    top_ambiguous = []
    for key, counts in sorted(
        ambiguous.items(), key=lambda item: (-len(item[1]), -sum(item[1].values()), item[0])
    )[:100]:
        top_ambiguous.append({
            "sum_key": list(key),
            "candidate_count": len(counts),
            "spectra": sum(counts.values()),
            "candidate_examples": [
                fine_label_examples.get(species, str(species))
                for species, _ in counts.most_common(20)
            ],
        })

    metrics = {
        "record_count": len(records),
        "status_counts": dict(status_counts),
        "class_counts": dict(class_counts),
        "fine_spectra": fine_spectra,
        "unique_fine_species": unique_species,
        "sum_composition_classes": len(by_sum),
        "ambiguous_sum_classes": len(ambiguous),
        "ambiguous_spectrum_fraction": ambiguous_spectra / fine_spectra if fine_spectra else 0.0,
        "ambiguous_candidate_count_median": statistics.median(candidate_counts) if candidate_counts else 0,
        "ambiguous_candidate_count_mean": statistics.mean(candidate_counts) if candidate_counts else 0,
        "ambiguous_candidate_count_max": max(candidate_counts, default=0),
        "weighted_conditional_entropy_bits": weighted_entropy,
        "coarse_spectra": sum(coarse_counts.values()),
        "coarse_sum_classes": len(coarse_counts),
        "coarse_spectra_with_known_candidate_set": matched_coarse_spectra,
        "instrument_counts_top20": instrument_counts.most_common(20),
        "ion_mode_counts": dict(ion_counts),
        "top_ambiguous_classes": top_ambiguous,
        "unparsed_examples": unparsed_examples,
    }
    gates = {
        "fine_spectra_ge_20000": fine_spectra >= 20000,
        "unique_fine_species_ge_750": unique_species >= 750,
        "ambiguous_sum_classes_ge_250": len(ambiguous) >= 250,
        "ambiguous_spectrum_fraction_ge_0_30": (
            ambiguous_spectra / fine_spectra >= 0.30 if fine_spectra else False
        ),
        "ambiguous_candidate_median_ge_2": (
            statistics.median(candidate_counts) >= 2 if candidate_counts else False
        ),
        "weighted_entropy_ge_0_50_bits": weighted_entropy >= 0.50,
    }
    result = {
        "audit_version": "0.1.0",
        "sources": source_manifest,
        "metrics": metrics,
        "gates": gates,
        "numeric_gate_pass": all(gates.values()),
        "literature_gate": "manual_pending",
    }
    (outdir / "audit.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (outdir / "source_manifest.json").write_text(
        json.dumps(source_manifest, indent=2), encoding="utf-8"
    )

    lines = [
        "# Lipid lattice feasibility audit",
        "",
        f"- Records downloaded: {len(records):,}",
        f"- Fine spectra: {fine_spectra:,}",
        f"- Unique fine species: {unique_species:,}",
        f"- Sum-composition classes: {len(by_sum):,}",
        f"- Ambiguous sum classes: {len(ambiguous):,}",
        f"- Fine spectra in ambiguous classes: {(ambiguous_spectra / fine_spectra if fine_spectra else 0):.4%}",
        f"- Candidate-count median/mean/max: {statistics.median(candidate_counts) if candidate_counts else 0} / {statistics.mean(candidate_counts) if candidate_counts else 0:.3f} / {max(candidate_counts, default=0)}",
        f"- Weighted H(fine | class,sum): {weighted_entropy:.4f} bits",
        f"- Native coarse spectra: {sum(coarse_counts.values()):,}",
        f"- Native coarse spectra with known candidate set: {matched_coarse_spectra:,}",
        "",
        "## Fixed numerical gates",
        "",
    ]
    lines.extend(f"- [{'x' if passed else ' '}] {name}" for name, passed in gates.items())
    lines.extend([
        "",
        f"Overall numerical gate: {'PASS' if all(gates.values()) else 'FAIL'}",
        "",
        "The separate literature novelty gate remains manual and is not implied by this audit.",
    ])
    (outdir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
