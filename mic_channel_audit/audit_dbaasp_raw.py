#!/usr/bin/env python3
"""Audit raw DBAASP MIC labels without collapsing bounds to point targets.

This is a viability audit, not a model-training script. It samples peptide IDs
uniformly across the public DBAASP identifier range and reports how often MIC
annotations preserve censoring, ranges, and assay context.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures
import datetime as dt
import json
import math
import random
import re
import statistics
import time
from pathlib import Path
from typing import Any

import requests

VERSION = "0.1.0"
BASE_URL = "https://dbaasp.org/peptides/{identifier}"
DEFAULT_MAX_ID = 24_207
EMPTY = {None, "", "NA", "N/A", "-", "null", "None"}


def stratified_ids(max_id: int, sample_n: int) -> list[int]:
    sample_n = min(max_id, max(1, sample_n))
    ids = {max(1, min(max_id, int((i + 0.5) * max_id / sample_n))) for i in range(sample_n)}
    return sorted(ids)


def fetch_one(identifier: int, retries: int = 4) -> tuple[int, dict[str, Any] | None, str | None]:
    headers = {"Accept": "application/json", "User-Agent": "MIC-channel-audit/0.1 research"}
    for attempt in range(retries):
        try:
            response = requests.get(
                BASE_URL.format(identifier=identifier),
                headers=headers,
                timeout=(10, 30),
            )
            if response.status_code == 404:
                return identifier, None, "404"
            if response.status_code == 429 or response.status_code >= 500:
                wait = min(12.0, 0.7 * (2**attempt)) + random.random() * 0.4
                time.sleep(wait)
                continue
            response.raise_for_status()
            payload = response.json()
            return identifier, payload, None
        except (requests.RequestException, ValueError) as exc:
            if attempt + 1 == retries:
                return identifier, None, f"{type(exc).__name__}: {exc}"
            time.sleep(min(12.0, 0.7 * (2**attempt)) + random.random() * 0.4)
    return identifier, None, "retry_exhausted"


def nonempty(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip() not in EMPTY
    return value is not None


def flatten_leaves(value: Any, prefix: str = ""):
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            yield from flatten_leaves(child, path)
    elif isinstance(value, list):
        for child in value:
            yield from flatten_leaves(child, f"{prefix}[]")
    elif nonempty(value):
        yield prefix, value


def first_matching_leaf(activity: dict[str, Any], tokens: tuple[str, ...]) -> str | None:
    for path, value in flatten_leaves(activity):
        lower = path.lower()
        if any(token in lower for token in tokens):
            text = str(value).strip()
            if text and text not in EMPTY:
                return text
    return None


def species_name(activity: dict[str, Any]) -> str:
    species = activity.get("targetSpecies")
    if isinstance(species, dict) and nonempty(species.get("name")):
        return str(species["name"]).strip()
    return first_matching_leaf(activity, ("targetspecies", "species", "strain")) or "UNKNOWN"


def is_mic(activity: dict[str, Any]) -> bool:
    value = activity.get("activityMeasureValue")
    if isinstance(value, dict):
        value = value.get("name") or value.get("value")
    return str(value).strip().upper() == "MIC"


RIGHT_RE = re.compile(r"(?:>=|≥|>)")
LEFT_RE = re.compile(r"(?:<=|≤|<|\bup\s+to\b)", re.IGNORECASE)
RANGE_RE = re.compile(r"[-–—]")
PLUS_MINUS_RE = re.compile(r"(?:±|\+/-)")


def label_kind(raw: Any) -> str:
    if raw is None:
        return "missing"
    text = str(raw).strip()
    if text in EMPTY:
        return "missing"
    has_right = bool(RIGHT_RE.search(text))
    has_left = bool(LEFT_RE.search(text))
    if has_right and has_left:
        return "two_sided_censored"
    if has_right:
        return "right_censored"
    if has_left:
        return "left_censored"
    if PLUS_MINUS_RE.search(text):
        return "plus_minus"
    if RANGE_RE.search(text):
        return "finite_range"
    try:
        float(text.replace(",", ".").replace(" ", ""))
        return "exact_numeric"
    except ValueError:
        return "other"


def quantiles(values: list[int]) -> dict[str, float | int] | None:
    if not values:
        return None
    ordered = sorted(values)
    def q(p: float):
        index = p * (len(ordered) - 1)
        lo, hi = math.floor(index), math.ceil(index)
        if lo == hi:
            return ordered[lo]
        return ordered[lo] * (hi - index) + ordered[hi] * (index - lo)
    return {
        "min": ordered[0],
        "q25": q(0.25),
        "median": q(0.5),
        "q75": q(0.75),
        "max": ordered[-1],
        "mean": statistics.fmean(ordered),
    }


def serializable_counter(counter: collections.Counter, limit: int | None = None):
    items = counter.most_common(limit)
    return [{"value": str(key), "count": int(count)} for key, count in items]


def audit(records: list[dict[str, Any]], requested: int, missing_404: int, errors: list[str]):
    kind_counts: collections.Counter[str] = collections.Counter()
    unit_counts: collections.Counter[str] = collections.Counter()
    concentration_counts: collections.Counter[str] = collections.Counter()
    leaf_path_counts: collections.Counter[str] = collections.Counter()
    condition_presence: collections.Counter[str] = collections.Counter()
    medium_values: collections.Counter[str] = collections.Counter()
    cfu_values: collections.Counter[str] = collections.Counter()
    method_values: collections.Counter[str] = collections.Counter()
    reference_values: collections.Counter[str] = collections.Counter()
    pair_counts: collections.Counter[tuple[str, str]] = collections.Counter()
    pair_contexts: dict[tuple[str, str], set[tuple[str, str, str]]] = collections.defaultdict(set)
    examples: list[dict[str, Any]] = []
    peptide_with_mic: set[str] = set()
    total_activities = 0
    mic_activities = 0

    for record in records:
        peptide_id = str(record.get("id", "UNKNOWN"))
        activities = record.get("targetActivities") or []
        if not isinstance(activities, list):
            continue
        total_activities += len(activities)
        for activity in activities:
            if not isinstance(activity, dict) or not is_mic(activity):
                continue
            mic_activities += 1
            peptide_with_mic.add(peptide_id)
            raw = activity.get("concentration")
            kind = label_kind(raw)
            kind_counts[kind] += 1
            concentration_counts[str(raw).strip()] += 1

            unit = activity.get("unit")
            if isinstance(unit, dict):
                unit = unit.get("name") or unit.get("value")
            unit_counts[str(unit or "MISSING")] += 1

            leaves = list(flatten_leaves(activity))
            for path, _ in leaves:
                leaf_path_counts[path] += 1

            medium = first_matching_leaf(activity, ("medium",))
            cfu = first_matching_leaf(activity, ("cfu", "inoculum"))
            method = first_matching_leaf(activity, ("method", "assay"))
            reference = first_matching_leaf(activity, ("reference", "article", "publication", "doi"))

            for name, value in (("medium", medium), ("cfu", cfu), ("method", method), ("reference", reference)):
                if value is not None:
                    condition_presence[name] += 1
            if medium:
                medium_values[medium] += 1
            if cfu:
                cfu_values[cfu] += 1
            if method:
                method_values[method] += 1
            if reference:
                reference_values[reference] += 1

            species = species_name(activity)
            pair = (peptide_id, species)
            pair_counts[pair] += 1
            pair_contexts[pair].add((medium or "MISSING", cfu or "MISSING", method or "MISSING"))

            if len(examples) < 750:
                examples.append({
                    "peptide_id": peptide_id,
                    "species": species,
                    "concentration": raw,
                    "kind": kind,
                    "unit": unit,
                    "medium": medium,
                    "cfu": cfu,
                    "method": method,
                    "activity_keys": sorted(activity.keys()),
                })

    pair_sizes = list(pair_counts.values())
    context_sizes = [len(contexts) for contexts in pair_contexts.values()]
    censored_or_range = sum(
        kind_counts[name]
        for name in ("right_censored", "left_censored", "two_sided_censored", "finite_range", "plus_minus")
    )

    summary = {
        "audit_version": VERSION,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source": "DBAASP public REST endpoint https://dbaasp.org/peptides/{id}",
        "sampling": {
            "strategy": "deterministic equal-width stratified identifier sample",
            "requested_ids": requested,
            "successful_records": len(records),
            "not_found_404": missing_404,
            "other_errors": len(errors),
            "max_identifier": DEFAULT_MAX_ID,
        },
        "counts": {
            "all_target_activities": total_activities,
            "mic_measurements": mic_activities,
            "peptides_with_mic": len(peptide_with_mic),
            "peptide_species_pairs": len(pair_counts),
            "pairs_with_at_least_2_measurements": sum(v >= 2 for v in pair_sizes),
            "pairs_with_at_least_3_measurements": sum(v >= 3 for v in pair_sizes),
            "pairs_with_multiple_observed_contexts": sum(v >= 2 for v in context_sizes),
            "censored_range_or_plus_minus_measurements": censored_or_range,
        },
        "fractions": {
            "mic_among_all_target_activities": mic_activities / total_activities if total_activities else None,
            "censored_range_or_plus_minus_among_mic": censored_or_range / mic_activities if mic_activities else None,
            "condition_field_presence_among_mic": {
                key: value / mic_activities if mic_activities else None
                for key, value in condition_presence.items()
            },
            "repeated_pair_fraction": sum(v >= 2 for v in pair_sizes) / len(pair_sizes) if pair_sizes else None,
            "multiple_context_pair_fraction": sum(v >= 2 for v in context_sizes) / len(context_sizes) if context_sizes else None,
        },
        "label_kinds": serializable_counter(kind_counts),
        "units": serializable_counter(unit_counts, 30),
        "top_raw_concentrations": serializable_counter(concentration_counts, 40),
        "top_leaf_paths": serializable_counter(leaf_path_counts, 80),
        "top_medium_values": serializable_counter(medium_values, 30),
        "top_cfu_values": serializable_counter(cfu_values, 30),
        "top_method_values": serializable_counter(method_values, 30),
        "repeat_count_distribution": quantiles(pair_sizes),
        "context_count_distribution": quantiles(context_sizes),
        "errors_first_50": errors[:50],
    }
    return summary, examples


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-id", type=int, default=DEFAULT_MAX_ID)
    parser.add_argument("--sample-n", type=int, default=3000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out-dir", type=Path, default=Path("audit_output"))
    args = parser.parse_args()

    ids = stratified_ids(args.max_id, args.sample_n)
    records: list[dict[str, Any]] = []
    errors: list[str] = []
    missing_404 = 0
    started = time.perf_counter()

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(fetch_one, identifier): identifier for identifier in ids}
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            identifier, payload, error = future.result()
            if payload is not None:
                records.append(payload)
            elif error == "404":
                missing_404 += 1
            else:
                errors.append(f"{identifier}\t{error}")
            if index % 250 == 0:
                elapsed = time.perf_counter() - started
                print(f"progress={index}/{len(ids)} successful={len(records)} elapsed_s={elapsed:.1f}", flush=True)

    summary, examples = audit(records, len(ids), missing_404, errors)
    summary["sampling"]["elapsed_seconds"] = time.perf_counter() - started
    summary["sampling"]["max_identifier"] = args.max_id

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    with (args.out_dir / "examples.jsonl").open("w", encoding="utf-8") as handle:
        for item in examples:
            handle.write(json.dumps(item, ensure_ascii=False) + "\n")
    (args.out_dir / "errors.tsv").write_text("\n".join(errors), encoding="utf-8")
    readme = (
        "# Raw DBAASP MIC label audit\n\n"
        "This artifact is a deterministic viability audit. It does not train a model and does not claim novelty.\n\n"
        f"- audit version: {VERSION}\n"
        f"- requested IDs: {len(ids)}\n"
        f"- successful records: {len(records)}\n"
        f"- MIC measurements: {summary['counts']['mic_measurements']}\n"
        f"- censored/range/plus-minus fraction: "
        f"{summary['fractions']['censored_range_or_plus_minus_among_mic']}\n"
    )
    (args.out_dir / "README.md").write_text(readme, encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
