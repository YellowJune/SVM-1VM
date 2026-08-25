#!/usr/bin/env python3
"""Audit graph-preserving page-label non-identifiability in bpRNA-1m."""

from __future__ import annotations

import hashlib
import json
import os
import random
import statistics
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

API = "https://datasets-server.huggingface.co/parquet?dataset=multimolecule/bprna"
OUT = Path("rna_page_gauge_audit/audit")
CACHE = OUT / "cache"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

PAIRS = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
PAIRS.extend((chr(ord("A") + i), chr(ord("a") + i)) for i in range(26))
OPEN_TO_PAGE = {a: i for i, (a, b) in enumerate(PAIRS)}
CLOSE_TO_PAGE = {b: i for i, (a, b) in enumerate(PAIRS)}


def fetch_bytes(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research-audit/0.2"})
    with urllib.request.urlopen(req, timeout=180) as response:
        return response.read()


def parse_structure(structure: str):
    stacks = [[] for _ in PAIRS]
    base_pairs = []
    used = set()
    unknown = []
    for idx, char in enumerate(structure):
        if char == ".":
            continue
        if char in OPEN_TO_PAGE:
            page = OPEN_TO_PAGE[char]
            stacks[page].append(idx)
            used.add(page)
        elif char in CLOSE_TO_PAGE:
            page = CLOSE_TO_PAGE[char]
            used.add(page)
            if not stacks[page]:
                return None, None, [f"unmatched-close:{idx}:{char}"]
            base_pairs.append((stacks[page].pop(), idx))
        else:
            unknown.append(f"unknown:{idx}:{char}")
    for page, stack in enumerate(stacks):
        if stack:
            unknown.append(f"unmatched-open:{page}:{len(stack)}")
    if unknown:
        return None, None, unknown
    return frozenset(base_pairs), tuple(sorted(used)), []


def permute_pages(structure: str, pages: tuple[int, ...], rng: random.Random) -> str:
    shuffled = list(pages)
    rng.shuffle(shuffled)
    if len(shuffled) > 1 and all(a == b for a, b in zip(pages, shuffled)):
        shuffled = shuffled[1:] + shuffled[:1]
    mapping = dict(zip(pages, shuffled))
    out = []
    for char in structure:
        if char in OPEN_TO_PAGE:
            out.append(PAIRS[mapping[OPEN_TO_PAGE[char]]][0])
        elif char in CLOSE_TO_PAGE:
            out.append(PAIRS[mapping[CLOSE_TO_PAGE[char]]][1])
        else:
            out.append(char)
    return "".join(out)


api_raw = fetch_bytes(API)
api_sha = hashlib.sha256(api_raw).hexdigest()
api_obj = json.loads(api_raw)
files = api_obj.get("parquet_files", [])
if not files:
    raise RuntimeError(f"No parquet files returned: {api_obj}")

selected = []
for item in files:
    if item.get("split") != "train":
        continue
    url = item["url"]
    name = Path(urllib.parse.urlparse(url).path).name or item.get("filename", "data.parquet")
    target = CACHE / name
    raw = fetch_bytes(url)
    target.write_bytes(raw)
    selected.append({
        "config": item.get("config"),
        "split": item.get("split"),
        "url": url,
        "filename": name,
        "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "path": str(target),
    })

if not selected:
    raise RuntimeError(f"No train parquet selected from {len(files)} files")

records = []
schemas = []
for meta in selected:
    table = pq.read_table(meta["path"])
    schemas.append(str(table.schema))
    names = set(table.column_names)
    required = {"id", "sequence", "secondary_structure"}
    if not required.issubset(names):
        raise RuntimeError(f"Missing columns {required - names}; found {table.column_names}")
    records.extend(zip(
        table["id"].to_pylist(),
        table["sequence"].to_pylist(),
        table["secondary_structure"].to_pylist(),
    ))

rng = random.Random(20260825)
page_counts = Counter()
source_counts = Counter()
lengths = []
pk_lengths = []
parse_errors = []
length_mismatch = 0
pk_records = 0
graph_preserved = 0
permutation_trials = 0
hamming_all = []
hamming_paired = []
pair_densities = []
page2_pair_shares = []

for rec_id, sequence, structure in records:
    sequence = str(sequence)
    structure = str(structure)
    lengths.append(len(sequence))
    source = str(rec_id).split("_")[1] if "_" in str(rec_id) else "unknown"
    source_counts[source] += 1
    if len(sequence) != len(structure):
        length_mismatch += 1
        continue
    base_pairs, pages, errors = parse_structure(structure)
    if errors:
        if len(parse_errors) < 20:
            parse_errors.append({"id": rec_id, "errors": errors[:5]})
        continue
    page_count = len(pages)
    page_counts[page_count] += 1
    pair_density = 2 * len(base_pairs) / max(1, len(structure))
    pair_densities.append(pair_density)
    if page_count >= 2:
        pk_records += 1
        pk_lengths.append(len(sequence))
        page2_pairs = sum(structure.count(PAIRS[p][0]) for p in pages[1:])
        page2_pair_shares.append(page2_pairs / max(1, len(base_pairs)))
        for _ in range(3):
            changed = permute_pages(structure, pages, rng)
            changed_pairs, changed_pages, changed_errors = parse_structure(changed)
            permutation_trials += 1
            if not changed_errors and changed_pairs == base_pairs:
                graph_preserved += 1
            diff = sum(a != b for a, b in zip(structure, changed))
            paired_positions = 2 * len(base_pairs)
            hamming_all.append(diff / max(1, len(structure)))
            hamming_paired.append(diff / max(1, paired_positions))

valid_records = sum(page_counts.values())
le2 = sum(n for k, n in page_counts.items() if k <= 2)
gate = {
    "minimum_pseudoknotted_records_5000": pk_records >= 5000,
    "at_least_95pct_valid_records_page_le_2": (le2 / valid_records >= 0.95) if valid_records else False,
    "mean_graph_preserving_full_token_hamming_at_least_0_20": (
        statistics.mean(hamming_all) >= 0.20 if hamming_all else False
    ),
}
summary = {
    "audit_version": "0.2.0",
    "source_api": API,
    "source_api_bytes": len(api_raw),
    "source_api_sha256": api_sha,
    "parquet_files": selected,
    "schemas": schemas,
    "records": len(records),
    "valid_records": valid_records,
    "length_mismatch": length_mismatch,
    "parse_error_records": len(records) - valid_records - length_mismatch,
    "parse_error_examples": parse_errors,
    "source_counts": dict(source_counts),
    "page_count_distribution": {str(k): v for k, v in sorted(page_counts.items())},
    "pseudoknotted_records": pk_records,
    "pseudoknotted_fraction": pk_records / valid_records if valid_records else None,
    "page_le_2_fraction": le2 / valid_records if valid_records else None,
    "permutation_trials": permutation_trials,
    "graph_preservation_fraction": graph_preserved / permutation_trials if permutation_trials else None,
    "mean_full_token_hamming_under_graph_preserving_page_permutation": (
        statistics.mean(hamming_all) if hamming_all else None
    ),
    "mean_paired_token_hamming_under_graph_preserving_page_permutation": (
        statistics.mean(hamming_paired) if hamming_paired else None
    ),
    "mean_pair_density": statistics.mean(pair_densities) if pair_densities else None,
    "mean_nonprimary_page_pair_share_in_pk_records": (
        statistics.mean(page2_pair_shares) if page2_pair_shares else None
    ),
    "length": {
        "min": min(lengths) if lengths else None,
        "median": statistics.median(lengths) if lengths else None,
        "max": max(lengths) if lengths else None,
        "pk_median": statistics.median(pk_lengths) if pk_lengths else None,
        "pk_max": max(pk_lengths) if pk_lengths else None,
    },
    "preregistered_gate": gate,
    "gate_pass": all(gate.values()),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
with (OUT / "manifest.sha256").open("w", encoding="utf-8") as handle:
    for path in sorted(OUT.glob("*")):
        if path.is_file() and path.name != "manifest.sha256":
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n")
print(json.dumps(summary, indent=2))
