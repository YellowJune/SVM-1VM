#!/usr/bin/env python3
"""Audit native co-optimal alignment ambiguity in CRISPR-Bulge datasets."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import math
import pathlib
import statistics
import sys
import zipfile

EXPECTED_BYTES = 524_400_344
EXPECTED_SHA256 = "f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0"
REQUIRED_COLUMNS = {"Align.sgRNA", "Align.off-target", "sgRNA"}
GROUP_SUFFIXES = {
    "CHANGEseq": "/CHANGEseq/include_on_targets/CHANGEseq_CR_Lazzarotto_2020_dataset.csv",
    "FullGUIDEseq": "/FullGUIDEseq/include_on_targets/FullGUIDEseq_CR_Lazzarotto_2020_dataset.csv",
    "Refined_TrueOT": "/Refined_TrueOT.csv",
}


def file_sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def resolve_member(names: list[str], suffix: str) -> str:
    exact = [name for name in names if ("/" + name.lstrip("/")).endswith(suffix)]
    if len(exact) != 1:
        raise RuntimeError(
            f"expected exactly one archive member ending {suffix!r}; found {exact}"
        )
    return exact[0]


def clean_sequence(value: str) -> str:
    return "".join(base for base in str(value).upper().strip() if base in "ACGTN-")


def mismatch_count(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("aligned sequences must have equal length")
    mismatches = 0
    for a, b in zip(left, right):
        if a == "-" or b == "-":
            continue
        if a == "N" or b == "N":
            continue
        mismatches += a != b
    return mismatches


def enumerate_cooptimal(raw_guide: str, raw_target: str):
    if abs(len(raw_guide) - len(raw_target)) != 1:
        return []
    if len(raw_guide) < len(raw_target):
        shorter = raw_guide
        longer = raw_target
        guide_shorter = True
    else:
        shorter = raw_target
        longer = raw_guide
        guide_shorter = False

    protospacer_length = max(0, len(shorter) - 3)
    candidates = []
    for gap_position in range(protospacer_length + 1):
        gapped = shorter[:gap_position] + "-" + shorter[gap_position:]
        if guide_shorter:
            aligned_guide, aligned_target = gapped, longer
        else:
            aligned_guide, aligned_target = longer, gapped
        candidates.append(
            {
                "guide": aligned_guide,
                "target": aligned_target,
                "gap_position": gap_position,
                "mismatches": mismatch_count(aligned_guide, aligned_target),
            }
        )
    if not candidates:
        return []
    best = min(item["mismatches"] for item in candidates)
    unique = {}
    for item in candidates:
        if item["mismatches"] == best:
            unique[(item["guide"], item["target"])] = item
    return sorted(
        unique.values(),
        key=lambda item: (item["gap_position"], item["guide"], item["target"]),
    )


def float_value(value):
    try:
        return float(str(value).strip())
    except Exception:
        return None


def int_label(value):
    number = float_value(value)
    if number is None:
        return None
    return int(number > 0)


def classify_row(group: str, row: dict):
    reads = float_value(row.get("reads"))
    label = int_label(row.get("label"))
    if group == "Refined_TrueOT":
        if label is not None:
            return ("active" if label == 1 else "inactive"), "binary_label"
        if reads is not None:
            return ("active" if reads > 0 else "inactive"), "positive_signal"
        return "unusable", "missing_activity"
    if reads is None:
        return "unusable", "missing_reads"
    if reads >= 100:
        return "active", "reads_ge_100"
    if reads == 0:
        return "inactive", "reads_zero"
    return "low_signal_excluded", "reads_between_0_and_100"


def audit_group(archive: zipfile.ZipFile, member: str, group: str):
    counts = collections.Counter()
    ambiguity_counts = collections.Counter()
    active_ambiguous_guides = set()
    active_ambiguity_entropies = []
    all_ambiguity_entropies = []
    active_gap_spans = []
    examples = []
    column_names = []
    activity_rule_counts = collections.Counter()
    canonical_not_in_cooptimal = 0
    canonical_checked = 0

    with archive.open(member, "r") as raw:
        text = io.TextIOWrapper(raw, encoding="utf-8-sig", errors="replace", newline="")
        reader = csv.DictReader(text)
        column_names = reader.fieldnames or []
        missing = REQUIRED_COLUMNS - set(column_names)
        if missing:
            raise RuntimeError(f"{member} missing required columns: {sorted(missing)}")
        if "reads" not in column_names and "label" not in column_names:
            raise RuntimeError(f"{member} lacks both reads and label activity fields")

        for row in reader:
            counts["rows"] += 1
            state, rule = classify_row(group, row)
            counts[state] += 1
            activity_rule_counts[rule] += 1
            if state not in {"active", "inactive"}:
                continue

            aligned_guide = clean_sequence(row.get("Align.sgRNA", ""))
            aligned_target = clean_sequence(row.get("Align.off-target", ""))
            raw_guide = aligned_guide.replace("-", "")
            raw_target = aligned_target.replace("-", "")
            if not raw_guide or not raw_target:
                counts["missing_sequence"] += 1
                continue
            if abs(len(raw_guide) - len(raw_target)) != 1:
                counts["non_single_bulge"] += 1
                continue

            counts["bulge_sites"] += 1
            counts[f"{state}_bulge_sites"] += 1
            cooptimal = enumerate_cooptimal(raw_guide, raw_target)
            if len(cooptimal) < 2:
                continue

            ambiguity_counts["ambiguous_sites"] += 1
            ambiguity_counts[f"{state}_ambiguous_sites"] += 1
            entropy = math.log2(len(cooptimal))
            all_ambiguity_entropies.append(entropy)
            positions = [item["gap_position"] for item in cooptimal]
            span = max(positions) - min(positions)
            if state == "active":
                active_ambiguous_guides.add(
                    clean_sequence(row.get("sgRNA", "")) or raw_guide
                )
                active_ambiguity_entropies.append(entropy)
                active_gap_spans.append(span)

            canonical_pair = (aligned_guide, aligned_target)
            canonical_pairs = {
                (item["guide"], item["target"]) for item in cooptimal
            }
            if len(aligned_guide) == len(aligned_target):
                canonical_checked += 1
                if canonical_pair not in canonical_pairs:
                    canonical_not_in_cooptimal += 1

            if len(examples) < 100:
                examples.append(
                    {
                        "state": state,
                        "sgRNA": clean_sequence(row.get("sgRNA", "")) or raw_guide,
                        "raw_guide": raw_guide,
                        "raw_target": raw_target,
                        "canonical_guide": aligned_guide,
                        "canonical_target": aligned_target,
                        "cooptimal_count": len(cooptimal),
                        "gap_positions": positions,
                        "gap_span": span,
                        "best_mismatches": cooptimal[0]["mismatches"],
                        "reads": row.get("reads"),
                        "label": row.get("label"),
                    }
                )

    return {
        "group": group,
        "archive_member": member,
        "columns": column_names,
        "counts": dict(counts),
        "ambiguity_counts": dict(ambiguity_counts),
        "active_ambiguous_distinct_sgRNAs": len(active_ambiguous_guides),
        "active_ambiguous_sgRNAs": sorted(active_ambiguous_guides),
        "active_ambiguity_entropy_mean_bits": (
            statistics.mean(active_ambiguity_entropies)
            if active_ambiguity_entropies
            else 0.0
        ),
        "all_ambiguity_entropy_mean_bits": (
            statistics.mean(all_ambiguity_entropies)
            if all_ambiguity_entropies
            else 0.0
        ),
        "active_gap_span_ge_2_fraction": (
            sum(span >= 2 for span in active_gap_spans) / len(active_gap_spans)
            if active_gap_spans
            else 0.0
        ),
        "canonical_alignment_checks": canonical_checked,
        "canonical_not_in_cooptimal": canonical_not_in_cooptimal,
        "activity_rule_counts": dict(activity_rule_counts),
        "examples": examples,
    }


def main():
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True)
    parser.add_argument("--outdir", default="results/crispr_alignment_audit")
    args = parser.parse_args()
    archive_path = pathlib.Path(args.archive)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    archive_bytes = archive_path.stat().st_size
    archive_sha = file_sha256(archive_path)
    with zipfile.ZipFile(archive_path) as archive:
        infos = archive.infolist()
        names = [item.filename for item in infos if not item.is_dir()]
        manifest = [
            {
                "path": item.filename,
                "bytes": item.file_size,
                "compressed_bytes": item.compress_size,
                "crc32": f"{item.CRC:08x}",
            }
            for item in infos
            if not item.is_dir()
        ]
        resolved = {
            group: resolve_member(names, suffix)
            for group, suffix in GROUP_SUFFIXES.items()
        }
        groups = [
            audit_group(archive, resolved[group], group)
            for group in GROUP_SUFFIXES
        ]

    total_active = sum(group["counts"].get("active", 0) for group in groups)
    total_inactive = sum(group["counts"].get("inactive", 0) for group in groups)
    total_ambiguous = sum(
        group["ambiguity_counts"].get("ambiguous_sites", 0) for group in groups
    )
    active_ambiguous = sum(
        group["ambiguity_counts"].get("active_ambiguous_sites", 0)
        for group in groups
    )
    active_bulge = sum(
        group["counts"].get("active_bulge_sites", 0) for group in groups
    )
    active_guides = set()
    for group in groups:
        active_guides.update(group["active_ambiguous_sgRNAs"])
    contributing = [
        group["group"]
        for group in groups
        if group["ambiguity_counts"].get("active_ambiguous_sites", 0) >= 10
    ]
    span_numerator = sum(
        group["active_gap_span_ge_2_fraction"]
        * group["ambiguity_counts"].get("active_ambiguous_sites", 0)
        for group in groups
    )
    span_fraction = span_numerator / active_ambiguous if active_ambiguous else 0.0
    active_ambiguous_fraction = (
        active_ambiguous / active_bulge if active_bulge else 0.0
    )

    metrics = {
        "total_active_sites": total_active,
        "total_inactive_sites": total_inactive,
        "total_cooptimal_ambiguous_sites": total_ambiguous,
        "active_cooptimal_ambiguous_sites": active_ambiguous,
        "active_bulge_sites": active_bulge,
        "active_ambiguous_fraction_of_active_bulge": active_ambiguous_fraction,
        "active_ambiguous_distinct_sgRNAs": len(active_guides),
        "partitions_with_ge_10_active_ambiguous": contributing,
        "active_gap_span_ge_2_fraction": span_fraction,
    }
    gates = {
        "source_bytes_and_sha_match": (
            archive_bytes == EXPECTED_BYTES and archive_sha == EXPECTED_SHA256
        ),
        "three_core_partitions_present": len(groups) == 3,
        "active_ge_50000_and_inactive_ge_500000": (
            total_active >= 50_000 and total_inactive >= 500_000
        ),
        "total_ambiguous_ge_1000": total_ambiguous >= 1_000,
        "active_ambiguous_ge_100": active_ambiguous >= 100,
        "active_ambiguous_guides_ge_15": len(active_guides) >= 15,
        "two_partitions_ge_10_active_ambiguous": len(contributing) >= 2,
        "active_bulge_ambiguity_fraction_ge_0_05": active_ambiguous_fraction >= 0.05,
        "active_gap_span_ge_2_fraction_ge_0_25": span_fraction >= 0.25,
        "archive_manifest_recorded": len(manifest) > 0,
    }
    result = {
        "audit_version": "0.1.0",
        "preregistration_commit": "fb900aa69344ad0ea805a1021df44871c1f80a85",
        "source": {
            "url": "https://media.githubusercontent.com/media/OrensteinLab/CRISPR-Bulge/main/files/datasets.zip",
            "repository": "https://github.com/OrensteinLab/CRISPR-Bulge",
            "license": "MIT",
            "bytes": archive_bytes,
            "sha256": archive_sha,
            "expected_bytes": EXPECTED_BYTES,
            "expected_sha256": EXPECTED_SHA256,
        },
        "resolved_members": resolved,
        "archive_file_count": len(manifest),
        "groups": groups,
        "metrics": metrics,
        "gates": gates,
        "numeric_gate_pass": all(gates.values()),
        "literature_gate": "manual_pending",
        "decision": "KEEP_FOR_METHOD_PROOF" if all(gates.values()) else "REJECT",
    }
    (outdir / "audit.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    (outdir / "archive_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    lines = [
        "# CRISPR co-optimal alignment feasibility audit",
        "",
        f"- Archive bytes: {archive_bytes:,}",
        f"- Archive SHA-256: `{archive_sha}`",
        f"- Core partitions: {len(groups)}",
        f"- Active / inactive sites: {total_active:,} / {total_inactive:,}",
        f"- Total co-optimal ambiguous sites: {total_ambiguous:,}",
        f"- Active co-optimal ambiguous sites: {active_ambiguous:,}",
        f"- Active bulge sites: {active_bulge:,}",
        f"- Ambiguous fraction among active bulge sites: {active_ambiguous_fraction:.4%}",
        f"- Distinct sgRNAs with active ambiguity: {len(active_guides):,}",
        f"- Partitions with >=10 active ambiguous sites: {', '.join(contributing) or 'none'}",
        f"- Active ambiguous sites with gap-position span >=2: {span_fraction:.4%}",
        "",
        "## Frozen numerical gates",
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
            "The separate no-direct-prior literature gate remains manual.",
        ]
    )
    (outdir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
