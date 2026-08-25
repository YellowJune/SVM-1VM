#!/usr/bin/env python3
"""Phase-A audit: prediction instability under co-optimal CRISPR alignments."""

from __future__ import annotations

import argparse
import collections
import csv
import gc
import hashlib
import heapq
import io
import json
import math
import os
import pathlib
import random
import statistics
import sys
import time
import zipfile

import numpy as np
import pandas as pd

EXPECTED_BYTES = 524_400_344
EXPECTED_SHA256 = "f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0"
SOURCE_REPO_COMMIT = "3eddcd5bfcaff00b2bdf29425116ec9756ace870"
GROUP_SUFFIXES = {
    "CHANGEseq": "/CHANGEseq/include_on_targets/CHANGEseq_CR_Lazzarotto_2020_dataset.csv",
    "FullGUIDEseq": "/FullGUIDEseq/include_on_targets/FullGUIDEseq_CR_Lazzarotto_2020_dataset.csv",
    "Refined_TrueOT": "/Refined_TrueOT.csv",
}
MODEL_TEMPLATE = (
    "files/bulges/1_folds/5_revision_ensemble_{}_exclude_RHAMPseq_continue_from_change_seq/"
    "read_ts_0/cleavage_models/aligned/FullGUIDEseq/classification/c_2/"
    "ln_x_plus_one_trans/model_fold_0"
)


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_sequence(value) -> str:
    return "".join(base for base in str(value).upper().strip() if base in "ACGTN-")


def float_value(value):
    try:
        return float(str(value).strip())
    except Exception:
        return None


def int_label(value):
    number = float_value(value)
    return None if number is None else int(number > 0)


def state_for(group: str, row: dict):
    reads = float_value(row.get("reads"))
    label = int_label(row.get("label"))
    if group == "Refined_TrueOT":
        if label is not None:
            return "active" if label else "inactive"
        if reads is not None:
            return "active" if reads > 0 else "inactive"
        return "unusable"
    if reads is None:
        return "unusable"
    if reads >= 100:
        return "active"
    if reads == 0:
        return "inactive"
    return "excluded"


def mismatch_count(left: str, right: str) -> int:
    result = 0
    for a, b in zip(left, right):
        if a == "-" or b == "-" or a == "N" or b == "N":
            continue
        result += a != b
    return result


def enumerate_cooptimal(raw_guide: str, raw_target: str):
    if abs(len(raw_guide) - len(raw_target)) != 1:
        return []
    if len(raw_guide) < len(raw_target):
        shorter, longer, guide_shorter = raw_guide, raw_target, True
    else:
        shorter, longer, guide_shorter = raw_target, raw_guide, False
    candidates = []
    # The last three bases are PAM. A gap may occur at spacer boundaries 0..L.
    protospacer_length = max(0, len(shorter) - 3)
    for gap_position in range(protospacer_length + 1):
        gapped = shorter[:gap_position] + "-" + shorter[gap_position:]
        guide, target = (gapped, longer) if guide_shorter else (longer, gapped)
        candidates.append({
            "guide": guide,
            "target": target,
            "gap_position": gap_position,
            "mismatches": mismatch_count(guide, target),
        })
    best = min(item["mismatches"] for item in candidates)
    unique = {
        (item["guide"], item["target"]): item
        for item in candidates if item["mismatches"] == best
    }
    return sorted(
        unique.values(),
        key=lambda item: (item["gap_position"], item["guide"], item["target"]),
    )


def resolve_member(names, suffix):
    matches = [name for name in names if ("/" + name.lstrip("/")).endswith(suffix)]
    if len(matches) != 1:
        raise RuntimeError(f"expected one member ending {suffix!r}, found {matches}")
    return matches[0]


def stable_site_id(group: str, row_number: int, row: dict, raw_guide: str, raw_target: str):
    location_fields = [
        row.get("chrom", ""),
        row.get("Align.chromStart", ""),
        row.get("Align.chromEnd", ""),
        row.get("strand", ""),
    ]
    payload = "|".join(
        [group, str(row_number), clean_sequence(row.get("sgRNA", "")), raw_guide, raw_target]
        + [str(item) for item in location_fields]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def compact_record(group, state, site_id, row, raw_guide, raw_target,
                   canonical_guide, canonical_target, candidates):
    return {
        "site_id": site_id,
        "partition": group,
        "state": state,
        "sgRNA": clean_sequence(row.get("sgRNA", "")) or raw_guide,
        "raw_guide": raw_guide,
        "raw_target": raw_target,
        "canonical_guide": canonical_guide,
        "canonical_target": canonical_target,
        "reads": row.get("reads", ""),
        "label": row.get("label", ""),
        "candidates": candidates,
    }


def collect_sites(archive_path: pathlib.Path, inactive_limit: int):
    active = []
    inactive_heap = []
    source_counts = collections.Counter()
    group_counts = collections.defaultdict(collections.Counter)
    with zipfile.ZipFile(archive_path) as archive:
        names = [item.filename for item in archive.infolist() if not item.is_dir()]
        resolved = {g: resolve_member(names, s) for g, s in GROUP_SUFFIXES.items()}
        for group, member in resolved.items():
            with archive.open(member, "r") as raw:
                text = io.TextIOWrapper(
                    raw, encoding="utf-8-sig", errors="replace", newline=""
                )
                reader = csv.DictReader(text)
                required = {"Align.sgRNA", "Align.off-target", "sgRNA"}
                missing = required - set(reader.fieldnames or [])
                if missing:
                    raise RuntimeError(f"{member} missing columns {sorted(missing)}")
                for row_number, row in enumerate(reader):
                    state = state_for(group, row)
                    source_counts[f"{group}:{state}"] += 1
                    if state not in {"active", "inactive"}:
                        continue
                    canonical_guide = clean_sequence(row.get("Align.sgRNA", ""))
                    canonical_target = clean_sequence(row.get("Align.off-target", ""))
                    raw_guide = canonical_guide.replace("-", "")
                    raw_target = canonical_target.replace("-", "")
                    if not raw_guide or not raw_target:
                        continue
                    candidates = enumerate_cooptimal(raw_guide, raw_target)
                    if len(candidates) < 2:
                        continue
                    group_counts[group][f"{state}_ambiguous"] += 1
                    site_id = stable_site_id(
                        group, row_number, row, raw_guide, raw_target
                    )
                    record = compact_record(
                        group, state, site_id, row, raw_guide, raw_target,
                        canonical_guide, canonical_target, candidates
                    )
                    if state == "active":
                        active.append(record)
                    else:
                        rank = int.from_bytes(
                            hashlib.sha256(
                                ("inactive-v1|" + site_id).encode("utf-8")
                            ).digest(),
                            "big",
                        )
                        item = (-rank, site_id, record)
                        if len(inactive_heap) < inactive_limit:
                            heapq.heappush(inactive_heap, item)
                        elif item[0] > inactive_heap[0][0]:
                            heapq.heapreplace(inactive_heap, item)
    inactive = [item[2] for item in inactive_heap]
    inactive.sort(key=lambda item: item["site_id"])
    active.sort(key=lambda item: item["site_id"])
    return active, inactive, dict(source_counts), {
        key: dict(value) for key, value in group_counts.items()
    }


def build_long_frame(sites):
    rows = []
    for site in sites:
        canonical_pair = (site["canonical_guide"], site["canonical_target"])
        candidate_pairs = {
            (item["guide"], item["target"]) for item in site["candidates"]
        }
        for candidate_index, item in enumerate(site["candidates"]):
            rows.append({
                "site_id": site["site_id"],
                "partition": site["partition"],
                "state": site["state"],
                "sgRNA": site["sgRNA"],
                "Align.sgRNA": item["guide"],
                "Align.off-target": item["target"],
                "candidate_index": candidate_index,
                "gap_position": item["gap_position"],
                "mismatches": item["mismatches"],
                "is_cooptimal": 1,
                "is_canonical": int(
                    (item["guide"], item["target"]) == canonical_pair
                ),
            })
        if canonical_pair not in candidate_pairs and len(canonical_pair[0]) == len(canonical_pair[1]):
            rows.append({
                "site_id": site["site_id"],
                "partition": site["partition"],
                "state": site["state"],
                "sgRNA": site["sgRNA"],
                "Align.sgRNA": canonical_pair[0],
                "Align.off-target": canonical_pair[1],
                "candidate_index": -1,
                "gap_position": -1,
                "mismatches": mismatch_count(*canonical_pair),
                "is_cooptimal": 0,
                "is_canonical": 1,
            })
    return pd.DataFrame(rows)


def model_files(repo: pathlib.Path):
    prefixes = [repo / MODEL_TEMPLATE.format(i) for i in range(5)]
    files = []
    for prefix in prefixes:
        for suffix in (".h5", ".pkl", "_model_parameters.json"):
            path = pathlib.Path(str(prefix) + suffix)
            if not path.is_file():
                raise FileNotFoundError(path)
            files.append(path)
    return prefixes, files


def score_frame(frame: pd.DataFrame, repo: pathlib.Path):
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    sys.path.insert(0, str(repo))
    old_cwd = pathlib.Path.cwd()
    os.chdir(repo)
    try:
        from OT_deep_score_src.data_processing_utilities import build_sequence_features
        from OT_deep_score_src.general_utilities import Encoding_type, Padding_type
        from OT_deep_score_src.models_inter import Model

        input_frame = frame[["Align.sgRNA", "Align.off-target"]].copy()
        features = build_sequence_features(
            input_frame,
            include_distance_feature=False,
            include_sequence_features=True,
            include_gmt_score=False,
            include_nuclea_seq_score=False,
            bulges=True,
            padding_type=Padding_type.GAP,
            aligned=True,
            encoding_type=Encoding_type.ONE_HOT,
            flat_encoding=False,
        )
        prefixes, _ = model_files(repo)
        prediction_columns = []
        load_and_predict_seconds = []
        for index, prefix in enumerate(prefixes):
            start = time.perf_counter()
            model = Model.load_model_instance(str(prefix))
            prediction = np.asarray(model.predict(features)).reshape(-1)
            if len(prediction) != len(frame) or not np.isfinite(prediction).all():
                raise RuntimeError(f"invalid predictions from ensemble component {index}")
            column = f"pred_component_{index}"
            frame[column] = prediction
            prediction_columns.append(column)
            load_and_predict_seconds.append(time.perf_counter() - start)
            del model
            gc.collect()
        frame["pred_ensemble_mean"] = frame[prediction_columns].mean(axis=1)
        return frame, prediction_columns, load_and_predict_seconds
    finally:
        os.chdir(old_cwd)


def site_summaries(sites, scored: pd.DataFrame, prediction_columns):
    by_site = {key: group for key, group in scored.groupby("site_id", sort=False)}
    summaries = []
    active_alignment_rows = []
    canonical_members = 0
    canonical_checks = 0
    for site in sites:
        group = by_site[site["site_id"]]
        candidates = group[group["is_cooptimal"] == 1]
        canonical = group[group["is_canonical"] == 1]
        canonical_checks += 1
        canonical_in = bool(
            ((group["is_cooptimal"] == 1) & (group["is_canonical"] == 1)).any()
        )
        canonical_members += canonical_in
        values = candidates["pred_ensemble_mean"].to_numpy(dtype=float)
        row = {
            "site_id": site["site_id"],
            "partition": site["partition"],
            "state": site["state"],
            "sgRNA": site["sgRNA"],
            "raw_guide": site["raw_guide"],
            "raw_target": site["raw_target"],
            "cooptimal_count": len(candidates),
            "gap_positions": ";".join(str(x) for x in candidates["gap_position"]),
            "gap_span": int(candidates["gap_position"].max() - candidates["gap_position"].min()),
            "canonical_in_cooptimal": int(canonical_in),
            "prediction_min": float(values.min()),
            "prediction_max": float(values.max()),
            "prediction_range": float(values.max() - values.min()),
            "prediction_std": float(values.std(ddof=0)),
            "canonical_prediction": (
                float(canonical.iloc[0]["pred_ensemble_mean"])
                if len(canonical) else math.nan
            ),
            "reads": site["reads"],
            "label": site["label"],
        }
        for column in prediction_columns:
            component_values = candidates[column].to_numpy(dtype=float)
            row[f"{column}_range"] = float(
                component_values.max() - component_values.min()
            )
        summaries.append(row)
        if site["state"] == "active":
            active_alignment_rows.append(group)
    return (
        pd.DataFrame(summaries),
        pd.concat(active_alignment_rows, ignore_index=True),
        canonical_members,
        canonical_checks,
    )


def ordering_flip_fraction(active_summary: pd.DataFrame, pair_count: int = 100_000):
    if len(active_summary) < 2:
        return 0.0, 0
    intervals = active_summary[["prediction_min", "prediction_max"]].to_numpy(float)
    rng = random.Random(20260825)
    flips = 0
    used = 0
    for _ in range(pair_count):
        i = rng.randrange(len(intervals))
        j = rng.randrange(len(intervals) - 1)
        if j >= i:
            j += 1
        i_min, i_max = intervals[i]
        j_min, j_max = intervals[j]
        can_reverse = (i_min < j_max) and (i_max > j_min)
        flips += can_reverse
        used += 1
    return flips / used if used else 0.0, used


def write_report(outdir, metrics, gates):
    lines = [
        "# Published-model co-optimal alignment instability audit",
        "",
        f"- Audited active ambiguous sites: {metrics['active_sites']:,}",
        f"- Deterministically sampled inactive ambiguous sites: {metrics['inactive_sites']:,}",
        f"- Scored co-optimal alignment rows: {metrics['cooptimal_alignment_rows']:,}",
        f"- Canonical alignment membership: {metrics['canonical_membership_fraction']:.4%}",
        f"- Active sites with prediction range >=0.05: {metrics['active_range_ge_0_05_fraction']:.4%}",
        f"- Median active prediction range: {metrics['active_prediction_range_median']:.6f}",
        f"- Pairwise ordering-reversal opportunity: {metrics['pair_order_flip_fraction']:.4%}",
        f"- Partitions meeting instability definition: {', '.join(metrics['unstable_partitions']) or 'none'}",
        "",
        "## Frozen Phase-A gates",
        "",
    ]
    lines.extend(f"- [{'x' if value else ' '}] {key}" for key, value in gates.items())
    lines.extend([
        "",
        f"Decision: **{'PASS_TO_PHASE_B' if all(gates.values()) else 'REJECT'}**",
        "",
        "Instability for a partition is defined before scoring as both n>=10 active ambiguous sites and >=10% with prediction range >=0.05.",
    ])
    (outdir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True)
    parser.add_argument("--source-repo", required=True)
    parser.add_argument("--outdir", default="results/crispr_alignment_phase_a")
    parser.add_argument("--inactive-limit", type=int, default=100_000)
    args = parser.parse_args()

    archive = pathlib.Path(args.archive).resolve()
    repo = pathlib.Path(args.source_repo).resolve()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    archive_bytes = archive.stat().st_size
    archive_sha = sha256_file(archive)
    if archive_bytes != EXPECTED_BYTES or archive_sha != EXPECTED_SHA256:
        raise RuntimeError("pinned source archive bytes/SHA-256 mismatch")

    prefixes, required_model_files = model_files(repo)
    model_hashes = {
        str(path.relative_to(repo)): {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in required_model_files
    }

    collect_start = time.perf_counter()
    active, inactive, source_counts, group_counts = collect_sites(
        archive, args.inactive_limit
    )
    selected = active + inactive
    collect_seconds = time.perf_counter() - collect_start
    if len(inactive) != args.inactive_limit:
        raise RuntimeError(
            f"expected {args.inactive_limit} inactive samples, obtained {len(inactive)}"
        )

    long_frame = build_long_frame(selected)
    score_start = time.perf_counter()
    scored, prediction_columns, model_seconds = score_frame(long_frame, repo)
    score_seconds = time.perf_counter() - score_start

    summary, active_predictions, canonical_members, canonical_checks = site_summaries(
        selected, scored, prediction_columns
    )
    active_summary = summary[summary["state"] == "active"].copy()
    inactive_summary = summary[summary["state"] == "inactive"].copy()
    flip_fraction, pair_count = ordering_flip_fraction(active_summary)

    partition_metrics = {}
    unstable_partitions = []
    for partition, group in active_summary.groupby("partition"):
        fraction = float((group["prediction_range"] >= 0.05).mean())
        item = {
            "active_ambiguous_sites": int(len(group)),
            "range_ge_0_05_fraction": fraction,
            "prediction_range_median": float(group["prediction_range"].median()),
            "unstable_definition_pass": bool(len(group) >= 10 and fraction >= 0.10),
        }
        partition_metrics[partition] = item
        if item["unstable_definition_pass"]:
            unstable_partitions.append(partition)

    active_range_fraction = float(
        (active_summary["prediction_range"] >= 0.05).mean()
    )
    median_range = float(active_summary["prediction_range"].median())
    canonical_fraction = canonical_members / canonical_checks if canonical_checks else 0.0
    finite_predictions = bool(
        scored[prediction_columns + ["pred_ensemble_mean"]]
        .apply(np.isfinite)
        .all()
        .all()
    )

    metrics = {
        "audit_version": "0.1.0",
        "preregistration_commit": "b7fb312d253605509c6cf4ba0cb00a55e5175a33",
        "source_repo_commit": SOURCE_REPO_COMMIT,
        "archive": {
            "bytes": archive_bytes,
            "sha256": archive_sha,
        },
        "active_sites": int(len(active_summary)),
        "inactive_sites": int(len(inactive_summary)),
        "cooptimal_alignment_rows": int(scored["is_cooptimal"].sum()),
        "all_scored_rows_including_noncooptimal_canonical": int(len(scored)),
        "canonical_members": int(canonical_members),
        "canonical_checks": int(canonical_checks),
        "canonical_membership_fraction": canonical_fraction,
        "active_range_ge_0_05_fraction": active_range_fraction,
        "active_prediction_range_median": median_range,
        "active_prediction_range_mean": float(active_summary["prediction_range"].mean()),
        "active_prediction_range_max": float(active_summary["prediction_range"].max()),
        "pair_order_flip_fraction": float(flip_fraction),
        "pair_order_comparisons": int(pair_count),
        "partition_metrics": partition_metrics,
        "unstable_partitions": sorted(unstable_partitions),
        "source_state_counts": source_counts,
        "ambiguity_group_counts": group_counts,
        "collection_seconds": collect_seconds,
        "scoring_seconds": score_seconds,
        "component_load_and_predict_seconds": model_seconds,
        "model_hashes": model_hashes,
    }
    gates = {
        "source_bytes_and_sha_match": (
            archive_bytes == EXPECTED_BYTES and archive_sha == EXPECTED_SHA256
        ),
        "five_models_finite": len(prediction_columns) == 5 and finite_predictions,
        "canonical_membership_ge_0_95": canonical_fraction >= 0.95,
        "active_range_ge_0_05_fraction_ge_0_10": active_range_fraction >= 0.10,
        "active_prediction_range_median_ge_0_01": median_range >= 0.01,
        "pair_order_flip_fraction_ge_0_05": flip_fraction >= 0.05,
        "instability_in_ge_2_partitions": len(unstable_partitions) >= 2,
    }
    result = {
        "metrics": metrics,
        "gates": gates,
        "decision": "PASS_TO_PHASE_B" if all(gates.values()) else "REJECT",
    }

    summary.to_csv(outdir / "site_summary.csv.gz", index=False, compression="gzip")
    active_predictions.to_csv(
        outdir / "active_alignment_predictions.csv.gz",
        index=False,
        compression="gzip",
    )
    (outdir / "metrics.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    (outdir / "model_manifest.json").write_text(
        json.dumps(model_hashes, indent=2), encoding="utf-8"
    )
    write_report(outdir, metrics, gates)


if __name__ == "__main__":
    main()
