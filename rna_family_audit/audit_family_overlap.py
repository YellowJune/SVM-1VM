#!/usr/bin/env python3
"""Audit whether ncRNA class benchmarks test unseen sequences or merely seen Rfam families."""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import re
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Iterable, Iterator

VERSION = "0.1.0"
DEFAULT_URL = "https://zenodo.org/records/17358094/files/seq_cls_data.zip?download=1"
RFAM_RE = re.compile(r"RF\d{5}", re.IGNORECASE)
FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas")

Record = tuple[str | None, str, str, str]


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "rna-family-audit/0.1 research"})
    with urllib.request.urlopen(request, timeout=120) as response, destination.open("wb") as handle:
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            handle.write(block)


def normalize_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper().replace("U", "T")


def reverse_complement(sequence: str) -> str:
    table = str.maketrans(
        "ACGTRYSWKMBDHVN",
        "TGCAYRSWMKVHDBN",
    )
    return sequence.translate(table)[::-1]


def sequence_digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii", "ignore")).hexdigest()


def parse_label(header: str) -> str:
    parts = header.strip().split()
    return parts[-1] if len(parts) >= 2 else "UNKNOWN"


def parse_family(header: str) -> str | None:
    match = RFAM_RE.search(header)
    return match.group(0).upper() if match else None


def read_fasta(lines: Iterable[bytes]) -> Iterator[Record]:
    header: str | None = None
    chunks: list[str] = []
    for raw in lines:
        line = raw.decode("utf-8", "replace").strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                sequence = normalize_sequence("".join(chunks))
                yield parse_family(header), parse_label(header), sequence, header
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        sequence = normalize_sequence("".join(chunks))
        yield parse_family(header), parse_label(header), sequence, header


def role_for(member: str) -> str | None:
    name = Path(member).name.lower()
    if "train" in name:
        return "train"
    if "test" in name:
        return "test"
    return None


def group_key(member: str) -> str:
    path = Path(member)
    stem = path.stem.lower()
    stem = re.sub(r"(?:^|[_\-.])(train|test)(?:$|[_\-.])", "_split_", stem)
    return str(path.parent / stem)


def majority(counter: collections.Counter[str]) -> str:
    return sorted(counter.items(), key=lambda item: (-item[1], item[0]))[0][0]


def summarize_pair(train: list[Record], test: list[Record], key: str) -> dict:
    train_family_labels: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    union_family_labels: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    class_families: dict[str, set[str]] = collections.defaultdict(set)
    class_family_sequence_counts: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)

    train_seq_labels: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    train_rc_labels: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)

    for family, label, sequence, _ in train:
        digest = sequence_digest(sequence)
        train_seq_labels[digest][label] += 1
        train_rc_labels[sequence_digest(reverse_complement(sequence))][label] += 1
        if family:
            train_family_labels[family][label] += 1

    for family, label, _, _ in train + test:
        if family:
            union_family_labels[family][label] += 1
            class_families[label].add(family)
            class_family_sequence_counts[label][family] += 1

    test_with_family = 0
    test_seen_family = 0
    family_correct = 0
    family_label_conflicts = 0
    exact_overlap = 0
    exact_correct = 0
    rc_overlap = 0
    rc_correct = 0

    examples_unseen: list[dict[str, str]] = []
    for family, label, sequence, header in test:
        digest = sequence_digest(sequence)
        if digest in train_seq_labels:
            exact_overlap += 1
            exact_correct += int(majority(train_seq_labels[digest]) == label)
        rc_digest = sequence_digest(reverse_complement(sequence))
        if rc_digest in train_seq_labels or digest in train_rc_labels:
            rc_overlap += 1
            candidates = train_seq_labels.get(rc_digest, collections.Counter()) + train_rc_labels.get(digest, collections.Counter())
            if candidates:
                rc_correct += int(majority(candidates) == label)

        if family:
            test_with_family += 1
            if family in train_family_labels:
                test_seen_family += 1
                predicted = majority(train_family_labels[family])
                family_correct += int(predicted == label)
                family_label_conflicts += int(len(train_family_labels[family]) > 1)
            elif len(examples_unseen) < 30:
                examples_unseen.append({"family": family, "label": label, "header": header})

    class_details = {}
    feasible_classes = 0
    for label in sorted(class_families):
        fam_counts = class_family_sequence_counts[label]
        families_ge5 = sum(count >= 5 for count in fam_counts.values())
        families_ge20 = sum(count >= 20 for count in fam_counts.values())
        total_sequences = sum(fam_counts.values())
        if families_ge5 >= 5 and total_sequences >= 100:
            feasible_classes += 1
        class_details[label] = {
            "families": len(fam_counts),
            "families_with_at_least_5_sequences": families_ge5,
            "families_with_at_least_20_sequences": families_ge20,
            "sequences_with_family_id": total_sequences,
            "largest_family_sequences": max(fam_counts.values(), default=0),
        }

    family_consistency = {
        family: dict(labels)
        for family, labels in union_family_labels.items()
        if len(labels) > 1
    }

    return {
        "dataset_key": key,
        "train_sequences": len(train),
        "test_sequences": len(test),
        "train_sequences_with_family_id": sum(record[0] is not None for record in train),
        "test_sequences_with_family_id": test_with_family,
        "train_families": len(train_family_labels),
        "test_families": len({record[0] for record in test if record[0]}),
        "shared_families": len({record[0] for record in test if record[0]} & set(train_family_labels)),
        "test_seen_family_sequences": test_seen_family,
        "test_seen_family_fraction_among_family_annotated": (
            test_seen_family / test_with_family if test_with_family else None
        ),
        "family_lookup_oracle_accuracy_among_seen": (
            family_correct / test_seen_family if test_seen_family else None
        ),
        "seen_family_rows_with_ambiguous_train_mapping": family_label_conflicts,
        "exact_sequence_overlap": exact_overlap,
        "exact_sequence_overlap_fraction": exact_overlap / len(test) if test else None,
        "exact_lookup_accuracy": exact_correct / exact_overlap if exact_overlap else None,
        "reverse_complement_overlap": rc_overlap,
        "reverse_complement_lookup_accuracy": rc_correct / rc_overlap if rc_overlap else None,
        "classes_feasible_for_family_holdout": feasible_classes,
        "class_feasibility_rule": ">=5 families with >=5 sequences and >=100 family-annotated sequences",
        "class_details": class_details,
        "families_assigned_to_multiple_classes": family_consistency,
        "unseen_family_examples": examples_unseen,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--zip-path", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("rna_family_audit_output"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.zip_path:
        archive_path = args.zip_path
    else:
        archive_path = Path(tempfile.gettempdir()) / "seq_cls_data.zip"
        download(args.url, archive_path)

    archive_sha256 = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    groups: dict[str, dict[str, list[Record]]] = collections.defaultdict(lambda: {"train": [], "test": []})
    members: list[dict[str, object]] = []

    with zipfile.ZipFile(archive_path) as archive:
        for info in archive.infolist():
            members.append({"name": info.filename, "uncompressed_bytes": info.file_size})
            if info.is_dir() or not info.filename.lower().endswith(FASTA_SUFFIXES):
                continue
            role = role_for(info.filename)
            if role is None:
                continue
            with archive.open(info) as handle:
                records = list(read_fasta(handle))
            groups[group_key(info.filename)][role].extend(records)

    pair_summaries = []
    for key, split in sorted(groups.items()):
        if split["train"] and split["test"]:
            pair_summaries.append(summarize_pair(split["train"], split["test"], key))

    summary = {
        "audit_version": VERSION,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_url": args.url,
        "archive_sha256": archive_sha256,
        "archive_bytes": archive_path.stat().st_size,
        "fasta_pairs_found": len(pair_summaries),
        "pairs": pair_summaries,
        "preregistered_gate": {
            "leakage": "test_seen_family_fraction >= 0.90 and family_lookup_oracle_accuracy >= 0.95",
            "feasibility": "at least 10 classes feasible for family holdout",
        },
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    (args.out_dir / "archive_members.json").write_text(json.dumps(members, indent=2), encoding="utf-8")
    (args.out_dir / "README.md").write_text(
        "# ncRNA family-overlap audit\n\n"
        "Deterministic audit of train/test Rfam-family overlap. No model is trained and no novelty is claimed.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
