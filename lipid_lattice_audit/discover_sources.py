#!/usr/bin/env python3
"""Discover public native coarse-resolution lipid MS2 sources."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import pathlib
import re
import urllib.request
import zipfile
import io


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "lipid-source-discovery/0.1"})
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read()


def lower_field(row: dict, names):
    values = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        if name.lower() in values:
            return values[name.lower()]
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/lipid_source_discovery")
    args = parser.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    library_url = "https://external.gnps2.org/gnpslibrary.json"
    raw = fetch(library_url)
    (outdir / "gnpslibrary.json").write_bytes(raw)
    payload = json.loads(raw)
    rows = payload if isinstance(payload, list) else payload.get("libraries", payload.get("data", []))
    if not isinstance(rows, list):
        raise TypeError(f"unexpected GNPS library list payload: {type(rows)!r}")
    lipid_rows = []
    for row in rows:
        if isinstance(row, str):
            name = row
            normalized = {"libraryname": row}
        elif isinstance(row, dict):
            name = str(lower_field(row, ["libraryname", "library_name", "name"]) or "")
            normalized = row
        else:
            continue
        if "LIPID" in name.upper():
            lipid_rows.append(normalized)
    (outdir / "lipid_libraries.json").write_text(
        json.dumps(lipid_rows, indent=2), encoding="utf-8"
    )

    repo_tree_url = "https://api.github.com/repos/LipiTUM/lipidetective/git/trees/main?recursive=1"
    tree_raw = fetch(repo_tree_url)
    (outdir / "lipidetective_tree.json").write_bytes(tree_raw)
    tree = json.loads(tree_raw).get("tree", [])
    data_like = [
        {"path": item.get("path"), "size": item.get("size"), "type": item.get("type")}
        for item in tree
        if re.search(r"(?:\.hdf5$|\.h5$|\.mzML$|\.mgf$|\.msp$|\.json$)", item.get("path", ""), re.I)
    ]
    substantive_data = [
        item for item in data_like
        if not item["path"].startswith("tests/") and not item["path"].startswith("models/")
    ]

    article_url = "https://academic.oup.com/bib/article/27/4/bbag378/8742281"
    supp_manifest = []
    supp_text_files = []
    zip_match = None
    article_error = None
    try:
        article_raw = fetch(article_url)
        article_text = article_raw.decode("utf-8", errors="replace")
        (outdir / "lipidetective_article.html").write_bytes(article_raw)
        zip_match = re.search(
            r'https?://[^"\\']+supplementary-material_bbag378\\.zip[^"\\']*',
            html.unescape(article_text),
        )
    except Exception as exc:
        article_error = f"{type(exc).__name__}: {exc}"
        (outdir / "lipidetective_article_error.txt").write_text(
            article_error + "\n", encoding="utf-8"
        )

    if zip_match:
        supp_url = zip_match.group(0).replace("&amp;", "&")
        try:
            supp_raw = fetch(supp_url)
            (outdir / "lipidetective_supplementary.zip").write_bytes(supp_raw)
            supp_manifest.append({
                "url": supp_url,
                "bytes": len(supp_raw),
                "sha256": hashlib.sha256(supp_raw).hexdigest(),
            })
            with zipfile.ZipFile(io.BytesIO(supp_raw)) as archive:
                for info in archive.infolist():
                    supp_manifest.append({
                        "path": info.filename,
                        "bytes": info.file_size,
                        "compressed_bytes": info.compress_size,
                    })
                    if info.filename.lower().endswith((".txt", ".csv", ".tsv", ".md")):
                        data = archive.read(info)
                        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", info.filename)
                        target = outdir / "supplementary_text" / safe
                        target.parent.mkdir(parents=True, exist_ok=True)
                        target.write_bytes(data)
                        supp_text_files.append(str(target.relative_to(outdir)))
        except Exception as exc:
            supp_manifest.append({
                "url": supp_url,
                "error": f"{type(exc).__name__}: {exc}",
            })
    else:
        supp_manifest.append({
            "error": "supplementary zip link not found or article fetch blocked",
            "article_error": article_error,
        })

    result = {
        "version": "0.1.0",
        "gnps_library_list": {
            "url": library_url,
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "total_rows": len(rows),
            "lipid_rows": len(lipid_rows),
            "lipid_library_names": [
                str(lower_field(row, ["libraryname", "library_name", "name"]) or row)
                if isinstance(row, dict) else str(row)
                for row in lipid_rows
            ],
        },
        "lipidetective_repository": {
            "tree_url": repo_tree_url,
            "data_like_files": data_like,
            "substantive_public_data_files": substantive_data,
        },
        "lipidetective_supplementary": supp_manifest,
        "extracted_plaintext_files": supp_text_files,
    }
    (outdir / "discovery.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    report = [
        "# Native coarse lipid-source discovery",
        "",
        f"- GNPS libraries total: {len(rows):,}",
        f"- GNPS library names containing LIPID: {len(lipid_rows):,}",
        f"- LipiDetective repository data-like files: {len(data_like):,}",
        f"- Substantive public training-data files in repository: {len(substantive_data):,}",
        f"- Supplementary archive found: {bool(zip_match)}",
        "",
        "## GNPS lipid libraries",
        "",
    ]
    report.extend(
        f"- {name}" for name in result["gnps_library_list"]["lipid_library_names"]
    )
    report.extend([
        "",
        "This discovery step does not yet count coarse spectra. It enumerates candidate sources for the next audit.",
    ])
    (outdir / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
