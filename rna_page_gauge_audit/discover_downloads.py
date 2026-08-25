#!/usr/bin/env python3
"""Discover the official bpRNA bulk-download URLs without trusting mirrors."""

from __future__ import annotations

import hashlib
import html
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

PAGE = "https://bprna.cgrb.oregonstate.edu/download.php"
OUT = Path("rna_page_gauge_audit/discovery")
OUT.mkdir(parents=True, exist_ok=True)

req = urllib.request.Request(PAGE, headers={"User-Agent": "Mozilla/5.0 research-audit/0.1"})
with urllib.request.urlopen(req, timeout=60) as response:
    raw = response.read()
    final_url = response.geturl()
page_text = raw.decode("utf-8", errors="replace")
(OUT / "official_download_page.html").write_text(page_text, encoding="utf-8")

anchors = []
for match in re.finditer(r"<a\b([^>]*)>(.*?)</a>", page_text, flags=re.I | re.S):
    attrs, body = match.groups()
    hm = re.search(r"""href\s*=\s*(['"])(.*?)\1""", attrs, flags=re.I | re.S)
    if not hm:
        continue
    href = html.unescape(hm.group(2).strip())
    label = re.sub(r"<[^>]+>", " ", body)
    label = " ".join(html.unescape(label).split())
    anchors.append({
        "label": label,
        "href": href,
        "absolute_url": urllib.parse.urljoin(final_url, href),
    })

payload = {
    "audit_version": "0.1.1",
    "page_requested": PAGE,
    "page_final_url": final_url,
    "page_bytes": len(raw),
    "page_sha256": hashlib.sha256(raw).hexdigest(),
    "anchors": anchors,
}
(OUT / "discovery.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
