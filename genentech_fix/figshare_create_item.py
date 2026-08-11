#!/usr/bin/env python
"""Create a DRAFT figshare item and print its id.

Deliberately does not publish. A draft is private and deletable; publishing
mints a DOI and is irreversible, so that stays a human action in the web UI.

Reads the token from ~/.figshare_token and never prints it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

BASE = "https://api.figshare.com/v2"
TOKEN_FILE = Path.home() / ".figshare_token"

TITLE = "GERP RS conservation scores for GRCh38 (lifted from hg19)"

DESCRIPTION = """\
<p>GERP RS (Genomic Evolutionary Rate Profiling, rejected substitutions)
conservation scores for GRCh38, lifted over from the hg19 track
(<code>All_hg19_RS.bw</code>) using UCSC's <code>hg19ToHg38.over.chain.gz</code>
with default settings.</p>

<p>Used as a moderator in the unified meta-regression model of
AJHG-D-24-00613, "A unified meta-regression model identifies genes associated
with epilepsy".</p>

<p><b>Verification.</b> Positional correctness was checked against the hg19
source through the chain: of 1,500 sampled positions where GERP varies sharply
between neighbouring bases, 1,498 carry the source value at the chain-mapped
coordinate at zero offset, and none at an offset of &plusmn;1 or &plusmn;2.</p>

<p><b>File properties.</b> bigWig, 16,358,613,092 bytes.
Coverage 2,864,325,273 bases. Values range &minus;12.3 to 6.17.
UCSC-style chromosome naming, including hg38 alt contigs.
md5 <code>10a03ee7969d3000ffd2e6e6f84f2453</code>.</p>

<p>Deposited so that reproductions of the published analysis use the identical
input rather than an independently derived liftover, which can differ in
coverage depending on the tool and chain used.</p>
"""

KEYWORDS = ["GERP", "conservation", "GRCh38", "hg38", "liftOver",
            "genomic constraint", "epilepsy", "bigWig"]

REFERENCES = ["https://doi.org/10.6084/m9.figshare.27184245"]


def read_token() -> str:
    if not TOKEN_FILE.exists():
        sys.exit(f"no token at {TOKEN_FILE}")
    if TOKEN_FILE.stat().st_mode & 0o077:
        sys.exit(f"{TOKEN_FILE} is readable by others; run chmod 600")
    return TOKEN_FILE.read_text().strip()


def call(token, method, url, data=None):
    if not url.startswith("http"):
        url = BASE + url
    body = json.dumps(data).encode() if data is not None else None
    req = Request(url, data=body, method=method)
    req.add_header("Authorization", f"token {token}")
    if body:
        req.add_header("Content-Type", "application/json")
    try:
        with urlopen(req) as r:
            raw = r.read()
            return json.loads(raw) if raw else {}
    except HTTPError as e:
        sys.exit(f"\n{method} {url}\nHTTP {e.code}: "
                 f"{e.read().decode(errors='replace')[:600]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--items", action="store_true",
                    help="use /account/items rather than /account/articles")
    ap.add_argument("--update-id", default=None,
                    help="populate an existing draft instead of creating a new "
                         "one, so a draft started in the UI is not duplicated")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    kind = "items" if args.items else "articles"
    token = read_token()

    payload = {
        "title": TITLE,
        "description": DESCRIPTION,
        "keywords": KEYWORDS,
        "references": REFERENCES,
        "defined_type": "dataset",
    }

    verb = "Updating" if args.update_id else "Creating"
    print(f"Would {verb.lower()}:" if args.dry_run else f"{verb} draft item:")
    print(f"  endpoint : /account/{kind}"
          f"{'/' + args.update_id if args.update_id else ''}")
    print(f"  title    : {TITLE}")
    print(f"  keywords : {', '.join(KEYWORDS)}")
    if args.dry_run:
        return

    if args.update_id:
        existing = call(token, "GET", f"/account/{kind}/{args.update_id}")
        if existing.get("published_date"):
            sys.exit(f"  {args.update_id} is already published; refusing to "
                     f"overwrite its metadata")
        call(token, "PUT", f"/account/{kind}/{args.update_id}", payload)
        article_id = args.update_id
        print(f"\n  updated  : {article_id}")
    else:
        created = call(token, "POST", f"/account/{kind}", payload)
        loc = created.get("location", "")
        article_id = loc.rstrip("/").split("/")[-1]
        print(f"\n  created  : {loc}")
    print(f"  ITEM ID  : {article_id}")

    detail = call(token, "GET", f"/account/{kind}/{article_id}")
    print(f"  status   : {detail.get('status')}   "
          f"(published: {bool(detail.get('published_date'))})")
    print(f"\n  Edit at: https://figshare.com/account/articles/{article_id}")
    print("  Not published. Set categories and licence in the UI, then publish "
          "there when ready.")


if __name__ == "__main__":
    main()
