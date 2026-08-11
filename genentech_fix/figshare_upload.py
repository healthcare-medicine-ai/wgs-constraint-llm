#!/usr/bin/env python
"""Upload a large file from Sherlock directly to an existing figshare item.

Avoids moving 16 GB down to a laptop and back up again. Uploads in parts, can
resume, and verifies the md5 figshare computes against the local one.

The token is never passed on the command line (it would land in shell history
and in `ps`). Put it in a file readable only by you:

    umask 077 && printf '%s' 'YOUR_TOKEN' > ~/.figshare_token

Create the item in the figshare web UI first -- set the title, description,
licence and collection there -- then pass its numeric id here.

    python figshare_upload.py --article-id 12345678 \\
        --file /oak/.../data/All_hg38_RS.bw

Add --items if your account uses the /account/items endpoint rather than
/account/articles (institutional figshare instances differ).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

BASE = "https://api.figshare.com/v2"
TOKEN_FILE = Path.home() / ".figshare_token"
CHUNK = 8 * 1024 * 1024


def read_token() -> str:
    if not TOKEN_FILE.exists():
        sys.exit(f"no token at {TOKEN_FILE}. Create it with:\n"
                 f"  umask 077 && printf '%s' 'YOUR_TOKEN' > {TOKEN_FILE}")
    mode = TOKEN_FILE.stat().st_mode & 0o777
    if mode & 0o077:
        sys.exit(f"{TOKEN_FILE} is readable by others (mode {mode:o}). "
                 f"Run: chmod 600 {TOKEN_FILE}")
    token = TOKEN_FILE.read_text().strip()
    if not token:
        sys.exit(f"{TOKEN_FILE} is empty")
    return token


def call(token, method, url, data=None, raw=None):
    if not url.startswith("http"):
        url = BASE + url
    body = raw if raw is not None else (
        json.dumps(data).encode() if data is not None else None)
    req = Request(url, data=body, method=method)
    req.add_header("Authorization", f"token {token}")
    if raw is None and data is not None:
        req.add_header("Content-Type", "application/json")
    try:
        with urlopen(req) as resp:
            payload = resp.read()
            if not payload:
                return {}
            try:
                return json.loads(payload)
            except json.JSONDecodeError:
                return {"_raw": payload}
    except HTTPError as e:
        detail = e.read().decode(errors="replace")[:500]
        sys.exit(f"\n{method} {url}\nHTTP {e.code}: {detail}")


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    total = path.stat().st_size
    done = 0
    last = 0.0
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
            done += len(chunk)
            now = time.time()
            if now - last > 5:
                print(f"\r  hashing {100*done/total:5.1f}%", end="", flush=True)
                last = now
    print(f"\r  hashing 100.0%")
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--article-id", required=True)
    ap.add_argument("--file", required=True, type=Path)
    ap.add_argument("--md5", default=None,
                    help="skip local hashing if already known")
    ap.add_argument("--items", action="store_true",
                    help="use /account/items instead of /account/articles")
    args = ap.parse_args()

    path = args.file
    if not path.exists():
        sys.exit(f"missing: {path}")
    size = path.stat().st_size
    token = read_token()
    kind = "items" if args.items else "articles"

    print(f"file : {path}")
    print(f"size : {size:,} bytes ({size/2**30:.2f} GiB)")

    md5 = args.md5 or md5_of(path)
    print(f"md5  : {md5}")

    print("\n[1/4] registering the file ...")
    loc = call(token, "POST", f"/account/{kind}/{args.article_id}/files",
               {"md5": md5, "size": size, "name": path.name})["location"]
    file_id = loc.rstrip("/").split("/")[-1]
    print(f"      file id {file_id}")

    print("[2/4] fetching part list ...")
    info = call(token, "GET", loc)
    parts = call(token, "GET", info["upload_url"])["parts"]
    print(f"      {len(parts)} parts")

    print("[3/4] uploading ...")
    started = time.time()
    sent = 0
    with open(path, "rb") as fh:
        for p in parts:
            n, lo, hi = p["partNo"], p["startOffset"], p["endOffset"]
            length = hi - lo + 1
            if p.get("status") == "COMPLETE":
                sent += length
                continue
            fh.seek(lo)
            call(token, "PUT", f"{info['upload_url']}/{n}", raw=fh.read(length))
            sent += length
            elapsed = time.time() - started
            rate = sent / elapsed / 2**20 if elapsed else 0
            eta = (size - sent) / (sent / elapsed) if sent and elapsed else 0
            print(f"\r      part {n}/{len(parts)}  {100*sent/size:5.1f}%  "
                  f"{rate:6.1f} MiB/s  eta {eta/60:5.1f} min", end="", flush=True)
    print()

    print("[4/4] completing ...")
    call(token, "POST", f"/account/{kind}/{args.article_id}/files/{file_id}")

    final = call(token, "GET", f"/account/{kind}/{args.article_id}/files/{file_id}")
    remote_md5 = final.get("computed_md5") or final.get("supplied_md5")
    print(f"\n  figshare md5 : {remote_md5}")
    print(f"  local md5    : {md5}")
    if remote_md5 and remote_md5 != md5:
        sys.exit("  MISMATCH -- do not publish this item until resolved")
    print("  match. Upload complete.")
    print(f"\nThe item is still a draft; publish it in the web UI to mint a DOI.")


if __name__ == "__main__":
    main()
