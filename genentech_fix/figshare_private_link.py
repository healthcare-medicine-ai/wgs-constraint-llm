#!/usr/bin/env python
"""Create, list or revoke a figshare private link for an UNPUBLISHED item.

A private link gives a named recipient access without publishing. No DOI is
minted, the item is not indexed, and the link can be revoked at any time. That
is the right mechanism for data whose redistribution terms are unsettled:
publishing mints a permanent DOI and cannot be walked back.

Reads the token from ~/.figshare_token and never prints it.

  python figshare_private_link.py --article-id 123 --create --expires 2026-11-30
  python figshare_private_link.py --article-id 123 --list
  python figshare_private_link.py --article-id 123 --revoke LINK_ID
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


def read_token() -> str:
    if not TOKEN_FILE.exists():
        sys.exit("no token at " + str(TOKEN_FILE))
    if TOKEN_FILE.stat().st_mode & 0o077:
        sys.exit(str(TOKEN_FILE) + " is readable by others; run chmod 600")
    return TOKEN_FILE.read_text().strip()


def call(token, method, url, data=None):
    body = json.dumps(data).encode() if data is not None else None
    req = Request(BASE + url, data=body, method=method)
    req.add_header("Authorization", "token " + token)
    req.add_header("Content-Type", "application/json")
    try:
        with urlopen(req) as response:
            raw = response.read()
    except HTTPError as exc:
        detail = exc.read().decode()[:400]
        sys.exit("HTTP " + str(exc.code) + " on " + method + " " + url + ": " + detail)
    return json.loads(raw) if raw else {}


def as_url(value):
    if not value:
        return ""
    if value.startswith("http"):
        return value
    return "https://figshare.com/s/" + value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--article-id", required=True)
    ap.add_argument("--create", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--revoke", default=None, metavar="LINK_ID")
    ap.add_argument("--expires", default=None,
                    help="YYYY-MM-DD; the link stops working after this date")
    args = ap.parse_args()

    token = read_token()
    path = "/account/articles/" + str(args.article_id) + "/private_links"

    detail = call(token, "GET", "/account/articles/" + str(args.article_id))
    published = bool(detail.get("published_date"))
    print("  item " + str(args.article_id) + ": " + str(detail.get("title", ""))[:60])
    print("  status: " + str(detail.get("status")) + "   published: " + str(published))
    if published:
        print("  NOTE: this item is already public; a private link adds nothing.")

    if args.create:
        payload = {"expires_date": args.expires} if args.expires else {}
        created = call(token, "POST", path, payload)
        # html_location is the shareable URL; "location" is the API endpoint
        # for the link object and is useless to a recipient.
        location = as_url(created.get("html_location") or created.get("location"))
        print()
        print("  PRIVATE LINK CREATED")
        print("  " + location)
        if args.expires:
            print("  expires: " + args.expires)

    if args.revoke:
        call(token, "DELETE", path + "/" + str(args.revoke))
        print()
        print("  revoked link " + str(args.revoke))

    if args.list or args.create or args.revoke:
        links = call(token, "GET", path)
        print()
        print("  " + str(len(links)) + " private link(s) now:")
        for link in links:
            print("    id=" + str(link.get("id"))
                  + "  expires=" + str(link.get("expires_date")))
            print("      " + as_url(link.get("html_location")))


if __name__ == "__main__":
    main()
