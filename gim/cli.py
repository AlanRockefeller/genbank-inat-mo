"""gim: find iNaturalist / Mushroom Observer numbers in GenBank and check them."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from . import db, harvest, report, verify
from .ncbi import Ncbi

log = logging.getLogger("gim")


def _ncbi(args) -> Ncbi:
    key = os.environ.get("NCBI_API_KEY")
    if not key and args.api_key_file:
        key = Path(args.api_key_file).expanduser().read_text().strip()
    if not key:
        log.info("No NCBI API key (NCBI_API_KEY): limited to 3 requests/second")
    return Ncbi(api_key=key or None, email=os.environ.get("NCBI_EMAIL") or None)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="gim", description=__doc__)
    p.add_argument("--db", default="data/gim.sqlite", help="SQLite state file (default: %(default)s)")
    p.add_argument("--api-key-file", help="file holding the NCBI API key (or set NCBI_API_KEY)")
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search", help="run the Entrez keyword and wildcard searches")
    s.add_argument("--refresh", action="store_true", help="rerun searches already done")
    s.add_argument("--only", choices=["keywords", "wildcards"], help="run one kind only")

    f = sub.add_parser("fetch", help="download candidate records as GenBank text and extract numbers")
    f.add_argument("--batch", type=int, default=200)
    f.add_argument("--workers", type=int, default=3)

    sc = sub.add_parser("scan", help="fetch EVERY record matching an Entrez query (deep recall check)")
    sc.add_argument("query")
    sc.add_argument("--batch", type=int, default=500)
    sc.add_argument("--workers", type=int, default=3)
    sc.add_argument("--refresh", action="store_true", help="start over instead of resuming")

    fa = sub.add_parser("fasta", help="cross-check: extract numbers from a FASTA file's titles and queue those records")
    fa.add_argument("path")

    v = sub.add_parser("verify", help="look up cited observations on iNaturalist and MO (read-only)")
    v.add_argument("--site", choices=["inat", "mo"], help="one site only")
    v.add_argument("--refresh", action="store_true", help="look up every observation again")
    v.add_argument("--max-age", type=float, metavar="DAYS", help="re-check lookups older than DAYS")

    r = sub.add_parser("report", help="write CSV/HTML report")
    r.add_argument("--out", default="report")

    a = sub.add_parser("run", help="search + fetch + verify + report")
    a.add_argument("--out", default="report")
    a.add_argument("--fasta", help="also cross-check against this FASTA file's titles")

    args = p.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname).1s %(message)s",
        datefmt="%H:%M:%S",
    )
    con = db.connect(args.db)

    if args.cmd == "search":
        h = harvest.Harvester(con, _ncbi(args), refresh=args.refresh)
        if args.only != "wildcards":
            for word, flt in harvest.KEYWORDS:
                h.run_keyword(word, flt)
        if args.only != "keywords":
            for root in harvest.ROOTS:
                h.run_root(root)
    elif args.cmd == "fetch":
        harvest.Harvester(con, _ncbi(args)).fetch_candidates(args.batch, args.workers)
    elif args.cmd == "scan":
        harvest.Harvester(con, _ncbi(args), refresh=args.refresh).scan(args.query, args.batch, args.workers)
    elif args.cmd == "fasta":
        n_titles, n_refs = harvest.load_fasta_titles(con, args.path)
        log.info("FASTA: %d titles cite iNat/MO (%d refs); run 'gim fetch' to download them", n_titles, n_refs)
    elif args.cmd == "verify":
        if args.site in (None, "inat"):
            verify.verify_inat(con, args.refresh, args.max_age)
        if args.site in (None, "mo"):
            verify.verify_mo(con, args.refresh, args.max_age)
    elif args.cmd == "report":
        report.write(con, Path(args.out))
        sys.stdout.write((Path(args.out) / "summary.txt").read_text())
    elif args.cmd == "run":
        h = harvest.Harvester(con, _ncbi(args))
        h.search_all()
        if args.fasta:
            harvest.load_fasta_titles(con, args.fasta)
        h.fetch_candidates()
        verify.verify_inat(con)
        verify.verify_mo(con)
        report.write(con, Path(args.out))
        sys.stdout.write((Path(args.out) / "summary.txt").read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
