"""Collect GenBank records that may cite an iNaturalist or Mushroom Observer
observation, fetch them as full GenBank text, and extract the references.

Why several kinds of search: Entrez indexes "MO # 346490", "MO_559147" and
"MUOB:123456" as the word "mo"/"muob" (plus the phrase "mo 346490"), so a
keyword search finds them. A number written straight after the letters
("MO411618", "iNat69337073") is a single word of its own, which only a
wildcard finds, and Entrez expands a wildcard to at most 600 words. So each
wildcard root is split into longer prefixes (mo4*, mo41*, mo411*, ...) until
Entrez stops reporting truncation."""

from __future__ import annotations

import gzip
import logging
import re
import sqlite3
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

from . import extract, genbank
from .db import store_records
from .ncbi import Ncbi, Search

log = logging.getLogger(__name__)

# Bare "MO" is ambiguous outside fungi, so MO-only searches are limited to
# fungi and slime moulds (Eumycetozoa), which is everything MO records.
MO_FILTER = "(fungi[organism] OR mycetozoa[organism])"
# WGS contigs repeat the voucher on every contig: keep the WGS master record only.
EXCLUDE = '(wgs[filter] NOT "wgs master"[prop])'

KEYWORDS: list[tuple[str, str | None]] = [
    ("inaturalist", None),
    ("inat", None),
    ("naturalista", None),
    ("mushroomobserver", None),
    ("mushroom AND observer", None),
    ("muob", None),
    ("mo", MO_FILTER),
]

DIGITS = "0123456789"
ANY = DIGITS + "abcdefghijklmnopqrstuvwxyz./:_-"


@dataclass(frozen=True)
class Root:
    word: str  # e.g. "mo": enumerates mo1*, mo2*, ... (or word* when first is "")
    filter: str | None
    lengths: range  # digit counts worth an exact-term query when a prefix is split
    charset: str = DIGITS
    first: str = "123456789"  # first characters after the word


ROOTS = [
    Root("mo", MO_FILTER, range(5, 7)),
    Root("muob", None, range(5, 7)),
    Root("inat", None, range(5, 11)),
    Root("inaturalist", None, range(5, 11)),
    Root("naturalista", None, range(5, 11)),
    Root("mushroomobserver", None, range(0), charset=ANY, first=""),
]


def full_term(term: str, organism_filter: str | None) -> str:
    t = f"({term})"
    if organism_filter:
        t += f" AND {organism_filter}"
    return f"{t} NOT {EXCLUDE}"


class Harvester:
    def __init__(self, con: sqlite3.Connection, ncbi: Ncbi, refresh: bool = False):
        self.con = con
        self.ncbi = ncbi
        self.refresh = refresh

    # ---- searches -------------------------------------------------------
    def _done(self, term: str) -> sqlite3.Row | None:
        if self.refresh:
            return None
        return self.con.execute(
            "SELECT * FROM searches WHERE term = ? AND done = 1", (term,)
        ).fetchone()

    def _collect(self, s: Search, purpose: str) -> int:
        new = 0
        if s.count:
            for chunk in _chunks(self.ncbi.uids(s), 5000):
                cur = self.con.executemany(
                    "INSERT OR IGNORE INTO candidates (id, source) VALUES (?, ?)",
                    [(u, s.term) for u in chunk],
                )
                new += cur.rowcount
        self.con.execute(
            "INSERT OR REPLACE INTO searches (term, purpose, count, truncated, done, updated_at)"
            " VALUES (?, ?, ?, ?, 1, datetime('now'))",
            (s.term, purpose, s.count, int(s.truncated)),
        )
        self.con.commit()
        return new

    def run_keyword(self, word: str, organism_filter: str | None) -> None:
        term = full_term(word, organism_filter)
        if self._done(term):
            return
        s = self.ncbi.esearch(term)
        new = self._collect(s, "keyword")
        log.info("keyword %-24s %8d records, %7d new", word, s.count, new)

    def _probe(self, term: str) -> tuple[Search, list[str]]:
        """Worker thread: one wildcard search plus its UIDs (unless truncated)."""
        s = self.ncbi.esearch(term, inline=5000)
        if s.truncated or not s.count:
            return s, []
        return s, list(self.ncbi.uids(s))

    def run_root(self, root: Root, workers: int = 4) -> None:
        """Enumerate root+prefix* until no wildcard is truncated, one prefix
        length at a time, with the searches of each level run in parallel
        (the shared throttle keeps the total within NCBI's rate)."""
        level = [root.word + c + "*" for c in root.first] if root.first else [root.word + "*"]
        n_queries = n_records = 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            while level:
                todo, nxt = [], []
                for word in level:
                    row = self._done(full_term(word, root.filter))
                    if row is None:
                        todo.append(word)
                    elif row["truncated"]:
                        nxt += self._children(root, word)
                if todo:
                    log.info("root %s: %d searches at prefix length %d",
                             root.word, len(todo), len(todo[0].rstrip("*")) - len(root.word))
                for word, (s, ids) in zip(todo, pool.map(lambda w: self._probe(full_term(w, root.filter)), todo)):
                    n_queries += 1
                    self.con.executemany(
                        "INSERT OR IGNORE INTO candidates (id, source) VALUES (?, ?)",
                        [(u, s.term) for u in ids],
                    )
                    # A truncated search's UIDs are a subset of its children's:
                    # record the split and collect nothing from it.
                    self.con.execute(
                        "INSERT OR REPLACE INTO searches (term, purpose, count, truncated, done, updated_at)"
                        " VALUES (?, 'wildcard', ?, ?, 1, datetime('now'))",
                        (s.term, s.count, int(s.truncated)),
                    )
                    if s.truncated:
                        nxt += self._children(root, word)
                    else:
                        n_records += s.count
                    if n_queries % 100 == 0:
                        self.con.commit()
                self.con.commit()
                level = nxt
        log.info("root %-18s %5d new queries, %7d records matched", root.word, n_queries, n_records)

    @staticmethod
    def _children(root: Root, word: str) -> list[str]:
        stem = word[:-1] if word.endswith("*") else word
        if not word.endswith("*"):
            return []
        ndigits = len(stem) - len(root.word)
        out = [stem] if ndigits in root.lengths else []  # the exact word, e.g. "mo41161"
        if ndigits >= max(root.lengths, default=12):
            log.warning("%s still truncated at full length; not splitting further", word)
            return out
        return out + [stem + c + "*" for c in root.charset]

    def search_all(self) -> None:
        for word, flt in KEYWORDS:
            self.run_keyword(word, flt)
        for root in ROOTS:
            self.run_root(root)
        n = self.con.execute("SELECT count(*) FROM candidates").fetchone()[0]
        log.info("candidates: %d records", n)

    # ---- fetching --------------------------------------------------------
    def fetch_candidates(self, batch: int = 200, workers: int = 3) -> None:
        self._fetch_pending(batch, workers)
        if self.con.execute("SELECT 1 FROM title_refs LIMIT 1").fetchone():
            n = queue_title_misses(self.con)
            if n:
                log.info("fetch: %d FASTA-title records were not found by the searches; fetching them", n)
                self._fetch_pending(batch, workers)

    def _fetch_pending(self, batch: int, workers: int) -> None:
        todo = [r[0] for r in self.con.execute("SELECT id FROM candidates WHERE fetched = 0 ORDER BY id")]
        if not todo:
            log.info("fetch: nothing to do")
            return
        log.info("fetch: %d records in batches of %d", len(todo), batch)
        batches = [todo[i : i + batch] for i in range(0, len(todo), batch)]
        t0 = time.monotonic()
        done_recs = done_refs = 0

        with ThreadPoolExecutor(max_workers=workers) as pool:
            for i, (ok, missing, recs) in enumerate(pool.map(self._fetch_batch, batches), 1):
                n_rec, n_ref = store_records(self.con, recs, extract.find_in_fields)
                self.con.executemany("UPDATE candidates SET fetched = 1 WHERE id = ?", [(x,) for x in ok])
                # NCBI returned nothing for these, even one at a time
                # (suppressed or withdrawn): recorded as unavailable, not fetched.
                self.con.executemany("UPDATE candidates SET fetched = 2 WHERE id = ?", [(x,) for x in missing])
                if missing:
                    log.warning("fetch: no record returned for %s", ", ".join(missing[:10]))
                self.con.commit()
                done_recs += n_rec
                done_refs += n_ref
                if i % 25 == 0 or i == len(batches):
                    rate = done_recs / max(1e-6, time.monotonic() - t0)
                    left = (len(todo) - i * batch) / max(rate, 1e-6)
                    log.info(
                        "fetch: %d/%d batches, %d records, %d refs (%.0f rec/s, ~%.0f min left)",
                        i, len(batches), done_recs, done_refs, rate, max(0, left) / 60,
                    )

    def _fetch_batch(self, ids: list[str]) -> tuple[list[str], list[str], list[genbank.Record]]:
        """Worker thread: fetch `ids`, returning (fetched ids, ids NCBI has no
        record for, records). Returned records carry accessions, not the
        requested ids, so a short batch is retried once and then split in half
        until the ids that return nothing are isolated."""
        recs = genbank.parse_text(self.ncbi.genbank_by_ids(ids))
        if len(recs) >= len(ids):
            return ids, [], recs
        recs = genbank.parse_text(self.ncbi.genbank_by_ids(ids))
        if len(recs) >= len(ids):
            return ids, [], recs
        if len(ids) == 1:
            return [], ids, []
        log.debug("batch of %d returned %d records; splitting", len(ids), len(recs))
        mid = len(ids) // 2
        a_ok, a_miss, a_recs = self._fetch_batch(ids[:mid])
        b_ok, b_miss, b_recs = self._fetch_batch(ids[mid:])
        return a_ok + b_ok, a_miss + b_miss, a_recs + b_recs

    # ---- deep scan of an arbitrary query ----------------------------------
    def scan(self, query: str, batch: int = 200, workers: int = 3) -> None:
        """Fetch every record matching `query` and extract refs.

        The query's UIDs are saved as candidates first, then fetched like any
        other candidate. Resuming therefore works from that saved list, not
        from an offset into a new search whose results may have changed."""
        row = self.con.execute("SELECT * FROM searches WHERE term = ?", (query,)).fetchone()
        if row is None or not row["done"] or self.refresh:
            s = self.ncbi.esearch(query)
            if s.truncated:
                raise SystemExit(
                    f"scan: Entrez truncated a wildcard in {query!r} ({'; '.join(s.warnings)}). "
                    "Lengthen the wildcard root, or split the query into several scans."
                )
            log.info("scan: %d records match; saving their UIDs", s.count)
            new = self._collect(s, "scan")
            log.info("scan: %d UIDs not seen before", new)
        self.fetch_candidates(batch, workers)


# ---- FASTA titles (cross-check) -------------------------------------------
_FASTA_ACC = re.compile(r"^>(\S+)")
# NCBI FASTA id tags whose next field is an accession(.version).
_SEQID_TAGS = {"gb", "emb", "dbj", "ref", "tpg", "tpe", "tpd", "gpp", "pdb", "lcl"}


def fasta_accession(seqid: str) -> str:
    """'OL342397.1' -> 'OL342397.1'; 'gb|OL342397.1|LOCUS' and
    'gi|123|gb|OL342397.1|' -> 'OL342397.1'. Empty if none is found."""
    if "|" not in seqid:
        return seqid
    parts = seqid.split("|")
    for i, p in enumerate(parts[:-1]):
        if p.lower() in _SEQID_TAGS and parts[i + 1]:
            return parts[i + 1]
    return ""


def load_fasta_titles(con: sqlite3.Connection, path: str) -> tuple[int, int]:
    """Extract refs from FASTA definition lines into title_refs. Nothing is
    queued here: `fetch` queues only the accessions the searches missed
    (queue_title_misses), which is what makes this a recall check.
    Returns (titles with refs, refs)."""
    opener = gzip.open if path.endswith(".gz") else open
    n_titles = n_refs = 0
    rows = []
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = _FASTA_ACC.match(line)
            if not m:
                continue
            refs = extract.find_in_text(line[m.end() :].strip(), "FASTA title")
            if not refs:
                continue
            acc = fasta_accession(m.group(1)).split(".")[0]
            if not acc:
                continue
            n_titles += 1
            for r in refs:
                rows.append((acc, r.site, r.obs_id, r.kind))
            n_refs += len(refs)
    con.executemany("INSERT OR IGNORE INTO title_refs VALUES (?,?,?,?)", rows)
    con.commit()
    return n_titles, n_refs


def queue_title_misses(con: sqlite3.Connection) -> int:
    """After the search candidates are fetched: queue FASTA-title accessions
    that no search reached (source 'fasta'). Returns how many."""
    cur = con.execute(
        """INSERT OR IGNORE INTO candidates (id, source)
           SELECT DISTINCT t.accession, 'fasta' FROM title_refs t
           WHERE NOT EXISTS (SELECT 1 FROM records r WHERE r.accession = t.accession)"""
    )
    con.commit()
    return cur.rowcount


def _chunks(it, n):
    buf = []
    for x in it:
        buf.append(x)
        if len(buf) == n:
            yield buf
            buf = []
    if buf:
        yield buf
