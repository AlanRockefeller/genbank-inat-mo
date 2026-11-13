"""SQLite store. Every step is resumable: searches, fetched batches and
verified observations are recorded, and a rerun skips what is done."""

from __future__ import annotations

import sqlite3
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS searches (
    term        TEXT PRIMARY KEY,   -- full Entrez query
    purpose     TEXT,               -- keyword | wildcard | wildcard-exact | scan
    count       INTEGER,
    truncated   INTEGER DEFAULT 0,  -- wildcard hit the 600-variation limit
    progress    INTEGER DEFAULT 0,  -- scan: records fetched so far
    done        INTEGER DEFAULT 0,
    updated_at  TEXT DEFAULT (datetime('now'))
);
-- Entrez UIDs (or accession.version strings) still to fetch as GenBank text.
CREATE TABLE IF NOT EXISTS candidates (
    id          TEXT PRIMARY KEY,
    source      TEXT,               -- first search term (or 'fasta') that found it
    fetched     INTEGER DEFAULT 0
);
CREATE TABLE IF NOT EXISTS records (
    accession   TEXT PRIMARY KEY,   -- without version
    version     TEXT,
    organism    TEXT,
    definition  TEXT,
    length      INTEGER,
    n_refs      INTEGER,
    scanned_at  TEXT DEFAULT (datetime('now'))
);
CREATE TABLE IF NOT EXISTS refs (
    accession   TEXT,
    site        TEXT,               -- inat | mo
    obs_id      INTEGER,
    kind        TEXT,               -- inat | mo_explicit | mo_prefix
    field       TEXT,
    snippet     TEXT,
    PRIMARY KEY (accession, site, obs_id, field)
);
CREATE INDEX IF NOT EXISTS refs_obs ON refs (site, obs_id);
-- Refs found in the titles of a local FASTA file (recall cross-check).
CREATE TABLE IF NOT EXISTS title_refs (
    accession   TEXT,
    site        TEXT,
    obs_id      INTEGER,
    kind        TEXT,
    PRIMARY KEY (accession, site, obs_id)
);
CREATE TABLE IF NOT EXISTS obs (
    site        TEXT,
    obs_id      INTEGER,
    status      TEXT,               -- found | not_found | not_returned | error
    taxon       TEXT,
    rank        TEXT,
    seq_accessions TEXT,            -- space-separated accessions in structured fields
    text_accessions TEXT,           -- accessions mentioned in notes/description/comments
    has_sequence INTEGER,
    detail      TEXT,               -- compact JSON for the report
    fetched_at  TEXT DEFAULT (datetime('now')),
    PRIMARY KEY (site, obs_id)
);
"""


def connect(path: str | Path) -> sqlite3.Connection:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Long busy timeout: a second gim command (e.g. `fasta` during `search`) waits
    # for the other's transaction instead of failing.
    con = sqlite3.connect(path, timeout=300)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.executescript(SCHEMA)
    return con


def store_records(con: sqlite3.Connection, records, find_refs) -> tuple[int, int]:
    """Store parsed records and their refs. Returns (records, refs)."""
    n_refs = 0
    for rec in records:
        acc = (rec.accession or rec.version).split(".")[0]
        if not acc:
            continue
        refs = find_refs(rec.fields)
        con.execute(
            "INSERT OR REPLACE INTO records (accession, version, organism, definition, length, n_refs)"
            " VALUES (?,?,?,?,?,?)",
            (acc, rec.version, rec.organism, rec.definition[:500], rec.length, len(refs)),
        )
        con.execute("DELETE FROM refs WHERE accession = ?", (acc,))
        con.executemany(
            "INSERT OR REPLACE INTO refs VALUES (?,?,?,?,?,?)",
            [(acc, r.site, r.obs_id, r.kind, r.field, r.snippet) for r in refs],
        )
        n_refs += len(refs)
    return len(records), n_refs
