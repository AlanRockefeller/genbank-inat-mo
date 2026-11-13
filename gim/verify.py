"""Look up the cited observations on iNaturalist and Mushroom Observer.

Read-only: only public GET endpoints are used, nothing is ever written to
either site. Results are cached in the `obs` table, so a rerun only asks for
observations it has not seen (or all of them with --refresh / --max-age)."""

from __future__ import annotations

import html
import json
import logging
import re
import sqlite3
import time

from .extract import accessions_in
from .net import HttpError, Throttle, request

log = logging.getLogger(__name__)

INAT_API = "https://api.inaturalist.org/v1/observations"
MO_API = "https://mushroomobserver.org/api2/observations"

_TAGS = re.compile(r"<[^>]+>")
_DNA = re.compile(r"^[ACGTURYKMSWBDHVN\s-]{80,}$", re.I)


def _plain(s) -> str:
    return html.unescape(_TAGS.sub(" ", str(s or ""))).strip()


def _todo(con: sqlite3.Connection, site: str, refresh: bool, max_age_days: float | None) -> list[int]:
    ids = [r[0] for r in con.execute("SELECT DISTINCT obs_id FROM refs WHERE site = ? ORDER BY obs_id", (site,))]
    if refresh:
        return ids
    if max_age_days is None:
        known = {r[0] for r in con.execute("SELECT obs_id FROM obs WHERE site = ? AND status != 'error'", (site,))}
    else:
        known = {
            r[0]
            for r in con.execute(
                "SELECT obs_id FROM obs WHERE site = ? AND status != 'error'"
                " AND fetched_at >= datetime('now', ?)",
                (site, f"-{max_age_days} days"),
            )
        }
    return [i for i in ids if i not in known]


def _save(con, site, obs_id, status, taxon="", rank="", seq_acc=(), text_acc=(), has_seq=0, detail=None):
    con.execute(
        "INSERT OR REPLACE INTO obs (site, obs_id, status, taxon, rank, seq_accessions, text_accessions,"
        " has_sequence, detail, fetched_at) VALUES (?,?,?,?,?,?,?,?,?, datetime('now'))",
        (
            site,
            obs_id,
            status,
            taxon,
            rank,
            " ".join(sorted(seq_acc)),
            " ".join(sorted(text_acc)),
            int(has_seq),
            json.dumps(detail or {}, ensure_ascii=False),
        ),
    )


# ---- iNaturalist ------------------------------------------------------------
def _inat_summary(r: dict) -> dict:
    taxon = r.get("taxon") or {}
    seq_acc: set[str] = set()
    text_acc: set[str] = set()
    has_seq = False
    ofvs = []
    for o in r.get("ofvs") or []:
        name, value = str(o.get("name") or ""), str(o.get("value") or "")
        if _DNA.match(value) or "dna barcode" in name.lower():
            has_seq = has_seq or bool(value.strip())
            ofvs.append((name, f"<{len(value)} bp>"))
            continue
        accs = accessions_in(value)
        # Any observation field holding an accession counts as structured:
        # "Genbank Accession Number", "GenBank number (URL)", "Accession", ...
        if accs and re.search(r"genbank|accession|ncbi|sequence", name, re.I):
            seq_acc |= accs
        else:
            text_acc |= accs
        ofvs.append((name, value[:200]))
    text_acc |= accessions_in(r.get("description") or "")
    for c in r.get("comments") or []:
        text_acc |= accessions_in(c.get("body") or "")
    return {
        "taxon": taxon.get("name") or "",
        "rank": taxon.get("rank") or "",
        "seq_acc": seq_acc,
        "text_acc": text_acc - seq_acc,
        "has_seq": has_seq or bool(seq_acc),
        "detail": {
            "user": (r.get("user") or {}).get("login"),
            "quality": r.get("quality_grade"),
            "observed_on": r.get("observed_on"),
            "place": r.get("place_guess"),
            "ofvs": ofvs,
        },
    }


def verify_inat(con: sqlite3.Connection, refresh=False, max_age_days=None, batch=200) -> None:
    todo = _todo(con, "inat", refresh, max_age_days)
    log.info("iNaturalist: %d observations to look up", len(todo))
    throttle = Throttle(1.1)  # iNat asks for about 1 request/second
    for i in range(0, len(todo), batch):
        ids = todo[i : i + batch]
        _, body = request(
            INAT_API,
            params={"id": ",".join(map(str, ids)), "per_page": str(len(ids))},
            throttle=throttle,
            timeout=120,
        )
        data = json.loads(body)
        found = {}
        for r in data.get("results") or []:
            found[int(r["id"])] = r
        for oid in ids:
            r = found.get(oid)
            if r is None:
                _save(con, "inat", oid, "not_found")
            else:
                s = _inat_summary(r)
                _save(con, "inat", oid, "found", s["taxon"], s["rank"], s["seq_acc"], s["text_acc"], s["has_seq"], s["detail"])
        con.commit()
        log.info("iNaturalist: %d/%d looked up", min(i + batch, len(todo)), len(todo))


# ---- Mushroom Observer -------------------------------------------------------
_MO_MISSING = re.compile(r"Observation #(\d+) does not exist")


def _mo_summary(r: dict) -> dict:
    cons = r.get("consensus") or {}
    if not isinstance(cons, dict):
        cons = {"name": str(cons)}
    seq_acc: set[str] = set()
    has_seq = False
    seqs = []
    for s in r.get("sequences") or []:
        has_seq = True
        acc = str(s.get("accession") or "")
        seq_acc |= accessions_in(acc.upper())
        seqs.append({k: s.get(k) for k in ("locus", "archive", "accession")})
    text = [_plain(r.get("notes"))]
    nf = r.get("notes_fields")
    if isinstance(nf, dict):
        text += [_plain(v) for v in nf.values()]
    elif isinstance(nf, list):
        text += [_plain(x) for x in nf]
    for s in r.get("sequences") or []:
        text.append(_plain(s.get("notes")))
    for c in r.get("comments") or []:
        text += [_plain(c.get("summary")), _plain(c.get("content"))]
    for x in r.get("external_links") or []:
        text.append(json.dumps(x))
    text_acc = set()
    for t in text:
        text_acc |= accessions_in(t)
    return {
        "taxon": cons.get("name") or "",
        "rank": str(cons.get("rank") or "").lower(),
        "seq_acc": seq_acc,
        "text_acc": text_acc - seq_acc,
        "has_seq": has_seq,
        "detail": {
            "user": (r.get("owner") or {}).get("login_name") if isinstance(r.get("owner"), dict) else r.get("owner"),
            "date": r.get("date"),
            "location": (r.get("location") or {}).get("name") if isinstance(r.get("location"), dict) else r.get("location"),
            "sequences": seqs,
            "specimen_available": r.get("specimen_available"),
        },
    }


def verify_mo(con: sqlite3.Connection, refresh=False, max_age_days=None, batch=100) -> None:
    todo = _todo(con, "mo", refresh, max_age_days)
    log.info("Mushroom Observer: %d observations to look up (~%.0f min at MO's 5 s/request)",
             len(todo), len(todo) / batch * 5 / 60)
    throttle = Throttle(5.0)  # MO allows 20 anonymous requests a minute
    i = 0
    while i < len(todo):
        ids = todo[i : i + batch]
        while ids:
            status, body = request(
                MO_API,
                params={"id": ",".join(map(str, ids)), "detail": "high", "format": "json"},
                throttle=throttle,
                timeout=180,
                ok_statuses=(400, 404, 422, 500),
            )
            try:
                data = json.loads(body)
            except ValueError:
                raise HttpError(status, body, MO_API) from None
            run_time = float(data.get("run_time") or 0)
            if run_time > 5:
                throttle.push_back(run_time)  # MO: wait at least the last run_time
            missing = set()
            other_errors = []
            for e in data.get("errors") or []:
                m = _MO_MISSING.search(str(e.get("details") or ""))
                if m:
                    missing.add(int(m.group(1)))
                else:
                    other_errors.append(e.get("details") or e.get("code"))
            if missing:
                for oid in missing:
                    _save(con, "mo", oid, "not_found")
                con.commit()
                ids = [x for x in ids if x not in missing]
                continue  # the whole request failed; ask again without them
            if other_errors:
                log.error("MO error for ids %s..: %s", ids[:3], other_errors[:2])
                for oid in ids:
                    _save(con, "mo", oid, "error", detail={"error": str(other_errors[:2])})
                con.commit()
                break
            results = {int(r["id"]): r for r in data.get("results") or []}
            pages = int(data.get("number_of_pages") or 1)
            for page in range(2, pages + 1):
                _, body = request(
                    MO_API,
                    params={"id": ",".join(map(str, ids)), "detail": "high", "format": "json", "page": str(page)},
                    throttle=throttle,
                    timeout=180,
                )
                for r in json.loads(body).get("results") or []:
                    results[int(r["id"])] = r
            for oid in ids:
                r = results.get(oid)
                if r is None:
                    # Exists but not returned: MO lists only the primary
                    # observation of a merged occurrence.
                    _save(con, "mo", oid, "not_returned")
                else:
                    s = _mo_summary(r)
                    _save(con, "mo", oid, "found", s["taxon"], s["rank"], s["seq_acc"], s["text_acc"], s["has_seq"], s["detail"])
            con.commit()
            break
        i += batch
        log.info("Mushroom Observer: %d/%d looked up", min(i, len(todo)), len(todo))
