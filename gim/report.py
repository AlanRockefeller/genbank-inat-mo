"""Join GenBank references with the observation lookups and write the report:
links.csv (every reference), problems.csv, missed_by_search.csv (when a FASTA
cross-check was loaded), summary.txt and index.html."""

from __future__ import annotations

import csv
import html
import json
import logging
import re
import sqlite3
from collections import Counter
from pathlib import Path

log = logging.getLogger(__name__)

LOW_RANKS = {
    "genus", "subgenus", "section", "subsection", "series", "complex", "species",
    "subspecies", "variety", "form", "hybrid", "infrahybrid", "stirps", "group",
}
_SKIP_WORDS = {"uncultured", "cf.", "aff.", "nr.", "unidentified", "fungal", "fungi"}

# Status codes (CSV `status` column) and the plain-English label shown in the
# HTML report, summary.txt and the CSV `status_text` column. Every code says
# what it is about: the accession, or the observation itself.
STATUS_TEXT = {
    "accession_listed": "OK: observation lists this accession",
    "accession_in_text": "OK: accession is in the observation's notes or comments",
    "accession_missing_has_sequence": "Observation exists and has a sequence, but not this GenBank accession",
    "accession_missing_no_sequence": "Observation exists, but has no sequence or GenBank accession",
    "observation_not_found": "Observation does not exist (wrong number in GenBank, or deleted)",
    "observation_not_returned": "MO knows the number but returned no observation (merged into another?)",
    "lookup_error": "Lookup failed; run verify again",
    "unchecked": "Not looked up yet; run verify",
}
PROBLEM = {
    "accession_missing_has_sequence",
    "accession_missing_no_sequence",
    "observation_not_found",
    "observation_not_returned",
    # Not checked yet, or the lookup failed: listed so a clean-looking report
    # never hides observations nobody has actually checked.
    "lookup_error",
    "unchecked",
}
# obs.status (as stored by verify) -> report status.
_OBS_STATUS = {"not_found": "observation_not_found", "not_returned": "observation_not_returned", "error": "lookup_error"}


def what_to_do(status: str, accession: str, site: str, taxon_check: str, doubtful: bool) -> str:
    """One-line suggested action for a report row."""
    acc = accession.split(".")[0]
    field = "a 'Genbank Accession Number' observation field" if site == "inat" else "the observation's sequence record"
    if doubtful:
        return "Probably not an MO number (MO is also a herbarium code); check the GenBank record"
    if status == "accession_missing_has_sequence":
        return f"Add {acc} to {field}"
    if status == "accession_missing_no_sequence":
        return f"Add the sequence and {acc} to {field}"
    if status == "observation_not_found":
        return "Check the observation number in the GenBank record (typo, or the observation was deleted)"
    if status == "observation_not_returned":
        return "Find the observation it was merged into and add the accession there"
    if taxon_check == "different_genus":
        return "Accession is linked; check the identification or the observation number"
    if status in ("lookup_error", "unchecked"):
        return "Run gim verify"
    return ""


def obs_url(site: str, oid: int) -> str:
    return (
        f"https://www.inaturalist.org/observations/{oid}"
        if site == "inat"
        else f"https://mushroomobserver.org/obs/{oid}"
    )


def genus(name: str) -> str:
    for w in re.split(r"\s+", (name or "").strip()):
        if not w or w.lower() in _SKIP_WORDS:
            continue
        w = w.strip("[]'\"")
        return w if w[:1].isupper() else ""
    return ""


def taxon_check(gb_organism: str, obs_taxon: str, rank: str) -> str:
    g1, g2 = genus(gb_organism), genus(obs_taxon)
    if not g1 or not g2:
        return "not_comparable"
    if rank and rank.lower() not in LOW_RANKS:
        return "obs_above_genus"
    return "same_genus" if g1 == g2 else "different_genus"


def build_rows(con: sqlite3.Connection) -> list[dict]:
    rows = []
    q = """
        SELECT r.accession, rec.version, rec.organism, rec.definition, r.site, r.obs_id,
               group_concat(DISTINCT r.kind) AS kinds,
               group_concat(r.field, ' | ') AS fields,
               min(r.snippet) AS snippet,
               o.status, o.taxon, o.rank, o.seq_accessions, o.text_accessions, o.has_sequence, o.detail
        FROM refs r
        JOIN records rec ON rec.accession = r.accession
        LEFT JOIN obs o ON o.site = r.site AND o.obs_id = r.obs_id
        GROUP BY r.accession, r.site, r.obs_id
        ORDER BY r.site, r.obs_id, r.accession
    """
    for x in con.execute(q):
        acc = x["accession"].upper()
        seq = set((x["seq_accessions"] or "").split())
        txt = set((x["text_accessions"] or "").split())
        st = x["status"]
        if st is None:
            status = "unchecked"
        elif st != "found":
            status = _OBS_STATUS.get(st, "lookup_error")
        elif acc in seq:
            status = "accession_listed"
        elif acc in txt:
            status = "accession_in_text"
        elif x["has_sequence"]:
            status = "accession_missing_has_sequence"
        else:
            status = "accession_missing_no_sequence"
        kinds = set((x["kinds"] or "").split(","))
        kind = "mo_explicit" if "mo_explicit" in kinds else ("inat" if "inat" in kinds else "mo_prefix")
        tcheck = taxon_check(x["organism"], x["taxon"] or "", x["rank"] or "") if st == "found" else ""
        # A bare "MO 123456" on a record whose MO observation is missing, or is
        # a different genus without the accession, is probably not an MO number
        # (e.g. a Missouri Botanical Garden herbarium number).
        doubtful = kind == "mo_prefix" and (
            status == "observation_not_found"
            or (status.startswith("accession_missing") and tcheck == "different_genus")
        )
        detail = json.loads(x["detail"] or "{}")
        accession = x["version"] or x["accession"]
        rows.append(
            {
                "accession": accession,
                "organism": x["organism"],
                "site": x["site"],
                "obs_id": x["obs_id"],
                "obs_url": obs_url(x["site"], x["obs_id"]),
                "status": status,
                "status_text": STATUS_TEXT[status],
                "what_to_do": what_to_do(status, accession, x["site"], tcheck, doubtful),
                "obs_taxon": x["taxon"] or "",
                "taxon_check": tcheck,
                "match_kind": kind,
                "doubtful": "yes" if doubtful else "",
                "obs_accessions": " ".join(sorted(seq)),
                "obs_text_accessions": " ".join(sorted(txt)),
                "obs_user": detail.get("user") or "",
                "genbank_fields": x["fields"],
                "genbank_text": x["snippet"],
                "definition": x["definition"],
            }
        )
    return rows


def write(con: sqlite3.Connection, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    rows = build_rows(con)
    fields = list(rows[0].keys()) if rows else ["accession"]

    def dump(name: str, rs: list[dict]) -> None:
        with open(out / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rs)

    problems = [r for r in rows if r["status"] in PROBLEM or r["taxon_check"] == "different_genus"]
    dump("links.csv", rows)
    dump("problems.csv", problems)

    # Recall cross-check against FASTA titles, if loaded: a title ref is
    # missed when no search reached its record (fetch then queued it as
    # 'fasta'), or the record was fetched but the full text didn't yield it.
    missed = [
        dict(r)
        for r in con.execute(
            """SELECT t.accession, t.site, t.obs_id, t.kind,
                      CASE WHEN EXISTS (SELECT 1 FROM candidates c
                                        WHERE c.source = 'fasta' AND c.id = t.accession)
                           THEN 'not found by any search'
                           ELSE 'record fetched, number not extracted' END AS reason
               FROM title_refs t
               WHERE EXISTS (SELECT 1 FROM candidates c WHERE c.source = 'fasta' AND c.id = t.accession)
                  OR NOT EXISTS (SELECT 1 FROM refs r WHERE r.accession = t.accession
                                 AND r.site = t.site AND r.obs_id = t.obs_id)"""
        )
    ]
    n_title = con.execute("SELECT count(*) FROM title_refs").fetchone()[0]
    if n_title:
        with open(out / "missed_by_search.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["accession", "site", "obs_id", "kind", "reason"])
            w.writeheader()
            w.writerows(missed)

    stats = {
        "records_scanned": con.execute("SELECT count(*) FROM records").fetchone()[0],
        "records_with_refs": con.execute("SELECT count(*) FROM records WHERE n_refs > 0").fetchone()[0],
        "links": len(rows),
        "inat_obs": len({r["obs_id"] for r in rows if r["site"] == "inat"}),
        "mo_obs": len({r["obs_id"] for r in rows if r["site"] == "mo"}),
        "by_status": {
            site: dict(Counter(r["status"] for r in rows if r["site"] == site)) for site in ("inat", "mo")
        },
        "different_genus": sum(r["taxon_check"] == "different_genus" for r in rows),
        "doubtful_mo_prefix": sum(bool(r["doubtful"]) for r in rows),
        "fasta_title_refs": n_title,
        "fasta_title_refs_missed": len(missed),
    }
    lines = [
        f"GenBank records scanned:      {stats['records_scanned']}",
        f"Records citing iNat/MO:       {stats['records_with_refs']}",
        f"Accession-observation links:  {stats['links']}",
        f"Distinct iNat observations:   {stats['inat_obs']}",
        f"Distinct MO observations:     {stats['mo_obs']}",
        "",
    ]
    for site, label in (("inat", "iNaturalist"), ("mo", "Mushroom Observer")):
        lines.append(f"{label}:")
        for k in STATUS_TEXT:
            n = stats["by_status"][site].get(k, 0)
            if n:
                lines.append(f"  {n:7d}  {STATUS_TEXT[k]}")
        lines.append("")
    lines.append(f"Taxon differs at genus level:     {stats['different_genus']}")
    lines.append(f"Doubtful bare-'MO' matches:       {stats['doubtful_mo_prefix']}")
    if n_title:
        lines.append(
            f"FASTA-title refs not found by the search: {len(missed)} of {n_title}"
            " (see missed_by_search.csv)"
        )
    (out / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out / "index.html").write_text(_html(stats, lines, problems), encoding="utf-8")
    return stats


SITE_TEXT = {"inat": "iNaturalist", "mo": "Mushroom Observer"}
TAXON_TEXT = {
    "same_genus": "Same genus",
    "different_genus": "Different genus",
    "obs_above_genus": "Observation only identified above genus",
    "not_comparable": "Can't compare",
    "": "",
}
KIND_TEXT = {
    "inat": "iNaturalist number",
    "mo_explicit": "Named as Mushroom Observer / MUOB",
    "mo_prefix": "Bare “MO” prefix",
}

# (header, function of row -> cell HTML)
def _link(url: str, text) -> str:
    return f'<a href="{html.escape(url)}">{html.escape(str(text))}</a>'


_COLUMNS = [
    ("GenBank record", lambda r: _link(f"https://www.ncbi.nlm.nih.gov/nuccore/{r['accession']}", r["accession"])),
    ("GenBank organism", lambda r: html.escape(r["organism"] or "")),
    ("Site", lambda r: SITE_TEXT[r["site"]]),
    ("Observation", lambda r: _link(r["obs_url"], r["obs_id"])),
    ("Problem", lambda r: html.escape(r["status_text"])),
    ("What to do", lambda r: html.escape(r["what_to_do"])),
    ("Observation taxon", lambda r: html.escape(r["obs_taxon"])),
    ("Genus check", lambda r: TAXON_TEXT.get(r["taxon_check"], r["taxon_check"])),
    ("Accessions on observation", lambda r: html.escape(r["obs_accessions"] or "none")),
    ("Observer", lambda r: html.escape(r["obs_user"])),
    ("How GenBank cites it", lambda r: html.escape(f"{KIND_TEXT.get(r['match_kind'], r['match_kind'])}: {r['genbank_text']}")),
]


def _html(stats: dict, summary_lines: list[str], problems: list[dict]) -> str:
    def row(r: dict) -> str:
        attrs = (
            f' data-site="{r["site"]}" data-status="{r["status"]}"'
            f' data-taxon="{r["taxon_check"]}" data-doubtful="{"1" if r["doubtful"] else ""}"'
        )
        return f"<tr{attrs}>" + "".join(f"<td>{fn(r)}</td>" for _, fn in _COLUMNS) + "</tr>"

    body = "\n".join(row(r) for r in problems)
    head = "".join(f"<th>{h}</th>" for h, _ in _COLUMNS)
    status_opts = "".join(
        f'<option value="{k}">{html.escape(STATUS_TEXT[k])}</option>' for k in STATUS_TEXT if k in PROBLEM
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GenBank iNat/MO check</title>
<style>
:root {{ --bg:#fff; --fg:#1d2125; --muted:#5f6b76; --line:#d9dee3; --head:#f3f5f7; --link:#1a5fb4; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#16191c; --fg:#e4e7ea; --muted:#9aa5b0; --line:#30363c; --head:#1f2428; --link:#78aeed; }} }}
body {{ background:var(--bg); color:var(--fg); font:14px/1.45 system-ui,sans-serif; margin:0 16px 40px; }}
h1 {{ font-size:20px; margin:20px 0 8px; }}
pre {{ color:var(--muted); font-size:13px; }}
a {{ color:var(--link); }}
input, select {{ font:inherit; padding:4px 6px; margin:0 8px 8px 0; background:var(--bg); color:var(--fg); border:1px solid var(--line); border-radius:4px; }}
.wrap {{ overflow-x:auto; }}
table {{ border-collapse:collapse; font-size:13px; }}
th, td {{ border-bottom:1px solid var(--line); padding:4px 8px; text-align:left; vertical-align:top; }}
th {{ background:var(--head); position:sticky; top:0; cursor:pointer; }}
td:last-child {{ color:var(--muted); max-width:420px; }}
td:nth-child(5), td:nth-child(6) {{ min-width:220px; }}
.note {{ color:var(--muted); max-width:760px; }}
</style></head><body>
<h1>GenBank records citing iNaturalist / Mushroom Observer</h1>
<pre>{html.escape(chr(10).join(summary_lines))}</pre>
<h1>Problems ({len(problems)})</h1>
<p class="note">Each row is one GenBank record that cites an observation. “Problem” says what is wrong
with the observation's link back to that record; a different genus is listed even when the
accession is linked. Click a column header to sort.</p>
<input id="q" placeholder="Filter text…" size="30">
<select id="site"><option value="">All sites</option><option value="inat">iNaturalist</option><option value="mo">Mushroom Observer</option></select>
<select id="st"><option value="">All problems</option>{status_opts}<option value="different_genus">Different genus</option></select>
<label><input type="checkbox" id="hide"> hide doubtful bare-“MO” matches</label>
<div class="wrap"><table id="t"><thead><tr>{head}</tr></thead><tbody>
{body}
</tbody></table></div>
<script>
const rows=[...document.querySelectorAll('#t tbody tr')];
function f(){{const q=document.getElementById('q').value.toLowerCase(),s=document.getElementById('site').value,
st=document.getElementById('st').value,h=document.getElementById('hide').checked;
for(const r of rows){{const d=r.dataset;let ok=(!s||d.site===s)&&(!h||!d.doubtful);
if(st)ok=ok&&(d.status===st||d.taxon===st);
if(q)ok=ok&&r.textContent.toLowerCase().includes(q);r.style.display=ok?'':'none';}}}}
for(const id of ['q','site','st','hide'])document.getElementById(id).addEventListener('input',f);
document.querySelectorAll('#t th').forEach((th,i)=>th.addEventListener('click',()=>{{
const tb=document.querySelector('#t tbody'),asc=th.dataset.asc!=='1';th.dataset.asc=asc?'1':'0';
rows.sort((a,b)=>{{const x=a.cells[i].textContent,y=b.cells[i].textContent,nx=+x,ny=+y;
const d=(!isNaN(nx)&&!isNaN(ny)&&x&&y)?nx-ny:x.localeCompare(y);return asc?d:-d;}});rows.forEach(r=>tb.appendChild(r));}}));
</script></body></html>
"""
