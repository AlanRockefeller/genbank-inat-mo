"""Find iNaturalist and Mushroom Observer observation numbers in GenBank text.

Each match carries a `kind` saying how sure it is:

  inat         "iNaturalist", "iNat", "naturalista", or an inaturalist.org URL
  mo_explicit  "Mushroom Observer", "MushroomObserver", "MUOB", or a
               mushroomobserver.org URL
  mo_prefix    bare "MO" before the number (MO123456, MO_123456, MO # 123456).
               Usually MO, but MO is also the Index Herbariorum code of the
               Missouri Botanical Garden herbarium, so verification decides.

MO observation numbers on sequenced collections are 6 digits, rarely 5; only
those lengths are accepted. iNat numbers are accepted at 5 to 10 digits
(observations below 10000 date from 2008 and are not plausible vouchers;
the lower bound also keeps years such as "iNat 2021" out).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Words and separators allowed between the site name and the number:
# "iNaturalist # 98565406", "iNat obs. no. 123456", "MUOB:123456",
# "mushroomobserver.org/observations/123456", "Mushroom Observer (MO) 123456".
_GAP_WORDS = (
    r"observations?|obs|show_observation|observer|records?|id|no|nr|num(?:ber)?|"
    r"www|org|https?|mo"
)
_GAP = rf"(?:[\s#:=_./,()\[\]\-]|(?:{_GAP_WORDS})(?![a-z])){{0,10}}"

_INAT = re.compile(
    rf"(?<![a-z0-9])(?:inaturalist|inat|naturalista)(?:\.[a-z]{{2,3}}){{0,3}}{_GAP}([1-9]\d{{4,9}})(?![\d])",
    re.I,
)
_MO_EXPLICIT = re.compile(
    rf"(?<![a-z0-9])(?:mushroom[\s_\-]*observer|muob)(?:\.org)?{_GAP}([1-9]\d{{4,5}})(?![\d])",
    re.I,
)
_MO_PREFIX = re.compile(
    r"(?<![a-z0-9])mo(?:[\s#:_\-.]{0,3}|\s*obs(?:ervation)?\.?\s*(?:#|no\.?)?\s*)([1-9]\d{4,5})(?![\d])",
    re.I,
)

# Nucleotide accessions: 1 letter + 5 digits, 2 letters + 6 or 8 digits,
# RefSeq (NR_123456), and WGS/TSA (4 letters + 8-10 digits, 6 letters + 9-11 digits).
ACCESSION = re.compile(
    r"\b([A-Z]{1,2}\d{5,8}|[A-Z]{2}_\d{6,9}|[A-Z]{4}\d{8,10}|[A-Z]{6}\d{9,11})(?:\.\d+)?\b"
)


@dataclass(frozen=True)
class Ref:
    site: str  # "inat" or "mo"
    obs_id: int
    kind: str  # inat | mo_explicit | mo_prefix
    field: str
    snippet: str


_ZIP_BEFORE = re.compile(r",\s*$")
_ZIP_AFTER = re.compile(r"^(?:-\d{4}\b|,?\s*(?:USA|U\.S\.A\.|United States)\b)")


def _is_zip(text: str, m: re.Match) -> bool:
    """'St. Louis, MO 63108, USA': a Missouri ZIP code in an address (NCBI puts
    submitter addresses in COMMENT and REFERENCE), not an MO observation."""
    if len(m.group(1)) != 5:
        return False
    return bool(_ZIP_BEFORE.search(text[max(0, m.start() - 3) : m.start()]) or _ZIP_AFTER.match(text[m.end() :]))


def _snippet(text: str, m: re.Match, width: int = 40) -> str:
    a = max(0, m.start() - width)
    b = min(len(text), m.end() + width)
    return ("…" if a else "") + text[a:b] + ("…" if b < len(text) else "")


def find_in_text(text: str, field: str = "") -> list[Ref]:
    refs: list[Ref] = []
    explicit_spans: list[tuple[int, int]] = []
    for m in _INAT.finditer(text):
        refs.append(Ref("inat", int(m.group(1)), "inat", field, _snippet(text, m)))
    for m in _MO_EXPLICIT.finditer(text):
        refs.append(Ref("mo", int(m.group(1)), "mo_explicit", field, _snippet(text, m)))
        explicit_spans.append(m.span(1))
    for m in _MO_PREFIX.finditer(text):
        if _is_zip(text, m):
            continue
        if any(a <= m.start(1) < b for a, b in explicit_spans):
            continue  # "Mushroom Observer MO 123456": already counted as explicit
        refs.append(Ref("mo", int(m.group(1)), "mo_prefix", field, _snippet(text, m)))
    return refs


def find_in_fields(fields: list[tuple[str, str]]) -> list[Ref]:
    """All refs in a record's fields, one per (site, obs_id, field); the
    strongest kind wins when a number is matched more than one way."""
    rank = {"mo_explicit": 0, "inat": 0, "mo_prefix": 1}
    best: dict[tuple[str, int, str], Ref] = {}
    for label, text in fields:
        for r in find_in_text(text, label):
            key = (r.site, r.obs_id, r.field)
            if key not in best or rank[r.kind] < rank[best[key].kind]:
                best[key] = r
    return list(best.values())


def accessions_in(text: str) -> set[str]:
    """Version-less, upper-case GenBank-style accessions found in free text."""
    return {m.group(1) for m in ACCESSION.finditer((text or "").upper())}
