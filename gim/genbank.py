"""Minimal streaming parser for GenBank flat files (the text efetch returns
for rettype=gb). Keeps only what the voucher search needs: identifiers,
organism, definition, comment, reference titles and feature qualifiers.
The sequence is skipped."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Iterator

_WS = re.compile(r"\s+")


@dataclass
class Record:
    accession: str = ""
    version: str = ""
    length: int = 0
    definition: str = ""
    organism: str = ""
    # (label, text) pairs to scan, e.g. ("DEFINITION", ...), ("source/specimen_voucher", ...)
    fields: list[tuple[str, str]] = field(default_factory=list)

    def qualifier(self, name: str, feature: str = "source") -> str:
        label = f"{feature}/{name}"
        for k, v in self.fields:
            if k == label:
                return v
        return ""


def _join(lines: list[str]) -> str:
    return _WS.sub(" ", " ".join(s.strip() for s in lines)).strip()


def parse(lines: Iterable[str]) -> Iterator[Record]:
    """Parse GenBank flat-file text, given as an iterable of lines."""
    rec: Record | None = None
    keyword = ""  # current top-level keyword (DEFINITION, COMMENT, ...)
    sub = ""  # current sub-keyword under REFERENCE / SOURCE
    buf: list[str] = []
    in_features = False
    in_origin = False
    feature = ""
    qual_name = ""
    qual_buf: list[str] = []
    in_quote = False

    def flush_block() -> None:
        nonlocal buf
        if rec is None or not buf:
            buf = []
            return
        text = _join(buf)
        if keyword == "DEFINITION":
            rec.definition = text.rstrip(".")
            rec.fields.append(("DEFINITION", rec.definition))
        elif keyword == "COMMENT":
            rec.fields.append(("COMMENT", text))
        elif keyword == "KEYWORDS" and text not in (".", ""):
            rec.fields.append(("KEYWORDS", text))
        elif keyword == "REFERENCE" and sub == "TITLE":
            rec.fields.append(("REFERENCE/TITLE", text))
        elif keyword == "SOURCE" and sub == "ORGANISM":
            rec.organism = buf[0].strip()
        elif keyword == "ACCESSION":
            rec.accession = text.split()[0] if text else ""
        elif keyword == "VERSION":
            rec.version = text.split()[0] if text else ""
        buf = []

    def flush_qual() -> None:
        nonlocal qual_name, qual_buf
        if rec is not None and qual_name and qual_name != "translation":
            value = _join(qual_buf)
            if value.startswith('"') and value.endswith('"') and len(value) >= 2:
                value = value[1:-1]
            rec.fields.append((f"{feature}/{qual_name}", value.replace('""', '"')))
        qual_name, qual_buf = "", []

    for raw in lines:
        line = raw.rstrip("\r\n")
        if line.startswith("//"):
            if rec is not None:
                flush_qual()
                flush_block()
                yield rec
            rec, keyword, sub, buf = None, "", "", []
            in_features = in_origin = in_quote = False
            continue
        if line.startswith("LOCUS"):
            rec = Record()
            m = re.search(r"\s(\d+) (?:bp|aa)", line)
            rec.length = int(m.group(1)) if m else 0
            keyword = "LOCUS"
            continue
        if rec is None or in_origin:
            continue
        if line.startswith("ORIGIN") or line.startswith("CONTIG"):
            flush_qual()
            flush_block()
            in_origin = True
            continue
        if in_features:
            if line and not line[0].isspace():
                flush_qual()
                in_features = False  # fall through to keyword handling
            else:
                key = line[5:21].strip() if len(line) > 5 else ""
                body = line[21:]
                if key and not in_quote:
                    flush_qual()
                    feature = key
                    continue
                if body.startswith("/") and not in_quote:
                    flush_qual()
                    name, _, value = body[1:].partition("=")
                    qual_name, qual_buf = name, [value]
                    in_quote = value.count('"') % 2 == 1
                else:
                    qual_buf.append(body)
                    if body.count('"') % 2 == 1:
                        in_quote = not in_quote
                continue
        if line.startswith("FEATURES"):
            flush_block()
            keyword, in_features, feature = "FEATURES", True, ""
            continue
        if line and not line[0].isspace():
            flush_block()
            keyword = line[:12].strip()
            sub = ""
            buf = [line[12:]]
        elif line[:12].strip():  # sub-keyword, e.g. "  ORGANISM", "  TITLE"
            flush_block()
            sub = line[:12].strip()
            buf = [line[12:]]
        else:
            buf.append(line[12:])
    # A record without its closing // was cut off mid-response: dropped, so a
    # partial record never replaces a stored one. Callers see a short batch
    # and retry.


def parse_text(text: str) -> list[Record]:
    return list(parse(text.splitlines()))
