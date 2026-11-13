"""NCBI E-utilities client: esearch (with history), UID lists, GenBank text."""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

from .net import Throttle, request

log = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "genbank-inat-mo"


@dataclass
class Search:
    term: str
    count: int
    webenv: str = ""
    query_key: str = ""
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    ids: list[str] | None = None  # set when every UID came back with the search

    @property
    def truncated(self) -> bool:
        return any("used only the first" in w for w in self.warnings)


class EutilsError(Exception):
    pass


class Ncbi:
    def __init__(self, api_key: str | None = None, email: str | None = None):
        self.api_key = api_key
        self.email = email
        # NCBI allows 3 requests/s without a key and 10/s with one.
        # Stay a little under both: NCBI counts per IP and answers bursts with 429.
        self.throttle = Throttle(0.12 if api_key else 0.4)

    def _call(self, endpoint: str, params: dict, timeout: float = 600) -> bytes:
        data = {"tool": TOOL, **params}
        if self.api_key:
            data["api_key"] = self.api_key
        if self.email:
            data["email"] = self.email
        _, body = request(EUTILS + endpoint, data=data, throttle=self.throttle, timeout=timeout)
        return body

    def esearch(self, term: str, history: bool = True, inline: int = 0) -> Search:
        """With `inline`, up to that many UIDs come back in the same call
        (Search.ids), saving a request for small result sets."""
        params = {"db": "nuccore", "term": term, "retmax": str(inline)}
        if history:
            params["usehistory"] = "y"
        for attempt in range(4):
            root = ET.fromstring(self._call("esearch.fcgi", params))
            err = root.findtext("ERROR")
            if err:
                if attempt < 3:
                    log.warning("esearch error %r for %s; retrying", err, term)
                    continue
                raise EutilsError(f"esearch {term!r}: {err}")
            s = Search(
                term=term,
                count=int(root.findtext("Count") or 0),
                webenv=root.findtext("WebEnv") or "",
                query_key=root.findtext("QueryKey") or "",
            )
            s.warnings = [e.text or "" for e in root.iterfind("WarningList/*")]
            s.errors = [f"{e.tag}: {e.text}" for e in root.iterfind("ErrorList/*")]
            ids = [e.text for e in root.iterfind("IdList/Id") if e.text]
            if inline and len(ids) == s.count:
                s.ids = ids
            return s
        raise AssertionError("unreachable")

    def uids(self, s: Search, page: int = 10000):
        """Yield every UID of a history search."""
        if s.ids is not None:
            yield from s.ids
            return
        got = 0
        while got < s.count:
            body = self._call(
                "efetch.fcgi",
                {
                    "db": "nuccore",
                    "WebEnv": s.webenv,
                    "query_key": s.query_key,
                    "rettype": "uilist",
                    "retmode": "text",
                    "retstart": str(got),
                    "retmax": str(page),
                },
            )
            ids = [x for x in body.decode().split() if x.isdigit()]
            if not ids:
                raise EutilsError(f"empty UID page at {got}/{s.count} for {s.term!r}")
            yield from ids
            got += len(ids)

    def genbank_by_ids(self, ids: list[str]) -> str:
        body = self._call(
            "efetch.fcgi",
            {"db": "nuccore", "id": ",".join(ids), "rettype": "gb", "retmode": "text"},
        )
        return self._check_text(body)

    def genbank_by_history(self, s: Search, start: int, size: int) -> str:
        body = self._call(
            "efetch.fcgi",
            {
                "db": "nuccore",
                "WebEnv": s.webenv,
                "query_key": s.query_key,
                "rettype": "gb",
                "retmode": "text",
                "retstart": str(start),
                "retmax": str(size),
            },
        )
        return self._check_text(body)

    @staticmethod
    def _check_text(body: bytes) -> str:
        text = body.decode("utf-8", errors="replace")
        head = text[:2000]
        if "<ERROR>" in head or head.lstrip().startswith("{\"error\""):
            raise EutilsError(head[:300])
        return text
