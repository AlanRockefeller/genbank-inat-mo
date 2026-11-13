"""HTTP with retries and per-host pacing (standard library only)."""

from __future__ import annotations

import logging
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

log = logging.getLogger(__name__)

USER_AGENT = "genbank-sequence-miner/0.1 (GenBank voucher cross-check; read-only; +https://github.com/AlanRockefeller/genbank-sequence-miner)"
RETRY_STATUSES = {408, 429, 500, 502, 503, 504}


class Throttle:
    """Keeps at least `interval` seconds between request starts, across threads."""

    def __init__(self, interval: float):
        self.interval = interval
        self._next = 0.0
        self._lock = threading.Lock()

    def wait(self, interval: float | None = None) -> None:
        with self._lock:
            now = time.monotonic()
            start = max(now, self._next)
            self._next = start + (self.interval if interval is None else interval)
        if start > now:
            time.sleep(start - now)

    def push_back(self, seconds: float) -> None:
        """Delay every later request, e.g. after a 429 or a slow server reply."""
        with self._lock:
            self._next = max(self._next, time.monotonic() + seconds)


class HttpError(Exception):
    def __init__(self, status: int, body: bytes, url: str):
        super().__init__(f"HTTP {status} from {url.split('?')[0]}: {body[:300]!r}")
        self.status = status
        self.body = body


def request(
    url: str,
    *,
    params: dict | None = None,
    data: dict | None = None,
    throttle: Throttle | None = None,
    retries: int = 6,
    timeout: float = 300,
    ok_statuses: tuple[int, ...] = (),
) -> tuple[int, bytes]:
    """GET (or POST when `data` is given). Returns (status, body).

    Retries network errors and 408/429/5xx with backoff. Any other non-2xx
    status raises HttpError unless listed in `ok_statuses`.
    """
    if params:
        url = url + ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    delay = 5.0
    for attempt in range(retries + 1):
        if throttle:
            throttle.wait()
        req = urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as e:
            payload = e.read()
            if e.code in ok_statuses:
                return e.code, payload
            if e.code not in RETRY_STATUSES or attempt == retries:
                raise HttpError(e.code, payload, url) from None
            wait = delay
            retry_after = e.headers.get("Retry-After")
            if retry_after and retry_after.isdigit():
                wait = max(wait, float(retry_after))
            log.warning("HTTP %s from %s; retrying in %.0fs", e.code, url.split("?")[0], wait)
            if throttle:
                throttle.push_back(wait)
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            if attempt == retries:
                raise
            wait = delay
            log.warning("%s from %s; retrying in %.0fs", e, url.split("?")[0], wait)
        time.sleep(wait)
        delay = min(delay * 2, 300)
    raise AssertionError("unreachable")
