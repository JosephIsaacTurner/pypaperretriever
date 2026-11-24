"""HTTP client with robots.txt, Crawl-Delay, backoff, 429 retry, and 403 suppression.

This module provides an :class:`HttpClient` that can optionally fetch and
respect robots.txt directives, including Crawl-Delay, for every domain it
contacts.  It can also automatically retry on HTTP 429 (Too Many Requests),
503, and 504 responses with exponential backoff, and suppress domains that
return HTTP 403 (Forbidden).  When all features are disabled the client is
a thin pass-through to :func:`requests.get`.
"""

from __future__ import annotations

import random
import time
from datetime import datetime
from email.utils import parsedate_to_datetime
from typing import Any, Dict, Optional
from urllib.parse import urlparse
from urllib.robotparser import RobotFileParser

import requests


class RobotsTxtCache:
    """Lazily fetch and cache :class:`RobotFileParser` instances per domain.

    Args:
        user_agent: User-Agent string used for ``can_fetch`` / ``crawl_delay``
            look-ups in the parsed robots.txt.
    """

    def __init__(self, user_agent: str) -> None:
        self.user_agent = user_agent
        self._cache: Dict[str, Optional[RobotFileParser]] = {}

    @staticmethod
    def _domain_key(url: str) -> str:
        """Return ``scheme://netloc`` for *url*."""
        parsed = urlparse(url)
        return f"{parsed.scheme}://{parsed.netloc}"

    def get_parser(self, url: str) -> Optional[RobotFileParser]:
        """Return a :class:`RobotFileParser` for the domain of *url*.

        The parser is fetched and cached on first access.  Returns ``None``
        when the robots.txt cannot be retrieved (treated as fully permissive).
        """
        key = self._domain_key(url)
        if key not in self._cache:
            self._cache[key] = self._fetch_robots_txt(key)
        return self._cache[key]

    @staticmethod
    def _fetch_robots_txt(domain_key: str) -> Optional[RobotFileParser]:
        """Fetch and parse robots.txt for *domain_key*."""
        robots_url = f"{domain_key}/robots.txt"
        rp = RobotFileParser()
        rp.set_url(robots_url)
        try:
            rp.read()
            return rp
        except Exception:
            return None

    def can_fetch(self, url: str) -> bool:
        """Check whether the user-agent is allowed to fetch *url*."""
        parser = self.get_parser(url)
        if parser is None:
            return True
        return parser.can_fetch(self.user_agent, url)

    def crawl_delay(self, url: str) -> Optional[float]:
        """Return the Crawl-Delay for the domain of *url*, or ``None``."""
        parser = self.get_parser(url)
        if parser is None:
            return None
        delay = parser.crawl_delay(self.user_agent)
        if delay is not None:
            return float(delay)
        return None


class HttpClient:
    """HTTP client with robots.txt, backoff, 429 retry, and 403 suppression.

    When *respect_robots_txt* is ``True`` the client:

    1. Checks ``can_fetch()`` before each request and returns ``None`` for
       disallowed URLs.
    2. Applies the Crawl-Delay directive for the target domain if present,
       otherwise falls back to a random delay between *rate_limit_min_s* and
       *rate_limit_max_s*.

    When *delay_min_s* or *delay_max_s* is set the client inserts a random
    per-domain delay between requests, independently of robots.txt.  If a
    Crawl-Delay is also present it takes precedence.

    The client maintains per-domain exponential backoff state.  When a request
    receives HTTP 503, 504, or 429 without delay headers, the backoff delay
    for that domain is increased.  On a successful (2xx) response the backoff
    delay is decreased linearly.  The backoff delay acts as a lower bound on
    the inter-request delay.

    When *honor_retry_after* is ``True`` the client automatically retries
    requests that receive an HTTP 429 (Too Many Requests) response.  The wait
    time is determined from ``Retry-After`` (RFC 7231), ``RateLimit-Reset`` or
    ``X-RateLimit-Reset`` (IETF draft) headers.  If none are present a
    *default_retry_after_s* fallback is used.  Requests that receive HTTP 503
    or 504 are also retried with exponential backoff.

    When *suppress_on_403* is ``True`` (the default) the client records any
    domain that returns HTTP 403 (Forbidden).  Subsequent requests to the same
    domain are short-circuited with a synthetic 403 response, avoiding
    redundant network round-trips during bulk downloads.

    Args:
        user_agent: Default User-Agent string.
        respect_robots_txt: Enable robots.txt checking and Crawl-Delay.
        delay_min_s: Minimum inter-request delay in seconds.  Setting this
            or *delay_max_s* enables per-domain delays.  Defaults to ``None``
            (off).  If only *delay_max_s* is given, this defaults to ``0``.
        delay_max_s: Maximum inter-request delay in seconds.  Defaults to
            ``None`` (off).  If only *delay_min_s* is given, this defaults
            to ``3 * delay_min_s``.
        honor_retry_after: Automatically retry on HTTP 429 responses.
        suppress_on_403: Automatically suppress domains that return 403.
        max_retries: Maximum number of retries per request (applies to 429,
            503, and 504 responses).
        default_retry_after_s: Fallback wait time when no retry header is
            present.
        max_retry_after_s: Cap on the wait time extracted from headers.
        backoff_base_s: Initial backoff delay in seconds.
        backoff_multiplier: Exponential multiplier applied to the backoff
            delay on each consecutive error.
        backoff_max_s: Maximum backoff delay in seconds.
        backoff_decay_s: Amount subtracted from the backoff delay on each
            successful (2xx) response.
        rate_limit_min_s: Minimum fallback delay in seconds (only used as
            a fallback when *respect_robots_txt* is ``True`` and no
            Crawl-Delay is declared).
        rate_limit_max_s: Maximum fallback delay in seconds (same caveat).
        verbosity: Logging verbosity (0=silent, 1=warnings, 2=info, 3=debug).
    """

    def __init__(
        self,
        user_agent: str = "PyPaperRetriever/1.0",
        respect_robots_txt: bool = False,
        delay_min_s: Optional[float] = None,
        delay_max_s: Optional[float] = None,
        honor_retry_after: bool = True,
        suppress_on_403: bool = True,
        max_retries: int = 3,
        default_retry_after_s: float = 60.0,
        max_retry_after_s: float = 300.0,
        backoff_base_s: float = 1.0,
        backoff_multiplier: float = 2.0,
        backoff_max_s: float = 300.0,
        backoff_decay_s: float = 1.0,
        rate_limit_min_s: float = 1.0,
        rate_limit_max_s: float = 3.0,
        verbosity: int = 1,
        polite_mode: bool = True,
    ) -> None:
        self.user_agent = user_agent
        self.respect_robots_txt = respect_robots_txt
        self.honor_retry_after = honor_retry_after
        self.suppress_on_403 = suppress_on_403
        self.max_retries = max_retries
        self.default_retry_after_s = default_retry_after_s
        self.max_retry_after_s = max_retry_after_s
        self.backoff_base_s = backoff_base_s
        self.backoff_multiplier = backoff_multiplier
        self.backoff_max_s = backoff_max_s
        self.backoff_decay_s = backoff_decay_s
        self.rate_limit_min_s = rate_limit_min_s
        self.rate_limit_max_s = rate_limit_max_s
        self.verbosity = verbosity
        self.polite_mode = polite_mode

        # Resolve inter-request delay bounds.
        if delay_min_s is not None or delay_max_s is not None:
            if delay_min_s is None:
                delay_min_s = 0.0
            if delay_max_s is None:
                delay_max_s = 3.0 * delay_min_s
            self.delay_min_s: Optional[float] = delay_min_s
            self.delay_max_s: Optional[float] = delay_max_s
        else:
            self.delay_min_s = None
            self.delay_max_s = None

        self._last_fetch_times: Dict[str, float] = {}
        self._suppressed_domains: Dict[str, str] = {}
        self._backoff_delays: Dict[str, float] = {}
        self._robots_cache: Optional[RobotsTxtCache] = None
        if self.respect_robots_txt:
            self._robots_cache = RobotsTxtCache(user_agent=self.user_agent)

    @staticmethod
    def _domain_key(url: str) -> str:
        """Return ``scheme://netloc`` for per-domain tracking."""
        parsed = urlparse(url)
        return f"{parsed.scheme}://{parsed.netloc}"

    def _check_suppression(self, url: str) -> Optional[requests.Response]:
        """Return a synthetic 403 if the domain of *url* is suppressed.

        Args:
            url: URL to check.

        Returns:
            A synthetic :class:`requests.Response` with status 403 if the
            domain is suppressed, otherwise ``None``.
        """
        domain = self._domain_key(url)
        if domain not in self._suppressed_domains:
            return None
        if self.verbosity >= 2:
            reason = self._suppressed_domains[domain]
            print(f"  [HttpClient] Skipping suppressed domain: {domain} ({reason})")
        response = requests.Response()
        response.status_code = 403
        response.url = url
        response._content = b""
        return response

    def _register_suppression(self, url: str) -> None:
        """Record the domain of *url* as suppressed after a 403 response.

        Args:
            url: URL whose domain should be suppressed.
        """
        domain = self._domain_key(url)
        if domain not in self._suppressed_domains:
            reason = f"403 Forbidden at {time.strftime('%Y-%m-%d %H:%M:%S')}"
            self._suppressed_domains[domain] = reason
            if self.verbosity >= 1:
                print(f"  [HttpClient] Domain suppressed due to 403 Forbidden: {domain}")

    @property
    def suppressed_domains(self) -> Dict[str, str]:
        """Return a copy of the currently suppressed domains and reasons."""
        return dict(self._suppressed_domains)

    @property
    def _delay_enabled(self) -> bool:
        """``True`` when an explicit inter-request delay has been configured."""
        return self.delay_min_s is not None

    def _increase_backoff(self, domain: str) -> float:
        """Increase the exponential backoff delay for *domain*.

        On first call for a domain the delay is set to *backoff_base_s*.
        On subsequent calls it is multiplied by *backoff_multiplier*, capped
        at *backoff_max_s*.

        Returns:
            The new backoff delay in seconds.
        """
        current = self._backoff_delays.get(domain)
        if current is None:
            new_delay = self.backoff_base_s
        else:
            new_delay = current * self.backoff_multiplier
        new_delay = min(new_delay, self.backoff_max_s)
        self._backoff_delays[domain] = new_delay
        if self.verbosity >= 2:
            print(f"  [HttpClient] Backoff for {domain}: {new_delay:.1f}s")
        return new_delay

    def _decrease_backoff(self, domain: str) -> None:
        """Linearly decrease the backoff delay for *domain* on success.

        Subtracts *backoff_decay_s* from the current backoff.  When the
        result reaches zero or below the backoff is cleared entirely.
        """
        current = self._backoff_delays.get(domain)
        if current is None:
            return
        new_delay = current - self.backoff_decay_s
        if new_delay <= 0:
            del self._backoff_delays[domain]
            if self.verbosity >= 3:
                print(f"  [HttpClient] Backoff cleared for {domain}")
        else:
            self._backoff_delays[domain] = new_delay
            if self.verbosity >= 3:
                print(
                    f"  [HttpClient] Backoff decreased for {domain}: "
                    f"{new_delay:.1f}s"
                )

    @property
    def backoff_delays(self) -> Dict[str, float]:
        """Return a copy of the current per-domain backoff delays."""
        return dict(self._backoff_delays)

    def _apply_rate_limit(self, url: str) -> None:
        """Sleep to honour Crawl-Delay, backoff, explicit delay, or fallback.

        Priority:

        1. If the user configured *delay_min_s* / *delay_max_s*, those form
           the base range.
        2. A robots.txt Crawl-Delay **overrides the lower bound** — if it is
           larger than the current lower bound it raises both the floor and,
           if necessary, the ceiling.
        3. An active exponential backoff **overrides the lower bound** using
           the same logic as Crawl-Delay.
        4. When *respect_robots_txt* is ``True`` but neither an explicit delay
           nor a Crawl-Delay is present, the *rate_limit_min_s* /
           *rate_limit_max_s* fallback is used.
        """
        domain = self._domain_key(url)
        last_time = self._last_fetch_times.get(domain)
        if last_time is None:
            return

        lower: Optional[float] = None
        upper: Optional[float] = None

        # 1. Start with user-configured bounds
        if self._delay_enabled:
            lower = self.delay_min_s
            upper = self.delay_max_s

        # 2. Crawl-Delay overrides the lower bound
        if self.respect_robots_txt and self._robots_cache is not None:
            crawl_delay = self._robots_cache.crawl_delay(url)
            if crawl_delay is not None:
                if lower is None:
                    # No explicit delay — use Crawl-Delay as a fixed value
                    lower = crawl_delay
                    upper = crawl_delay
                else:
                    lower = max(lower, crawl_delay)
                    upper = max(upper, lower)  # type: ignore[arg-type]

        # 3. Backoff overrides the lower bound
        backoff = self._backoff_delays.get(domain)
        if backoff is not None:
            if lower is None:
                lower = backoff
                upper = backoff
            else:
                lower = max(lower, backoff)
                upper = max(upper, lower)  # type: ignore[arg-type]

        # 4. Fallback for robots.txt mode
        if lower is None and self.respect_robots_txt:
            lower = self.rate_limit_min_s
            upper = self.rate_limit_max_s

        if lower is None or upper is None:
            return

        delay_seconds = random.uniform(lower, upper)

        if self.verbosity >= 3:
            print(
                f"  [HttpClient] Delay for {domain}: "
                f"{delay_seconds:.2f}s (range {lower:.2f}–{upper:.2f}s)"
            )

        elapsed = time.time() - last_time
        sleep_time = delay_seconds - elapsed
        if sleep_time > 0:
            if self.verbosity >= 3:
                print(f"  [HttpClient] Sleeping {sleep_time:.2f}s for {domain}")
            time.sleep(sleep_time)

    def _handle_429(
        self, response: requests.Response, url: str
    ) -> tuple[float, bool]:
        """Parse a 429 response to determine how long to wait before retrying.

        The method checks, in order:

        1. ``Retry-After`` header (RFC 7231) — integer seconds or HTTP date.
        2. ``RateLimit-Reset`` / ``X-RateLimit-Reset`` (IETF draft) — Unix
           timestamp (>1 000 000 000) or seconds-to-reset.
        3. Falls back to *default_retry_after_s*.

        The returned value is clamped to *max_retry_after_s*.

        Args:
            response: The 429 response.
            url: The URL that was requested (for logging).

        Returns:
            A ``(wait_seconds, header_found)`` tuple.  *header_found* is
            ``True`` when a rate-limit header was present, ``False`` when the
            default fallback was used.
        """
        wait_seconds: Optional[float] = None
        header_used: Optional[str] = None
        header_found = False

        # 1. Retry-After (RFC 7231)
        retry_after = response.headers.get("Retry-After")
        if retry_after:
            try:
                wait_seconds = float(retry_after)
                header_used = "Retry-After"
                header_found = True
            except ValueError:
                try:
                    retry_time = parsedate_to_datetime(retry_after)
                    wait_seconds = (
                        retry_time - datetime.now(retry_time.tzinfo)
                    ).total_seconds()
                    wait_seconds = max(0.0, wait_seconds)
                    header_used = "Retry-After"
                    header_found = True
                except (ValueError, TypeError):
                    pass

        # 2. RateLimit-Reset / X-RateLimit-Reset (IETF draft)
        if wait_seconds is None:
            reset_header = (
                response.headers.get("RateLimit-Reset")
                or response.headers.get("X-RateLimit-Reset")
            )
            if reset_header:
                try:
                    reset_value = int(reset_header)
                    if reset_value > 1_000_000_000:
                        wait_seconds = max(0.0, reset_value - time.time())
                    else:
                        wait_seconds = float(reset_value)
                    header_used = (
                        "RateLimit-Reset"
                        if "RateLimit-Reset" in response.headers
                        else "X-RateLimit-Reset"
                    )
                    header_found = True
                except (ValueError, TypeError):
                    pass

        # 3. Default fallback
        if wait_seconds is None:
            wait_seconds = self.default_retry_after_s
            header_used = "default"
            if self.verbosity >= 1:
                print(
                    "  [HttpClient] No rate-limit headers found, "
                    f"using {self.default_retry_after_s}s default"
                )
                if self.verbosity >= 2:
                    print("  All response headers:")
                    for name, value in response.headers.items():
                        print(f"    {name}: {value}")

        # Clamp to max
        wait_seconds = min(wait_seconds, self.max_retry_after_s)

        if self.verbosity >= 1:
            print(
                f"  [HttpClient] HTTP 429 — waiting {wait_seconds:.0f}s "
                f"before retry (from {header_used})"
            )

        return wait_seconds, header_found

    #: Domains that act as redirectors and should never be suppressed, because
    #: a 403 from these domains originates at the *destination* site.
    SUPPRESSION_EXEMPT_DOMAINS: set[str] = {
        "doi.org",
        "dx.doi.org",
    }

    def _is_suppression_exempt(self, url: str) -> bool:
        """Return ``True`` if the domain of *url* must never be suppressed."""
        netloc = urlparse(url).netloc.lower()
        return netloc in self.SUPPRESSION_EXEMPT_DOMAINS

    def _request(
        self, method: str, url: str, retry_on_504: bool, **kwargs: Any
    ) -> Optional[requests.Response]:
        """Send one request with optional robots.txt, backoff, retry, and 403.

        Shared by :meth:`get` and :meth:`post`; both use the same
        per-domain clock, backoff, and suppression list.

        Args:
            method: ``"GET"`` or ``"POST"``.
            url: The URL to request.
            retry_on_504: Whether a 504 is retried.  False for requests
                that are not idempotent: after a gateway timeout the
                origin may already have acted on the request.
            **kwargs: Passed to :func:`requests.request`.

        Returns:
            A :class:`requests.Response`, or ``None`` if the URL is
            disallowed by robots.txt.
        """
        if not self.polite_mode:
            headers = kwargs.get("headers", {})
            if "User-Agent" not in headers:
                headers["User-Agent"] = self.user_agent
                kwargs["headers"] = headers
            return self._send(method, url, **kwargs)

        # 1. Suppression check (before any network I/O)
        if self.suppress_on_403:
            suppressed = self._check_suppression(url)
            if suppressed is not None:
                return suppressed

        # 2. robots.txt check
        if self.respect_robots_txt and self._robots_cache is not None:
            if not self._robots_cache.can_fetch(url):
                if self.verbosity >= 1:
                    print(f"  [HttpClient] Blocked by robots.txt: {url}")
                return None

        self._apply_rate_limit(url)

        headers = kwargs.get("headers", {})
        if "User-Agent" not in headers:
            headers["User-Agent"] = self.user_agent
            kwargs["headers"] = headers

        domain = self._domain_key(url)
        self._last_fetch_times[domain] = time.time()

        response = self._send(method, url, **kwargs)

        # 3. Retry loop for 429, 503, 504
        for _attempt in range(self.max_retries):
            if response.status_code == 429 and self.honor_retry_after:
                wait, header_found = self._handle_429(response, url)
                if not header_found:
                    backoff_wait = self._increase_backoff(domain)
                    wait = max(wait, backoff_wait)
                time.sleep(wait)
            elif response.status_code == 503 or (
                response.status_code == 504 and retry_on_504
            ):
                wait = self._increase_backoff(domain)
                if self.verbosity >= 1:
                    print(
                        f"  [HttpClient] HTTP {response.status_code} — "
                        f"backing off {wait:.0f}s before retry"
                    )
                time.sleep(wait)
            else:
                break

            self._last_fetch_times[domain] = time.time()
            response = self._send(method, url, **kwargs)

        # 4. Decrease backoff on success
        if 200 <= response.status_code < 300:
            self._decrease_backoff(domain)

        # 5. 403 suppression registration
        if (
            self.suppress_on_403
            and response.status_code == 403
            and not self._is_suppression_exempt(url)
        ):
            self._register_suppression(url)

        return response

    @staticmethod
    def _send(method: str, url: str, **kwargs: Any) -> requests.Response:
        """Dispatch to :func:`requests.get` or :func:`requests.post`."""
        if method == "POST":
            return requests.post(url, **kwargs)
        return requests.get(url, **kwargs)

    def get(self, url: str, **kwargs: Any) -> Optional[requests.Response]:
        """Make a GET request with optional robots.txt, backoff, retry, and 403.

        Args:
            url: The URL to fetch.
            **kwargs: Additional keyword arguments passed to
                :func:`requests.get`.

        Returns:
            A :class:`requests.Response`, or ``None`` if the URL is disallowed
            by robots.txt.
        """
        return self._request("GET", url, retry_on_504=True, **kwargs)

    def post(self, url: str, **kwargs: Any) -> Optional[requests.Response]:
        """Make a POST request under the same politeness rules as :meth:`get`.

        Robots.txt, 403 suppression, the per-domain delay (one clock per
        domain shared with :meth:`get`) and backoff all apply.  Because
        POST is not idempotent it is retried only on 429 and 503, which
        mean the server did not act on the request, and never on 504.

        Args:
            url: The URL to post to.
            **kwargs: Additional keyword arguments passed to
                :func:`requests.post` (``json=``, ``data=``, ``timeout=``...).

        Returns:
            A :class:`requests.Response`, or ``None`` if the URL is disallowed
            by robots.txt.
        """
        return self._request("POST", url, retry_on_504=False, **kwargs)
