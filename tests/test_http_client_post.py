"""Tests for HttpClient.post().

post() shares get()'s polite-mode machinery: 403 suppression, robots.txt,
the per-domain rate limiter (one clock per domain for both methods), and
backoff.  It differs in one place: POST is not idempotent, so it is
retried only on responses that mean the server did not process the
request (429, 503), never on 504, where the request may have run.
"""

from unittest.mock import Mock, patch

from pypaperretriever.http_client import HttpClient, RobotsTxtCache

URL = "https://example.com/search/query"


def _polite(can_fetch=True, crawl_delay=None, **kwargs):
    client = HttpClient(user_agent="TestBot/1.0", respect_robots_txt=True,
                        verbosity=0, **kwargs)
    cache = Mock(spec=RobotsTxtCache)
    cache.can_fetch.return_value = can_fetch
    cache.crawl_delay.return_value = crawl_delay
    client._robots_cache = cache
    return client


class TestPostPassthrough:

    @patch("pypaperretriever.http_client.requests.post")
    def test_forwards_body_and_sets_user_agent(self, mock_post):
        mock_post.return_value = Mock(status_code=200)
        client = HttpClient(user_agent="TestBot/1.0", polite_mode=False)
        resp = client.post(URL, json={"q": 1}, timeout=5)
        assert resp.status_code == 200
        args, kwargs = mock_post.call_args
        assert args[0] == URL
        assert kwargs["json"] == {"q": 1}
        assert kwargs["timeout"] == 5
        assert kwargs["headers"]["User-Agent"] == "TestBot/1.0"

    @patch("pypaperretriever.http_client.requests.post")
    def test_preserves_caller_headers(self, mock_post):
        mock_post.return_value = Mock(status_code=200)
        client = HttpClient(user_agent="TestBot/1.0", polite_mode=False)
        client.post(URL, headers={"User-Agent": "Mine/2.0"})
        assert mock_post.call_args[1]["headers"]["User-Agent"] == "Mine/2.0"


class TestPostPolite:

    @patch("pypaperretriever.http_client.requests.post")
    def test_returns_none_when_robots_disallow(self, mock_post):
        client = _polite(can_fetch=False)
        assert client.post(URL, json={}) is None
        mock_post.assert_not_called()

    @patch("pypaperretriever.http_client.requests.post")
    def test_suppressed_domain_short_circuits(self, mock_post):
        client = _polite()
        client._register_suppression(URL)
        resp = client.post(URL, json={})
        assert resp.status_code == 403
        mock_post.assert_not_called()

    @patch("pypaperretriever.http_client.requests.post")
    def test_403_registers_suppression(self, mock_post):
        mock_post.return_value = Mock(status_code=403, headers={})
        client = _polite()
        client.post(URL, json={})
        assert client._check_suppression(URL) is not None

    @patch("pypaperretriever.http_client.time.sleep")
    @patch("pypaperretriever.http_client.time.time")
    @patch("pypaperretriever.http_client.requests.get")
    @patch("pypaperretriever.http_client.requests.post")
    def test_shares_the_domain_clock_with_get(
        self, mock_post, mock_get, mock_time, mock_sleep
    ):
        mock_post.return_value = Mock(status_code=200, headers={})
        mock_get.return_value = Mock(status_code=200, headers={})
        mock_time.side_effect = [100.0, 100.5, 100.5]
        client = _polite(crawl_delay=2.0)

        client.post(URL, json={})
        client.get("https://example.com/pub/1/")

        mock_sleep.assert_called_once()
        assert abs(mock_sleep.call_args[0][0] - 1.5) < 0.01


class TestPostRetry:

    @patch("pypaperretriever.http_client.time.sleep")
    @patch("pypaperretriever.http_client.time.time", return_value=100.0)
    @patch("pypaperretriever.http_client.requests.post")
    def test_retries_on_429(self, mock_post, _mock_time, mock_sleep):
        mock_post.side_effect = [
            Mock(status_code=429, headers={"Retry-After": "2"}),
            Mock(status_code=200, headers={}),
        ]
        client = HttpClient(honor_retry_after=True, verbosity=0)
        resp = client.post(URL, json={})
        assert resp.status_code == 200
        assert mock_post.call_count == 2
        mock_sleep.assert_called_once_with(2.0)

    @patch("pypaperretriever.http_client.time.sleep")
    @patch("pypaperretriever.http_client.time.time", return_value=100.0)
    @patch("pypaperretriever.http_client.requests.post")
    def test_retries_on_503(self, mock_post, _mock_time, _mock_sleep):
        mock_post.side_effect = [
            Mock(status_code=503, headers={}),
            Mock(status_code=200, headers={}),
        ]
        client = HttpClient(verbosity=0)
        assert client.post(URL, json={}).status_code == 200
        assert mock_post.call_count == 2

    @patch("pypaperretriever.http_client.time.sleep")
    @patch("pypaperretriever.http_client.time.time", return_value=100.0)
    @patch("pypaperretriever.http_client.requests.post")
    def test_does_not_retry_on_504(self, mock_post, _mock_time, mock_sleep):
        mock_post.return_value = Mock(status_code=504, headers={})
        client = HttpClient(verbosity=0)
        assert client.post(URL, json={}).status_code == 504
        mock_post.assert_called_once()
        mock_sleep.assert_not_called()


class TestGetUnchanged:

    @patch("pypaperretriever.http_client.time.sleep")
    @patch("pypaperretriever.http_client.time.time", return_value=100.0)
    @patch("pypaperretriever.http_client.requests.get")
    def test_get_still_retries_on_504(self, mock_get, _mock_time, _mock_sleep):
        # Regression guard for the shared request path.
        mock_get.side_effect = [
            Mock(status_code=504, headers={}),
            Mock(status_code=200, headers={}),
        ]
        client = HttpClient(verbosity=0)
        assert client.get(URL).status_code == 200
        assert mock_get.call_count == 2
