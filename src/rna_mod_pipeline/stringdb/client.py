from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Mapping


class StringAPIError(RuntimeError):
    pass


def _normalise(parameters: Mapping[str, object]) -> dict[str, object]:
    return {str(key): parameters[key] for key in sorted(parameters)}


class CachedStringClient:
    def __init__(
        self,
        cache_dir: str | Path,
        *,
        offline: bool = False,
        retries: int = 3,
        delay_seconds: float = 1.0,
        base_url: str = "https://string-db.org/api",
    ):
        self.cache_dir = Path(cache_dir).resolve()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.retries = retries
        self.delay_seconds = delay_seconds
        self.base_url = base_url.rstrip("/")
        self.request_log: list[dict[str, object]] = []
        self._last_network_call = 0.0

    def cache_key(
        self, endpoint: str, response_format: str, parameters: Mapping[str, object]
    ) -> str:
        request = {
            "endpoint": endpoint,
            "format": response_format,
            "parameters": _normalise(parameters),
        }
        raw = json.dumps(request, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(raw).hexdigest()

    def call(
        self,
        endpoint: str,
        parameters: Mapping[str, object],
        *,
        response_format: str = "json",
    ):
        params = _normalise(parameters)
        key = self.cache_key(endpoint, response_format, params)
        extension = ".json" if response_format == "json" else f".{response_format}"
        cache_path = self.cache_dir / f"{key}{extension}"
        request_path = self.cache_dir / f"{key}.request.json"
        source = "cache"
        if cache_path.is_file():
            raw = cache_path.read_bytes()
        else:
            if self.offline:
                raise StringAPIError(f"Offline STRING cache miss: {endpoint} ({key})")
            encoded = urllib.parse.urlencode(params).encode()
            url = f"{self.base_url}/{response_format}/{endpoint}"
            error = None
            for attempt in range(1, self.retries + 1):
                elapsed = time.monotonic() - self._last_network_call
                if elapsed < self.delay_seconds:
                    time.sleep(self.delay_seconds - elapsed)
                try:
                    request = urllib.request.Request(url, data=encoded, method="POST")
                    with urllib.request.urlopen(request, timeout=90) as response:
                        raw = response.read()
                    self._last_network_call = time.monotonic()
                    break
                except (urllib.error.URLError, TimeoutError, OSError) as exc:
                    error = exc
                    self._last_network_call = time.monotonic()
                    if attempt < self.retries:
                        time.sleep(min(2**attempt, 10))
            else:
                raise StringAPIError(f"STRING request failed: {endpoint}: {error}")
            temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
            temporary.write_bytes(raw)
            temporary.replace(cache_path)
            request_path.write_text(
                json.dumps(
                    {
                        "endpoint": endpoint,
                        "format": response_format,
                        "parameters": params,
                        "cache_key": key,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            source = "network"
        digest = hashlib.sha256(raw).hexdigest()
        self.request_log.append(
            {
                "endpoint": endpoint,
                "response_format": response_format,
                "cache_key": key,
                "response_source": source,
                "response_bytes": len(raw),
                "response_sha256": digest,
                "cache_file": cache_path.name,
            }
        )
        if response_format == "json":
            try:
                return json.loads(raw)
            except json.JSONDecodeError as exc:
                raise StringAPIError(f"Invalid JSON from STRING endpoint {endpoint}") from exc
        return raw
