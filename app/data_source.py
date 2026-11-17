"""Where the app's data comes from: local files natively, HTTP fetches under Shinylive.

Under Shinylive the app runs in a web worker at ``<site>/shinylive/pyodide-worker.js``; the
data files are served from ``<site>/data/...`` (``scripts/build_site.py`` copies them there).
"""

from __future__ import annotations

import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
IN_BROWSER = sys.platform == "emscripten"
_CACHE: dict[str, bytes] = {}


def _repo_root() -> Path:
    for p in (APP_DIR, APP_DIR.parent):
        if (p / "data").is_dir():
            return p
    return APP_DIR


def site_base() -> str:
    """Site root URL under Shinylive."""
    import js  # type: ignore[import-not-found]

    href = str(js.location.href)
    return href.split("/shinylive/")[0] if "/shinylive/" in href else href.rsplit("/", 1)[0]


async def read_bytes(rel: str) -> bytes:
    """Contents of ``<root>/<rel>`` (e.g. ``data/precomputed/manifest.json``), cached."""
    if rel in _CACHE:
        return _CACHE[rel]
    if IN_BROWSER:
        from pyodide.http import pyfetch  # type: ignore[import-not-found]

        resp = await pyfetch(f"{site_base()}/{rel}")
        if not resp.ok:
            raise OSError(f"could not fetch {rel}: HTTP {resp.status}")
        data = await resp.bytes()
    else:
        data = (_repo_root() / rel).read_bytes()
    _CACHE[rel] = data
    return data


def histone_dir() -> Path:
    """Histone inputs: bundled with the app (Shinylive) or the repository's data/histones."""
    bundled = APP_DIR / "data" / "histones"
    return bundled if bundled.is_dir() else _repo_root() / "data" / "histones"
