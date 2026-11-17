"""Build the static Shinylive site in site/ (ready for GitHub Pages).

Steps: bundle app/ + the core package + the histone inputs into build/app, run
``shinylive export``, copy data/precomputed into site/data/precomputed, add the crawler
blocks (a ``noindex`` meta tag in every HTML page and a ``robots.txt``), and check that no
file exceeds the 25 MB limit (CLAUDE.md section 2).

Usage:
    uv run python scripts/build_site.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build" / "app"
SITE = ROOT / "site"
LIMIT = 25 * 1024 * 1024
NOINDEX = '<meta name="robots" content="noindex" />'
ROBOTS_TXT = "User-agent: *\nDisallow: /\n"
HISTONE_FILES = ["human_histones.fasta", "histone_isoforms.fasta", "histone_annotation.tsv",
                 "known_marks.tsv", "metadata.json"]  # fmt: skip


def main() -> None:
    shutil.rmtree(BUILD, ignore_errors=True)
    shutil.rmtree(SITE, ignore_errors=True)
    BUILD.mkdir(parents=True)
    for f in (ROOT / "app").iterdir():
        if f.suffix in (".py", ".txt") and not f.name.startswith("test_"):
            shutil.copy(f, BUILD / f.name)
    shutil.copytree(ROOT / "src" / "proteasescope", BUILD / "proteasescope",
                    ignore=shutil.ignore_patterns("__pycache__"))  # fmt: skip
    (BUILD / "data" / "histones").mkdir(parents=True)
    for name in HISTONE_FILES:
        shutil.copy(ROOT / "data" / "histones" / name, BUILD / "data" / "histones" / name)
    shinylive = Path(sys.executable).with_name("shinylive")
    subprocess.run([str(shinylive), "export", str(BUILD), str(SITE)],
                   check=True)  # fmt: skip
    shutil.copytree(ROOT / "data" / "precomputed", SITE / "data" / "precomputed")
    (SITE / ".nojekyll").write_text("")  # GitHub Pages: serve files/dirs starting with _
    block_crawlers(SITE)
    big = [(p, p.stat().st_size) for p in SITE.rglob("*") if p.is_file()
           and p.stat().st_size > LIMIT]  # fmt: skip
    total = sum(p.stat().st_size for p in SITE.rglob("*") if p.is_file())
    print(f"\nsite: {SITE}  total {total / 1e6:.1f} MB")
    for p, size in big:
        print(f"  OVER 25 MB: {p.relative_to(SITE)} ({size / 1e6:.1f} MB)")
    if big:
        sys.exit(1)


def block_crawlers(site: Path) -> None:
    """Add a noindex meta tag to every HTML page and write a blocking robots.txt.

    On a GitHub Pages project site, crawlers only read robots.txt at the domain root, so the
    meta tag is what keeps the pages out of search results; robots.txt is kept for hosts
    that serve the site at a domain root.
    """
    pages = list(site.rglob("*.html"))
    for page in pages:
        html = page.read_text()
        if NOINDEX not in html:
            if "<head>" not in html:
                sys.exit(f"no <head> in {page}; cannot add the noindex tag")
            html = html.replace("<head>", f"<head>\n    {NOINDEX}", 1)
            page.write_text(html)
    (site / "robots.txt").write_text(ROBOTS_TXT)
    print(f"noindex added to {len(pages)} HTML page(s); robots.txt written")


if __name__ == "__main__":
    main()
