# ProteaseScope

In silico protease digestion and comparison: what each alternative digest lets you see that a baseline digest cannot, and what it loses. It covers the whole reviewed human proteome (with isoforms) and a histone module (known marks, isoform discrimination, lysine co-occurrence).

> **All results are predictions** from cleavage rules. They are not observed data. Every table, figure and export carries the settings that produced it.

**Try it in your browser:** https://crogers14.github.io/ProteaseScope/ (beta). The app runs entirely in the browser with Shinylive; the first load downloads about 150 MB and can take a minute.

## Layout

| Path | What it is |
|---|---|
| `src/proteasescope/` | Core library (numpy + pandas only; no UI code) |
| `app/` | Shiny for Python app; runs natively or in the browser via Shinylive |
| `scripts/` | `fetch_data.py`, `build_known_marks.py`, `precompute.py`, `build_site.py`, exports |
| `data/` | Reference data (with `metadata.json` provenance), histone inputs, precomputed app data |
| `tests/` | pytest suite, including golden tests against Protein Cleaver |
| `bench/pyodide/` | Pyodide benchmarks, plus a runner that executes the test suite inside Pyodide |

Decisions and their reasons are recorded in `CHANGELOG.md`.

## Develop

```bash
uv sync --all-groups
```

```bash
uv run pytest
```

Run the app natively at http://localhost:8010:

```bash
uv run shiny run app/app.py --port 8010
```

## Rebuild data

Run each command when its input changes:

```bash
uv run python scripts/fetch_data.py --force
```

```bash
uv run python scripts/build_known_marks.py
```

```bash
uv run python scripts/precompute.py
```

`build_known_marks.py` marks `known_marks.tsv` as DRAFT; it must be reviewed again before use. `precompute.py` takes about 16 minutes and needs about 4.5 GB of memory.

## Build the static site (Shinylive)

```bash
uv run python scripts/build_site.py
```

This writes `site/`: the Shinylive app plus `data/precomputed/`, about 150 MB, with no file over 25 MB. Every page carries a `noindex` tag. Serve `site/` from any static host. To preview locally:

```bash
python3 -m http.server --directory site 8020
```

## Test inside Pyodide

The app runs on Pyodide 0.27.7 (bundled with Shinylive 0.10.16). This runs the core test suite on that runtime:

```bash
cd bench/pyodide && npm install && PYODIDE_PACKAGES=pytest node run.mjs run_tests.py ../..
```

## Deployment

Every push to `main` runs `.github/workflows/pages.yml`: ruff and pytest first, then `scripts/build_site.py`, then a deploy to GitHub Pages. Nothing is deployed if a test fails. The repository's Pages source must be set to **GitHub Actions** (Settings → Pages).
