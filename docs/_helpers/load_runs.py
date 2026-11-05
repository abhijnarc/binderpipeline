"""Read the run manifests under runs/ and expose them to the Quarto pages.

Every number on this site is read from a file: either from a run's output
files, or from the `parameters` block of its manifest. Nothing is typed into
a page by hand. These helpers are the only way the pages reach the data.

Usage from a .qmd page:

    import sys, pathlib
    for candidate in [pathlib.Path.cwd(), *pathlib.Path.cwd().parents]:
        helpers = candidate / "_helpers"
        if (helpers / "load_runs.py").exists():
            sys.path.insert(0, str(helpers))
            break
    from load_runs import load_runs, params, runs_table
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path

import pandas as pd
import yaml

# ---------------------------------------------------------------------------
# Locating the repository
# ---------------------------------------------------------------------------


def find_root(start: Path | None = None) -> Path:
    """Walk up from this file (or `start`) to the directory holding runs/."""
    here = (start or Path(__file__).resolve()).resolve()
    for candidate in [here, *here.parents]:
        if (candidate / "runs").is_dir() and (candidate / "docs").is_dir():
            return candidate
    raise RuntimeError(
        "could not locate the repository root (a directory containing both "
        "runs/ and docs/)"
    )


ROOT = find_root()
RUNS_DIR = ROOT / "runs"

TRACK_LABELS = {
    "interface-kd": "Interface parameters and Kd",
    "physchem": "Physicochemical characterisation",
    "shared": "Shared",
}

STATUS_LABELS = {"ok": "ok", "partial": "partial", "failed": "failed"}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _read_manifest(run_dir: Path) -> dict | None:
    path = run_dir / "manifest.yml"
    if not path.is_file():
        return None
    with path.open() as handle:
        loaded = yaml.safe_load(handle)
    return loaded if isinstance(loaded, dict) else None


def load_manifests() -> dict[str, dict]:
    """Every manifest, keyed by run_id."""
    manifests: dict[str, dict] = {}
    if not RUNS_DIR.is_dir():
        return manifests
    for run_dir in sorted(RUNS_DIR.iterdir()):
        if not run_dir.is_dir():
            continue
        manifest = _read_manifest(run_dir)
        if manifest is None:
            continue
        manifests[str(manifest.get("run_id", run_dir.name))] = manifest
    return manifests


def load_runs() -> pd.DataFrame:
    """One row per run, newest first, with supersession resolved.

    Columns: run_id, timestamp_utc, date, track, track_label, stage, status,
    description, backfilled, code_version, supersedes, superseded_by,
    is_superseded, n_inputs, n_outputs, notes, path.
    """
    manifests = load_manifests()

    superseded_by: dict[str, list[str]] = {}
    for run_id, manifest in manifests.items():
        target = manifest.get("supersedes")
        if isinstance(target, str) and target:
            superseded_by.setdefault(target, []).append(run_id)

    rows = []
    for run_id, manifest in manifests.items():
        replacements = sorted(superseded_by.get(run_id, []))
        rows.append(
            {
                "run_id": run_id,
                "timestamp_utc": str(manifest.get("timestamp_utc", "")),
                "date": str(manifest.get("timestamp_utc", ""))[:10],
                "track": manifest.get("track", ""),
                "track_label": TRACK_LABELS.get(
                    manifest.get("track", ""), manifest.get("track", "")
                ),
                "stage": manifest.get("stage", ""),
                "status": manifest.get("status", ""),
                "description": str(manifest.get("description", "")).strip(),
                "backfilled": bool(manifest.get("backfilled", False)),
                "code_version": manifest.get("code_version"),
                "supersedes": manifest.get("supersedes"),
                "superseded_by": replacements[-1] if replacements else None,
                "is_superseded": bool(replacements),
                "n_inputs": len(manifest.get("inputs") or []),
                "n_outputs": len(manifest.get("outputs") or []),
                "notes": str(manifest.get("notes", "")).strip(),
                "path": f"runs/{run_id}",
            }
        )

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values(
        "timestamp_utc", ascending=False, kind="mergesort"
    ).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Reaching into a single run
# ---------------------------------------------------------------------------


def manifest(run_id: str) -> dict:
    """One run's manifest. Raises if the run_id is unknown."""
    manifests = load_manifests()
    if run_id not in manifests:
        raise KeyError(
            f"no run {run_id!r} under runs/. Known: "
            f"{', '.join(sorted(manifests)) or 'none'}"
        )
    return manifests[run_id]


def params(run_id: str) -> dict:
    """The `parameters` block: where a backfilled scalar result lives."""
    return dict(manifest(run_id).get("parameters") or {})


def param(run_id: str, key: str, default=None):
    """One recorded parameter value, read from the manifest file."""
    values = params(run_id)
    if key not in values:
        if default is not None:
            return default
        raise KeyError(
            f"run {run_id!r} records no parameter {key!r}. Recorded: "
            f"{', '.join(sorted(values)) or 'none'}"
        )
    return values[key]


def num(run_id: str, key: str) -> float:
    """A recorded parameter as a float.

    Manifests are immutable, and some older ones hold numbers written in
    exponent notation that YAML 1.1 parsed as strings ("1e-15"). Coercing on
    read keeps those runs usable without editing what was written.
    """
    value = param(run_id, key)
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(
            f"run {run_id!r} parameter {key!r} is {value!r}, which is not "
            "numeric"
        ) from None


def as_list(run_id: str, key: str) -> list[str]:
    """A recorded parameter as a list of strings.

    A multi-valued parameter may be stored either as a YAML list or as one
    comma-separated string, depending on how it was passed to new_run.py.
    Both are accepted so a page never has to guess.
    """
    value = param(run_id, key)
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value]
    return [item.strip() for item in str(value).split(",") if item.strip()]


def tool_versions(run_id: str) -> dict:
    return dict(manifest(run_id).get("tool_versions") or {})


def output_path(run_id: str, filename: str) -> Path:
    """Filesystem path to one of a run's output files."""
    path = RUNS_DIR / run_id / filename
    if not path.is_file():
        available = sorted(
            entry.get("path", "")
            for entry in (manifest(run_id).get("outputs") or [])
        )
        raise FileNotFoundError(
            f"run {run_id!r} has no output {filename!r}. Available: "
            f"{', '.join(available) or 'none'}"
        )
    return path


def read_output(run_id: str, filename: str, **kwargs) -> pd.DataFrame:
    """Read a run's CSV or TSV output into a DataFrame."""
    path = output_path(run_id, filename)
    separator = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    return pd.read_csv(path, sep=kwargs.pop("sep", separator), **kwargs)


def latest(track: str | None = None, stage: str | None = None) -> str | None:
    """run_id of the newest run matching track and stage, superseded ones
    excluded. Returns None when nothing matches."""
    frame = load_runs()
    if frame.empty:
        return None
    frame = frame[~frame.is_superseded]
    if track is not None:
        frame = frame[frame.track == track]
    if stage is not None:
        frame = frame[frame.stage == stage]
    if frame.empty:
        return None
    return str(frame.iloc[0].run_id)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def run_link(run_id: str, label: str | None = None, prefix: str = "") -> str:
    """Markdown link to a run directory, for citing a number's source.

    `prefix` is the path from the page to the site root ("" for a page at the
    root, "../" for one in a subdirectory).
    """
    return f"[{label or run_id}]({prefix}runs/{run_id}/)"


def cite(run_id: str, prefix: str = "") -> str:
    """Inline provenance marker: a linked run_id in brackets."""
    return f"([{run_id}]({prefix}runs/{run_id}/))"


def _escape(text: str) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def runs_table(frame: pd.DataFrame | None = None, prefix: str = "") -> str:
    """The full run table as Markdown, newest first.

    Superseded runs are struck through and name what replaced them, so the
    lineage of a corrected result stays visible rather than being
    overwritten.
    """
    frame = load_runs() if frame is None else frame
    if frame.empty:
        return "_No runs recorded yet._"

    lines = [
        "| Run | Date | Track | Stage | Status | Source | Description |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in frame.itertuples():
        link = run_link(row.run_id, prefix=prefix)
        if row.is_superseded:
            replacement = run_link(row.superseded_by, prefix=prefix)
            identifier = f"~~{link}~~<br/>replaced by {replacement}"
            description = f"<span class='superseded'>{_escape(row.description)}</span>"
        else:
            identifier = link
            description = _escape(row.description)

        status = row.status
        if status != "ok":
            status = f"**{status}**"

        source = "backfilled" if row.backfilled else "executed"
        if row.backfilled:
            source = "**backfilled**"

        lines.append(
            f"| {identifier} | {row.date} | {row.track} | {row.stage} "
            f"| {status} | {source} | {description} |"
        )
    return "\n".join(lines)


def status_table(track: str, stages: list[tuple[str, str, str]],
                 prefix: str = "") -> str:
    """Per-track status table: stage, state, last run, link.

    `stages` is a list of (stage_key, label, planned_state). A stage with a
    matching run reports that run's status and date; one without reports the
    planned state and no run.
    """
    lines = [
        "| Stage | State | Last run | Source |",
        "|---|---|---|---|",
    ]
    frame = load_runs()
    for stage_key, label, planned in stages:
        matching = (
            frame[(frame.track == track) & (frame.stage == stage_key)
                  & (~frame.is_superseded)]
            if not frame.empty
            else frame
        )
        if matching.empty:
            lines.append(f"| {label} | {planned} | — | — |")
            continue
        row = matching.iloc[0]
        state = row.status if row.status == "ok" else f"**{row.status}**"
        source = "**backfilled**" if row.backfilled else "executed"
        lines.append(
            f"| {label} | {state} | {run_link(row.run_id, row.date, prefix)} "
            f"| {source} |"
        )
    return "\n".join(lines)


def fmt_rho(rho: float, n: int, p: float | None = None,
            threshold: float = 0.05) -> str:
    """Format a rank correlation with its sample size and significance.

    Site convention: a coefficient never appears without n, and one that
    does not reach significance says so in the same breath.
    """
    text = f"rho = {float(rho):+.2f} (n = {n}"
    if p is not None:
        text += f", p = {float(p):.2f}"
        text += (
            ", not significant)"
            if float(p) >= threshold
            else f", p < {threshold})"
        )
    else:
        text += ")"
    return text


def today() -> str:
    return _dt.date.today().isoformat()


def counts_by(frame: pd.DataFrame, column: str) -> dict:
    if frame.empty:
        return {}
    return frame[column].value_counts().to_dict()
