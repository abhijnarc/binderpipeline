#!/usr/bin/env python3
"""Validate every run manifest under runs/. Run this before every render.

Checks, per run:

* every required manifest field is present, of the right type, non-empty
  where emptiness would be meaningless;
* every file listed in `outputs` actually exists in the run directory;
* every file listed in `inputs` still hashes to the recorded sha256, so a
  result computed from a since-changed dataset is visible as such;
* `description` is not missing or empty;
* `supersedes` names a run that exists, and no run is superseded twice;
* `run_id` matches the directory name and is unique;
* no oversized or binary artefact has crept into runs/.

Exit status is 0 when there are no errors, 1 otherwise. Stale inputs and
other interpretive problems are warnings; --strict turns them into errors.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent
RUNS_DIR = ROOT / "runs"

TRACKS = ("interface-kd", "physchem", "shared")
STATUSES = ("ok", "failed", "partial")

REQUIRED_FIELDS = (
    "run_id",
    "timestamp_utc",
    "track",
    "stage",
    "description",
    "status",
    "inputs",
    "outputs",
    "tool_versions",
    "code_version",
    "parameters",
    "supersedes",
    "backfilled",
    "notes",
)

MAX_RUN_FILE_BYTES = 5 * 1024 * 1024
BINARY_SUFFIXES = {
    ".pdb", ".cif", ".mmcif", ".ent", ".pt", ".pth", ".npy", ".npz",
    ".h5", ".hdf5", ".pkl", ".pickle", ".bam", ".sqlite", ".safetensors",
}


class Report:
    """Collected findings for the whole runs/ tree."""

    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def error(self, run_id: str, message: str) -> None:
        self.errors.append(f"{run_id}: {message}")

    def warn(self, run_id: str, message: str) -> None:
        self.warnings.append(f"{run_id}: {message}")


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_fields(run_id: str, manifest: dict, report: Report) -> None:
    """Required fields, with the types and emptiness rules that matter."""
    for field in REQUIRED_FIELDS:
        if field not in manifest:
            report.error(run_id, f"manifest is missing required field {field!r}")

    description = manifest.get("description")
    if not isinstance(description, str) or not description.strip():
        report.error(
            run_id,
            "description is missing or empty; a reader cannot tell what this "
            "run was",
        )

    track = manifest.get("track")
    if track is not None and track not in TRACKS:
        report.error(
            run_id, f"track {track!r} is not one of {', '.join(TRACKS)}"
        )

    status = manifest.get("status")
    if status is not None and status not in STATUSES:
        report.error(
            run_id, f"status {status!r} is not one of {', '.join(STATUSES)}"
        )

    if "backfilled" in manifest and not isinstance(
        manifest["backfilled"], bool
    ):
        report.error(
            run_id,
            f"backfilled must be true or false, got "
            f"{manifest['backfilled']!r}",
        )

    for field in ("inputs", "outputs"):
        value = manifest.get(field)
        if field in manifest and not isinstance(value, list):
            report.error(run_id, f"{field} must be a list")

    for field in ("tool_versions", "parameters"):
        value = manifest.get(field)
        if field in manifest and not isinstance(value, dict):
            report.error(run_id, f"{field} must be a mapping")

    if isinstance(manifest.get("tool_versions"), dict) and not manifest[
        "tool_versions"
    ]:
        report.warn(run_id, "tool_versions is empty")

    if isinstance(manifest.get("parameters"), dict) and not manifest[
        "parameters"
    ]:
        report.warn(
            run_id,
            "parameters is empty; record the actual values used, not a "
            "reference to defaults",
        )

    if manifest.get("code_version") is None:
        report.warn(
            run_id,
            "code_version is null; the code that produced this cannot be "
            "identified",
        )
    if manifest.get("code_version_dirty"):
        report.warn(
            run_id,
            "working tree was dirty when this ran, so the recorded SHA does "
            "not fully describe the code",
        )

    if manifest.get("backfilled") is True:
        report.warn(
            run_id,
            "backfilled: numbers were transcribed from a previous report, "
            "not produced by executing code",
        )

    if manifest.get("status") in ("failed", "partial"):
        report.warn(
            run_id, f"status is {manifest['status']!r}"
        )


def check_outputs(
    run_id: str, run_dir: Path, manifest: dict, report: Report
) -> None:
    outputs = manifest.get("outputs")
    if not isinstance(outputs, list):
        return
    if not outputs:
        report.warn(run_id, "no outputs listed")
    for index, entry in enumerate(outputs):
        if not isinstance(entry, dict) or "path" not in entry:
            report.error(
                run_id, f"outputs[{index}] must be a mapping with a 'path'"
            )
            continue
        target = run_dir / entry["path"]
        if not target.is_file():
            report.error(
                run_id,
                f"listed output {entry['path']!r} does not exist in the run "
                "directory",
            )
        if not str(entry.get("description", "")).strip():
            report.warn(
                run_id, f"output {entry['path']!r} has no description"
            )


def check_inputs(run_id: str, manifest: dict, report: Report) -> int:
    """Re-hash every recorded input. Returns the number that have changed."""
    inputs = manifest.get("inputs")
    if not isinstance(inputs, list):
        return 0
    changed = 0
    for index, entry in enumerate(inputs):
        if not isinstance(entry, dict) or "path" not in entry:
            report.error(
                run_id, f"inputs[{index}] must be a mapping with a 'path'"
            )
            continue
        recorded = entry.get("sha256")
        if not recorded:
            report.error(
                run_id, f"input {entry['path']!r} has no recorded sha256"
            )
            continue

        source = ROOT / entry["path"]
        if not source.is_file():
            report.warn(
                run_id,
                f"input {entry['path']!r} no longer exists, so its hash "
                "cannot be re-checked",
            )
            continue

        actual = sha256_of(source)
        if actual != recorded:
            changed += 1
            report.warn(
                run_id,
                f"STALE INPUT: {entry['path']!r} has changed since this run "
                f"(recorded {recorded[:12]}..., now {actual[:12]}...). This "
                "result was computed from a different version of that file.",
            )
    return changed


def check_external(run_id: str, manifest: dict, report: Report) -> None:
    for index, entry in enumerate(manifest.get("external_artefacts", []) or []):
        if not isinstance(entry, dict) or "path" not in entry:
            report.error(
                run_id,
                f"external_artefacts[{index}] must be a mapping with a 'path'",
            )
            continue
        if not entry.get("sha256"):
            report.error(
                run_id,
                f"external artefact {entry['path']!r} has no recorded sha256",
            )
        source = ROOT / entry["path"]
        if not source.is_file():
            report.warn(
                run_id,
                f"external artefact {entry['path']!r} is not present at the "
                "recorded location",
            )


def check_directory_size(run_id: str, run_dir: Path, report: Report) -> None:
    """runs/ must stay small, so the site opens quickly."""
    for path in sorted(run_dir.rglob("*")):
        if not path.is_file():
            continue
        size = path.stat().st_size
        relative = path.relative_to(run_dir)
        if size > MAX_RUN_FILE_BYTES:
            report.error(
                run_id,
                f"{relative} is {size} bytes, over the "
                f"{MAX_RUN_FILE_BYTES} byte limit; large artefacts belong "
                "outside runs/, recorded by path and hash",
            )
        if path.suffix.lower() in BINARY_SUFFIXES:
            report.error(
                run_id,
                f"{relative} looks like a binary artefact; record it in the "
                "manifest and leave the file where it lives",
            )


def load_manifest(run_dir: Path, report: Report) -> dict | None:
    manifest_path = run_dir / "manifest.yml"
    if not manifest_path.is_file():
        report.error(run_dir.name, "no manifest.yml")
        return None
    try:
        with manifest_path.open() as handle:
            manifest = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        report.error(run_dir.name, f"manifest.yml is not valid YAML: {exc}")
        return None
    if not isinstance(manifest, dict):
        report.error(run_dir.name, "manifest.yml is not a mapping")
        return None
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate every run manifest under runs/."
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="treat warnings (stale inputs, missing parameters) as errors",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="print only problems and the summary",
    )
    args = parser.parse_args(argv)

    if not RUNS_DIR.is_dir():
        print(f"no runs directory at {RUNS_DIR}; nothing to validate")
        return 0

    report = Report()
    run_dirs = sorted(d for d in RUNS_DIR.iterdir() if d.is_dir())
    manifests: dict[str, dict] = {}
    stale_total = 0

    for run_dir in run_dirs:
        manifest = load_manifest(run_dir, report)
        if manifest is None:
            continue

        run_id = str(manifest.get("run_id") or run_dir.name)
        if manifest.get("run_id") != run_dir.name:
            report.error(
                run_dir.name,
                f"run_id {manifest.get('run_id')!r} does not match the "
                "directory name",
            )
        if run_id in manifests:
            report.error(run_id, "duplicate run_id")
        manifests[run_id] = manifest

        check_fields(run_id, manifest, report)
        check_outputs(run_id, run_dir, manifest, report)
        stale_total += check_inputs(run_id, manifest, report)
        check_external(run_id, manifest, report)
        check_directory_size(run_id, run_dir, report)

        if not args.quiet:
            flags = []
            if manifest.get("backfilled"):
                flags.append("backfilled")
            if manifest.get("supersedes"):
                flags.append(f"supersedes {manifest['supersedes']}")
            suffix = f"  [{'; '.join(flags)}]" if flags else ""
            print(
                f"{run_id}  track={manifest.get('track')} "
                f"stage={manifest.get('stage')} "
                f"status={manifest.get('status')}{suffix}"
            )

    # Supersession is checked across runs, once all manifests are loaded.
    superseded_by: dict[str, list[str]] = {}
    for run_id, manifest in manifests.items():
        target = manifest.get("supersedes")
        if target is None:
            continue
        if not isinstance(target, str):
            report.error(run_id, f"supersedes must be a run_id or null")
            continue
        if target == run_id:
            report.error(run_id, "supersedes itself")
        if target not in manifests:
            report.error(
                run_id,
                f"supersedes {target!r}, which is not a run under runs/",
            )
        superseded_by.setdefault(target, []).append(run_id)

    for target, replacements in sorted(superseded_by.items()):
        if len(replacements) > 1:
            report.warn(
                target,
                "superseded by more than one run "
                f"({', '.join(sorted(replacements))}); which one is current "
                "is ambiguous",
            )

    print()
    if report.warnings:
        print(f"{len(report.warnings)} warning(s):")
        for message in report.warnings:
            print(f"  WARN  {message}")
        print()
    if report.errors:
        print(f"{len(report.errors)} error(s):")
        for message in report.errors:
            print(f"  ERROR {message}")
        print()

    print(
        f"checked {len(manifests)} run(s); "
        f"{len(report.errors)} error(s), {len(report.warnings)} warning(s)"
    )
    if stale_total:
        print(
            f"{stale_total} recorded input(s) have changed since their run — "
            "those results came from a different version of the data"
        )

    if report.errors:
        print("FAILED")
        return 1
    if args.strict and report.warnings:
        print("FAILED (--strict: warnings count as errors)")
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
