#!/usr/bin/env python3
"""Scaffold an immutable run directory under runs/.

Every pipeline execution gets its own directory, named

    runs/<UTC-timestamp>_<track>_<short-slug>/

holding its output files plus a manifest.yml describing what was run, from
what inputs, with which tool versions and parameters. A run directory is
append-only: once written it is never edited. A corrected or re-executed
analysis becomes a NEW run directory naming the one it replaces, via
--supersedes. That is what makes the result history exist independently of
git.

Examples
--------
Scaffold a run, hashing the input it consumed and copying in its outputs:

    python3 new_run.py \
        --track physchem --stage smoke-test --slug smoketest \
        --description "Smoke test of the sequence-only descriptors." \
        --input physchem/data/opmd_binders_clean.csv \
        --copy physchem/out/physchem_properties.tsv:"Main property table" \
        --param ph_values=7.0,8.0 --param window=9

Record a large artefact without copying it into the site:

    --external bench20/work/1avx_0_0002.pdb:"Relaxed model, kept in place"

Replace an earlier run:

    --supersedes 2026-11-05T0840Z_physchem_smoketest
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent
RUNS_DIR = ROOT / "runs"

TRACKS = ("interface-kd", "physchem", "shared")
STATUSES = ("ok", "failed", "partial")

TIMESTAMP_FORMAT = "%Y-%m-%dT%H%MZ"
SLUG_PATTERN = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")

# runs/ holds text-sized results only. Anything bigger is recorded by path
# and hash via --external and left where it lives, so the site stays small.
MAX_COPY_BYTES = 5 * 1024 * 1024
BINARY_SUFFIXES = {
    ".pdb", ".cif", ".mmcif", ".ent", ".sc", ".pt", ".pth", ".npy", ".npz",
    ".h5", ".hdf5", ".pkl", ".pickle", ".bam", ".sqlite", ".zip", ".tar",
    ".gz", ".bz2", ".xz", ".safetensors",
}

# Libraries whose versions are worth recording when importable. Absent ones
# are simply omitted rather than recorded as a guess.
CANDIDATE_LIBRARIES = (
    "Bio", "peptides", "pandas", "numpy", "scipy", "sklearn", "yaml",
    "freesasa", "matplotlib", "torch", "esm",
)
LIBRARY_DISPLAY_NAMES = {
    "Bio": "biopython",
    "sklearn": "scikit-learn",
    "yaml": "pyyaml",
}


def sha256_and_size(path: Path) -> tuple[str, int]:
    """Hash a file in chunks so large artefacts do not have to be read whole."""
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def relative_to_root(path: Path) -> str:
    """Path as recorded in a manifest: repo-relative where possible."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def git_code_version() -> tuple[str | None, bool]:
    """Return (commit SHA, working-tree-dirty). SHA is None if unavailable.

    The manifest records null rather than a guess when the repo is absent,
    so a result is never attributed to code that cannot be identified.
    """
    def git(*args: str) -> str | None:
        try:
            done = subprocess.run(
                ("git", *args),
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        if done.returncode != 0:
            return None
        return done.stdout.strip()

    sha = git("rev-parse", "HEAD")
    if not sha:
        return None, False
    status = git("status", "--porcelain")
    return sha, bool(status)


def collect_tool_versions() -> dict[str, str]:
    """Versions of the interpreter and every importable library of interest."""
    import importlib
    import platform

    versions = {"python": platform.python_version()}
    for module_name in CANDIDATE_LIBRARIES:
        try:
            module = importlib.import_module(module_name)
        except Exception:
            continue
        version = getattr(module, "__version__", None)
        if version is None:
            try:
                from importlib.metadata import version as dist_version

                version = dist_version(
                    LIBRARY_DISPLAY_NAMES.get(module_name, module_name)
                )
            except Exception:
                continue
        display = LIBRARY_DISPLAY_NAMES.get(module_name, module_name)
        versions[display] = str(version)

    quarto = shutil.which("quarto")
    if quarto:
        try:
            done = subprocess.run(
                (quarto, "--version"), capture_output=True, text=True,
                timeout=30,
            )
            if done.returncode == 0:
                versions["quarto"] = done.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass

    return versions


def split_annotated(value: str, what: str) -> tuple[str, str]:
    """Split a "path:description" argument. The description may be empty."""
    if ":" in value:
        path_part, description = value.split(":", 1)
        return path_part, description.strip()
    print(
        f"note: no description given for {what} {value!r}; "
        "use PATH:\"description\" to add one",
        file=sys.stderr,
    )
    return value, ""


def parse_parameter(value: str) -> tuple[str, object]:
    """Parse KEY=VALUE, interpreting the value as YAML scalar where possible.

    So --param window=9 records the integer 9, --param ph=7.0,8.0 records the
    list [7.0, 8.0], and --param scale=KyteDoolitle records the string.
    """
    if "=" not in value:
        raise SystemExit(
            f"--param must be KEY=VALUE, got {value!r}"
        )
    key, raw = value.split("=", 1)
    key = key.strip()
    raw = raw.strip()
    if not key:
        raise SystemExit(f"--param has an empty key: {value!r}")
    if "," in raw:
        items = [item.strip() for item in raw.split(",") if item.strip()]
        return key, [_scalar(item) for item in items]
    return key, _scalar(raw)


def _scalar(raw: str) -> object:
    try:
        loaded = yaml.safe_load(raw)
    except yaml.YAMLError:
        return raw
    if isinstance(loaded, str):
        # YAML 1.1 needs a decimal point to recognise exponent notation, so
        # "1e-15" would otherwise be stored as a string and break arithmetic
        # on the page that reads it back.
        try:
            return float(loaded)
        except ValueError:
            return loaded
    return raw if loaded is None and raw != "" else loaded


def build_manifest(args, run_id: str, timestamp: datetime) -> dict:
    """Assemble the manifest in the documented key order."""
    inputs = []
    for raw in args.input:
        path = Path(raw)
        if not path.is_file():
            raise SystemExit(f"--input file not found: {raw}")
        digest, size = sha256_and_size(path)
        inputs.append(
            {"path": relative_to_root(path), "sha256": digest, "bytes": size}
        )

    externals = []
    for raw in args.external:
        path_part, description = split_annotated(raw, "--external")
        path = Path(path_part)
        if not path.is_file():
            raise SystemExit(f"--external file not found: {path_part}")
        digest, size = sha256_and_size(path)
        externals.append(
            {
                "path": relative_to_root(path),
                "sha256": digest,
                "bytes": size,
                "description": description,
            }
        )

    if args.backfilled and not args.tool:
        print(
            "warning: --backfilled without --tool; the auto-captured versions "
            "are today's, not those that produced the numbers. Record the "
            "originals with --tool NAME=VERSION where they are known.",
            file=sys.stderr,
        )

    if args.no_autotools:
        tool_versions: dict[str, str] = {}
    else:
        tool_versions = collect_tool_versions()
    for entry in args.tool:
        if "=" not in entry:
            raise SystemExit(f"--tool must be NAME=VERSION, got {entry!r}")
        name, version = entry.split("=", 1)
        tool_versions[name.strip()] = version.strip()

    parameters: dict[str, object] = {}
    for entry in args.param:
        key, value = parse_parameter(entry)
        parameters[key] = value

    if args.code_version:
        sha, dirty = args.code_version.strip(), False
    elif args.backfilled:
        # A backfilled run was produced by code we cannot identify, and
        # today's HEAD is not it. Recording null is the honest answer; pass
        # --code-version if the historical SHA is actually known.
        sha, dirty = None, False
        print(
            "note: --backfilled, so code_version is null rather than today's "
            "HEAD, which did not produce these numbers. Pass --code-version "
            "SHA if the original commit is known.",
            file=sys.stderr,
        )
    else:
        sha, dirty = git_code_version()

    manifest = {
        "run_id": run_id,
        "timestamp_utc": timestamp.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": args.track,
        "stage": args.stage,
        "description": args.description.strip(),
        "status": args.status,
        "inputs": inputs,
        "outputs": [],  # filled in below, after files are copied
        "tool_versions": tool_versions,
        "code_version": sha,
        "parameters": parameters,
        "supersedes": args.supersedes,
        "backfilled": bool(args.backfilled),
        "notes": args.notes.strip() if args.notes else "",
    }
    if sha is not None and dirty:
        # Kept separate from code_version so that field stays a bare SHA.
        manifest["code_version_dirty"] = True
    if externals:
        manifest["external_artefacts"] = externals
    return manifest


def copy_outputs(args, run_dir: Path) -> list[dict]:
    """Copy declared outputs into the run directory and describe them."""
    outputs = []
    for raw in args.copy:
        path_part, description = split_annotated(raw, "--copy")
        source = Path(path_part)
        if not source.is_file():
            raise SystemExit(f"--copy file not found: {path_part}")

        _, size = sha256_and_size(source)
        if size > MAX_COPY_BYTES and not args.allow_large:
            raise SystemExit(
                f"refusing to copy {path_part} into runs/: {size} bytes "
                f"exceeds the {MAX_COPY_BYTES} byte limit. Large artefacts "
                "belong outside the site: record them with --external, or "
                "pass --allow-large if this really is a small text result."
            )
        if source.suffix.lower() in BINARY_SUFFIXES and not args.allow_large:
            raise SystemExit(
                f"refusing to copy {path_part} into runs/: {source.suffix} "
                "looks like a binary artefact (structure, model, embedding). "
                "Record it with --external and leave the file where it lives."
            )

        destination = run_dir / source.name
        if destination.exists():
            raise SystemExit(
                f"two outputs would both be written to {destination.name}; "
                "rename one at the source"
            )
        shutil.copy2(source, destination)
        outputs.append({"path": source.name, "description": description})
    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Scaffold an immutable run directory under runs/.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--track", required=True, choices=TRACKS)
    parser.add_argument(
        "--stage",
        required=True,
        help="short label, e.g. smoke-test, characterisation, correlations",
    )
    parser.add_argument(
        "--description",
        required=True,
        help="one or two sentences on what was run and why",
    )
    parser.add_argument(
        "--slug",
        help="short kebab-case suffix for the directory name "
        "(default: derived from --stage)",
    )
    parser.add_argument(
        "--input",
        action="append",
        default=[],
        metavar="PATH",
        help="an input file consumed; hashed into the manifest (repeatable)",
    )
    parser.add_argument(
        "--copy",
        action="append",
        default=[],
        metavar='PATH[:"description"]',
        help="an output file to copy into the run directory (repeatable)",
    )
    parser.add_argument(
        "--external",
        action="append",
        default=[],
        metavar='PATH[:"description"]',
        help="a large artefact to record by path and hash without copying "
        "(repeatable)",
    )
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="an actual parameter value used, not a reference to a default "
        "(repeatable)",
    )
    parser.add_argument(
        "--tool",
        action="append",
        default=[],
        metavar="NAME=VERSION",
        help="record or override a tool version (repeatable)",
    )
    parser.add_argument("--status", choices=STATUSES, default="ok")
    parser.add_argument(
        "--supersedes",
        default=None,
        metavar="RUN_ID",
        help="run_id this run replaces",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="allow a run alongside an existing one with the same track and "
        "slug, instead of superseding it (for a genuinely separate result)",
    )
    parser.add_argument(
        "--backfilled",
        action="store_true",
        help="the numbers were transcribed from a previous report rather "
        "than produced by executing code now",
    )
    parser.add_argument("--notes", default="", help="caveats and anomalies")
    parser.add_argument(
        "--code-version",
        default=None,
        metavar="SHA",
        help="record this commit SHA instead of discovering it, for "
        "backfilling a run whose original commit is known",
    )
    parser.add_argument(
        "--timestamp",
        default=None,
        metavar="YYYY-MM-DDTHHMMZ",
        help="override the UTC timestamp, for backfilling a past run",
    )
    parser.add_argument(
        "--no-autotools",
        action="store_true",
        help="do not auto-capture local tool versions (use with "
        "--backfilled, where today's versions would be misleading)",
    )
    parser.add_argument(
        "--allow-large",
        action="store_true",
        help="permit --copy of a file the size or suffix guard rejected",
    )
    args = parser.parse_args(argv)

    if not args.description.strip():
        raise SystemExit("--description must not be empty")

    slug = args.slug or args.stage
    slug = slug.strip().lower().replace("_", "-")
    slug = re.sub(r"[^a-z0-9-]+", "-", slug).strip("-")
    if not SLUG_PATTERN.match(slug):
        raise SystemExit(
            f"could not derive a usable slug from {args.slug or args.stage!r}"
        )

    if args.timestamp:
        try:
            timestamp = datetime.strptime(
                args.timestamp, TIMESTAMP_FORMAT
            ).replace(tzinfo=timezone.utc)
        except ValueError:
            raise SystemExit(
                f"--timestamp must look like 2026-11-05T0840Z, got "
                f"{args.timestamp!r}"
            ) from None
    else:
        timestamp = datetime.now(timezone.utc)

    run_id = f"{timestamp.strftime(TIMESTAMP_FORMAT)}_{args.track}_{slug}"
    run_dir = RUNS_DIR / run_id

    if run_dir.exists():
        raise SystemExit(
            f"{run_dir} already exists. Run directories are immutable and "
            "are never edited or reused: choose a different --slug, or make "
            "a new run that --supersedes this one."
        )

    if args.supersedes and not (RUNS_DIR / args.supersedes).is_dir():
        print(
            f"warning: --supersedes names {args.supersedes}, which is not a "
            "directory under runs/",
            file=sys.stderr,
        )

    # A re-execution of the same track/stage/slug must say what it replaces.
    # Without this, re-running a stage a minute later would quietly create a
    # second, parallel run and the history would not say which is current.
    if not args.supersedes and not args.parallel:
        siblings = sorted(
            path.name
            for path in RUNS_DIR.glob(f"*_{args.track}_{slug}")
            if path.is_dir()
        )
        if siblings:
            raise SystemExit(
                f"runs/ already contains {len(siblings)} run(s) for track "
                f"{args.track!r} and slug {slug!r}:\n"
                + "\n".join(f"  {name}" for name in siblings)
                + "\n\nA re-executed or corrected analysis must name the run "
                "it replaces:\n"
                f"  --supersedes {siblings[-1]}\n"
                "If this is genuinely a separate, parallel result rather "
                "than a replacement, pass --parallel, or give it a distinct "
                "--slug."
            )

    manifest = build_manifest(args, run_id, timestamp)

    run_dir.mkdir(parents=True)
    manifest["outputs"] = copy_outputs(args, run_dir)

    manifest_path = run_dir / "manifest.yml"
    with manifest_path.open("w") as handle:
        yaml.dump(
            manifest, handle, sort_keys=False, default_flow_style=False,
            allow_unicode=True, width=78,
        )

    print(f"created {relative_to_root(run_dir)}/")
    print(f"  manifest  manifest.yml")
    print(f"  inputs    {len(manifest['inputs'])} hashed")
    print(f"  outputs   {len(manifest['outputs'])} copied in")
    if manifest.get("external_artefacts"):
        print(
            f"  external  {len(manifest['external_artefacts'])} recorded "
            "by path and hash"
        )
    if manifest["code_version"] is None:
        print("  warning   no git commit SHA available; code_version is null")
    elif manifest.get("code_version_dirty"):
        print(
            "  warning   working tree was dirty; the recorded SHA does not "
            "fully describe the code that ran"
        )
    if not manifest["parameters"]:
        print(
            "  warning   no --param values recorded; a reader cannot tell "
            "what settings produced this"
        )
    print("\nThis directory is now append-only. Validate with check_runs.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
