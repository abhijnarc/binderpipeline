"""Copy the run directories into the rendered site, then verify the links.

Quarto will not follow a symlink out of the project directory, so `runs/`
cannot simply be declared as a resource. This runs as a post-render step: it
copies the immutable run directories into `_site/runs/` and then checks that
every `runs/<id>/` link appearing in the rendered HTML actually resolves.

Only the one copy under `runs/` is authoritative. What lands in `_site/` is
disposable and git-ignored, like the rest of the rendered output.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent
ROOT = DOCS.parent
SOURCE = ROOT / "runs"
SITE = DOCS / "_site"
DESTINATION = SITE / "runs"

# Matches the same limit new_run.py and check_runs.py enforce.
MAX_BYTES = 5 * 1024 * 1024


def main() -> int:
    if not SITE.is_dir():
        print(f"publish_runs: no rendered site at {SITE}", file=sys.stderr)
        return 1
    if not SOURCE.is_dir():
        print(f"publish_runs: no runs directory at {SOURCE}", file=sys.stderr)
        return 1

    if DESTINATION.exists():
        shutil.rmtree(DESTINATION)
    shutil.copytree(SOURCE, DESTINATION, symlinks=False)

    copied = sorted(p for p in DESTINATION.iterdir() if p.is_dir())
    total = sum(
        p.stat().st_size for p in DESTINATION.rglob("*") if p.is_file()
    )
    oversized = [
        p.relative_to(DESTINATION)
        for p in DESTINATION.rglob("*")
        if p.is_file() and p.stat().st_size > MAX_BYTES
    ]

    print(
        f"publish_runs: copied {len(copied)} run directories "
        f"({total / 1024:.0f} KiB) into _site/runs/"
    )
    if oversized:
        print(
            "publish_runs: WARNING oversized files in runs/: "
            + ", ".join(str(p) for p in oversized),
            file=sys.stderr,
        )

    # Every run link in the rendered pages must resolve to something real.
    # A broken link here means a page cites a run_id that does not exist.
    available = {path.name for path in copied}
    referenced: set[str] = set()
    pattern = re.compile(r'href="(?:\.\./)*runs/([^/"]+)/')
    for page in SITE.rglob("*.html"):
        for match in pattern.finditer(page.read_text(errors="replace")):
            referenced.add(match.group(1))

    missing = sorted(referenced - available)
    if missing:
        print(
            "publish_runs: ERROR pages link to run ids that do not exist: "
            + ", ".join(missing),
            file=sys.stderr,
        )
        return 1

    print(
        f"publish_runs: {len(referenced)} distinct run links checked, all "
        "resolve"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
