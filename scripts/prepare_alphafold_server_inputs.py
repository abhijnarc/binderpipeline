#!/usr/bin/env python3
"""Prepare AlphaFold Server JSON jobs from an explicit chain-sequence CSV.

Input columns: complex_id, crystal_structure, chain_id, sequence.
Rows are kept in input order; sequences are never trimmed or edited.
"""

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
REQUIRED = {"complex_id", "crystal_structure", "chain_id", "sequence"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chain_csv", type=Path, help="CSV with one row per protein chain")
    parser.add_argument("--outdir", type=Path, default=Path("alphafold3_server"))
    args = parser.parse_args()

    with args.chain_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            parser.error("CSV missing required columns: " + ", ".join(sorted(missing)))
        rows = list(reader)

    groups = defaultdict(list)
    errors = []
    seen_pairs = set()
    for lineno, raw in enumerate(rows, start=2):
        row = {key: (raw.get(key) or "").strip() for key in REQUIRED}
        cid, structure, chain, seq = (row[k] for k in ("complex_id", "crystal_structure", "chain_id", "sequence"))
        if not all((cid, structure, chain, seq)):
            errors.append(f"line {lineno}: missing complex_id, crystal_structure, chain_id, or sequence")
            continue
        pair = (cid, chain)
        if pair in seen_pairs:
            errors.append(f"line {lineno}: duplicate chain ID {chain!r} for complex {cid!r}")
            continue
        seen_pairs.add(pair)
        invalid = sorted(set(seq) - STANDARD_AA)
        if invalid:
            errors.append(f"line {lineno}: invalid amino-acid character(s) for {cid}/{chain}: {''.join(invalid)}")
            continue
        groups[cid].append(row)

    for cid, chains in groups.items():
        structures = {r["crystal_structure"] for r in chains}
        if len(chains) < 2:
            errors.append(f"{cid}: expected at least two protein chains; found {len(chains)}")
        if len(structures) != 1:
            errors.append(f"{cid}: rows refer to multiple crystal_structure values")

    if errors:
        print("Input validation failed; no jobs were written:", file=sys.stderr)
        print("\n".join(f"- {e}" for e in errors), file=sys.stderr)
        return 2
    if not groups:
        parser.error("CSV contains no valid chain rows")

    root = args.outdir
    inputs = root / "inputs"
    for name in ("inputs", "predictions", "metadata", "logs"):
        (root / name).mkdir(parents=True, exist_ok=True)
    manifest_path = inputs / "submission_manifest.csv"
    destinations = {}
    for cid in groups:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", cid)
        path = inputs / f"{safe}.json"
        if path in destinations:
            parser.error(f"complex IDs {destinations[path]!r} and {cid!r} map to the same output filename")
        destinations[path] = cid
    existing = [p for p in [manifest_path, *destinations] if p.exists()]
    if existing:
        parser.error("refusing to overwrite existing output(s): " + ", ".join(map(str, existing)))

    manifest = []
    for cid, chains in groups.items():
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", cid)
        job_path = inputs / f"{safe}.json"
        # AlphaFold Server assigns output chain IDs by entity order. It has no
        # input field for preserving PDB chain IDs, so retain an explicit map.
        job = [{
            "name": f"PDBBind_{safe}", "modelSeeds": [],
            "sequences": [{"proteinChain": {"sequence": r["sequence"], "count": 1}} for r in chains],
            "dialect": "alphafoldserver", "version": 1,
        }]
        job_path.write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")
        for index, row in enumerate(chains):
            manifest.append({
                "complex_id": cid, "crystal_structure": row["crystal_structure"],
                "input_chain_id": row["chain_id"], "server_entity_order": index + 1,
                "expected_server_chain_id": chr(ord("A") + index) if index < 26 else "see_server_output",
                "sequence": row["sequence"], "job_json": str(job_path),
                "server_job_id": "", "submission_date": "", "download_archive": "",
            })
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Prepared {len(groups)} jobs and {len(manifest)} chain mappings under {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
