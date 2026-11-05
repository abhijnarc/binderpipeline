"""CLI entry point for the sequence-only physicochemical pipeline.

    python3 physchem/run_physchem.py \
        --input physchem/data/opmd_binders_clean.csv \
        --outdir physchem/out

Reads one row per binder, validates every sequence, emits two rows per
binder (as given and tag stripped), and writes the main table, a provenance
record, and any validation failures.
"""

from __future__ import annotations

import argparse
import logging
import platform
import sys
from datetime import date
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import physchem_props as pc  # noqa: E402

LOG = logging.getLogger("physchem")

# Carried through unchanged from the input, grouped after the computed
# properties. Order is the input's order.
METADATA_COLUMNS = (
    "sequence",
    "sequence_length",
    "n_term_first3",
    "has_his_tag",
    "host_strain",
    "vector",
    "temperature_C",
    "iptg_mM",
    "duration_h",
    "multiple_conditions_tested",
    "expression",
    "soluble_fraction",
    "expressed_binary",
    "soluble_binary",
    "elisa_result",
    "source_sheets",
    "data_quality_flag",
    "notes",
)

REQUIRED_COLUMNS = ("binder_id", "target", "sequence")


def read_input(path: Path) -> pd.DataFrame:
    """Read the input as strings, preserving every field exactly as written.

    Everything is read as text on purpose. Experimental metadata is carried
    through unchanged, so no numeric coercion is applied to it: sentinel
    values such as LCN2_B1's CONFLICT in host_strain and iptg_mM, and
    multi-valued fields such as "0.3; 0.5; 1.0", survive verbatim and cannot
    disturb any other row. The only field parsed as a number is
    sequence_length, and only to compare it against the real length.
    """
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        raise SystemExit(
            f"input {path} is missing required column(s): {', '.join(missing)}"
        )
    return frame


def build_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (properties table, rejected table)."""
    rows: list[dict] = []
    rejected: list[dict] = []
    seen: set[str] = set()

    for record in frame.to_dict(orient="records"):
        binder_id = record["binder_id"]
        sequence = record["sequence"]

        if binder_id in seen:
            LOG.error("%s: duplicate binder_id", binder_id)
            rejected.append(
                {
                    "binder_id": binder_id,
                    "sequence": sequence,
                    "reason": "duplicate binder_id",
                }
            )
            continue
        seen.add(binder_id)

        try:
            pc.validate_sequence(
                binder_id, sequence, record.get("sequence_length")
            )
        except pc.SequenceValidationError as exc:
            LOG.error("%s: REJECTED: %s", binder_id, exc.reason)
            rejected.append(
                {
                    "binder_id": binder_id,
                    "sequence": sequence,
                    "reason": exc.reason,
                }
            )
            continue

        LOG.info(
            "%s: validated, %d residues, all standard amino acids",
            binder_id,
            len(sequence),
        )

        for computed in pc.binder_rows(
            binder_id, record.get("target", ""), sequence
        ):
            row = dict(computed)
            for column in METADATA_COLUMNS:
                if column in record:
                    row[column] = record[column]
            rows.append(row)

    properties = pd.DataFrame(rows)
    if not properties.empty:
        computed_columns = [
            c
            for c in properties.columns
            if c not in pc.IDENTITY_COLUMNS and c not in METADATA_COLUMNS
        ]
        present_metadata = [
            c for c in METADATA_COLUMNS if c in properties.columns
        ]
        properties = properties[
            list(pc.IDENTITY_COLUMNS) + computed_columns + present_metadata
        ]
        # Deterministic row order, independent of input order.
        properties = properties.sort_values(
            ["binder_id", "construct_form"], kind="mergesort"
        ).reset_index(drop=True)

    rejects = pd.DataFrame(
        rejected, columns=["binder_id", "sequence", "reason"]
    )
    return properties, rejects


def write_provenance(path: Path, input_path: Path, n_binders: int,
                     n_rows: int, n_rejected: int) -> None:
    """Record tools, versions, date, and the exact parameters used."""
    import Bio

    lines = [
        "Sequence-only physicochemical characterisation: provenance",
        "=" * 58,
        "",
        f"Run date:            {date.today().isoformat()}",
        f"Input file:          {input_path}",
        f"Binders read:        {n_binders}",
        f"Rows written:        {n_rows} (two per accepted binder)",
        f"Rejected sequences:  {n_rejected}",
        "",
        "Tools and versions",
        "-" * 18,
        f"Python               {platform.python_version()}",
        f"Biopython            {Bio.__version__}  "
        "(Bio.SeqUtils.ProtParam, Bio.SeqUtils.ProtParamData)",
        f"peptides             {pc.peptides.__version__}",
        f"pandas               {pd.__version__}",
        "",
        "Which tool computes what",
        "-" * 24,
        "Biopython ProtParam:  molecular weight, isoelectric point, net",
        "                      charge, GRAVY, instability index,",
        "                      aromaticity, amino-acid composition,",
        "                      extinction coefficients",
        "peptides.py:          aliphatic index (ProtParam has no equivalent)",
        "This module:          Shannon entropy, longest hydrophobic run,",
        "                      sliding-window hydropathy, charged fraction,",
        "                      cysteine parity, tag stripping",
        "",
        "Parameters",
        "-" * 10,
        f"Net charge pH values         {', '.join(f'{p:g}' for p in pc.CHARGE_PH_VALUES)}",
        f"pKa table for pI and charge  {pc.PKA_SOURCE}",
        f"Hydrophobicity scale         Kyte-Doolittle "
        f"(Bio.SeqUtils.ProtParamData.kd; ProtParam scale name "
        f"'{pc.HYDROPHOBICITY_SCALE_NAME}')",
        f"Sliding window size          {pc.HYDROPHOBICITY_WINDOW} residues",
        f"Hydrophobic residue set      {pc.HYDROPHOBIC_RESIDUES} "
        "(longest-run definition)",
        f"Charged residue set          {pc.CHARGED_RESIDUES} "
        "(fraction_charged_DEKR)",
        f"Aromatic residues counted    {pc.AROMATIC_RESIDUES}",
        "Shannon entropy              base 2 (bits), over the composition of",
        "                             the 20 standard residues; order-free",
        f"Float rounding               {pc.FLOAT_PRECISION} decimal places",
        "",
        "Tag-stripping rules",
        "-" * 19,
        f"N-terminal leaders stripped  {', '.join(pc.N_TERMINAL_LEADERS)} "
        "(longest match wins)",
        "His tag                      a run of 6 or more H anchored at",
        "                             either terminus; internal runs are kept",
        f"Adjacent GS linker           a run of >= {pc.GS_LINKER_MIN_LENGTH} "
        "residues drawn only from",
        "                             G and S, removed only when directly",
        "                             adjacent to a tag that was stripped;",
        "                             internal GS linkers are structural and",
        "                             are never removed",
        f"Minimum length kept          {pc.MIN_LENGTH_AFTER_STRIPPING} "
        "residues; stripping is abandoned",
        "                             and recorded if it would leave fewer",
        "",
        "Stages registered this run",
        "-" * 26,
    ]
    for stage in pc.STAGES:
        lines.append(f"  {stage.name} (requires: {', '.join(stage.requires)})")
    lines += [
        "",
        "A structural stage (SASA, secondary structure, structure-corrected",
        "solubility) and an ESM2 embedding stage are planned and are not",
        "part of this run. No network access is used at any point.",
        "",
        "This pipeline is descriptive. Nothing here is fitted, trained, or",
        "correlated against the experimental outcome columns.",
        "",
    ]
    path.write_text("\n".join(lines))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compute sequence-only physicochemical descriptors "
        "for protein binder sequences."
    )
    default_root = Path(__file__).resolve().parent
    parser.add_argument(
        "--input",
        type=Path,
        default=default_root / "data" / "opmd_binders_clean.csv",
        help="input CSV, one row per binder (default: %(default)s)",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=default_root / "out",
        help="directory for output files (default: %(default)s)",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="console log verbosity (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    frame = read_input(args.input)
    LOG.info("read %d binders from %s", len(frame), args.input)

    properties, rejects = build_tables(frame)

    properties_path = args.outdir / "physchem_properties.tsv"
    rejects_path = args.outdir / "rejected_sequences.tsv"
    provenance_path = args.outdir / "physchem_provenance.txt"

    properties.to_csv(properties_path, sep="\t", index=False)
    rejects.to_csv(rejects_path, sep="\t", index=False)
    write_provenance(
        provenance_path, args.input, len(frame), len(properties), len(rejects)
    )

    LOG.info(
        "wrote %s (%d rows x %d columns)",
        properties_path,
        len(properties),
        len(properties.columns),
    )
    LOG.info("wrote %s (%d rejected)", rejects_path, len(rejects))
    LOG.info("wrote %s", provenance_path)

    if len(rejects):
        # Fail loudly, but only after the valid rows and the rejection
        # reasons have been written, so nothing is lost.
        print(
            f"\nERROR: {len(rejects)} sequence(s) failed validation:",
            file=sys.stderr,
        )
        for reject in rejects.to_dict(orient="records"):
            print(
                f"  {reject['binder_id']}: {reject['reason']}", file=sys.stderr
            )
        print(f"See {rejects_path}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
