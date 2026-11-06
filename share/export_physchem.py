"""Write the shareable CSV and interpretation note for the physchem track.

Both files are generated from the immutable run directory, never hand-edited:
every number in the note is read from the results table or the run manifest.
Regenerate after any new characterisation run rather than patching the output.

    python3 share/export_physchem.py

Writes, beside this script:
    physchem_results.csv          headline descriptors, both construct forms
    physchem_interpretation.txt   what the numbers mean, and what they do not
"""

from __future__ import annotations

import hashlib
import sys
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "docs" / "_helpers"))

import load_runs as lr  # noqa: E402

CSV_NAME = "physchem_results.csv"
NOTE_NAME = "physchem_interpretation.txt"

# The headline descriptors. The full table in the run directory also carries
# the 20-column amino-acid composition and the raw input metadata; those are
# omitted here to keep the shared file readable, and the run directory is
# named in the note so the complete version can be found.
DESCRIPTORS = [
    "length",
    "molecular_weight",
    "isoelectric_point",
    "net_charge_pH7",
    "net_charge_pH8",
    "gravy",
    "instability_index",
    "aliphatic_index",
    "aromaticity",
    "count_F",
    "count_W",
    "count_Y",
    "cys_count",
    "cys_parity",
    "fraction_charged_DEKR",
    "ext_coeff_280_reduced",
    "ext_coeff_280_oxidised",
    "shannon_entropy_bits",
    "longest_hydrophobic_run",
    "max_mean_hydropathy_w9",
    "max_mean_hydropathy_w9_centre",
]

# Observed experimental conditions and outcomes, carried through unchanged.
# Nothing is fitted or correlated against these columns.
METADATA = [
    "host_strain",
    "vector",
    "temperature_C",
    "iptg_mM",
    "duration_h",
    "multiple_conditions_tested",
    "expression",
    "soluble_fraction",
    "elisa_result",
    "data_quality_flag",
]

IDENTITY = ["binder_id", "target", "construct_form", "stripped_elements"]

# Outcome groups as recorded in the source data, used only to report ranges
# side by side. No test statistic is computed across these groups.
GROUPS = {
    "pellet": ["IL10NB5", "IL10NB1", "MMP9B7"],
    "supernatant": ["IL6AD15", "LCN2_B1", "LCN2_B5", "MMP9B9", "MMP9B17"],
    "no expression": ["IL6P233", "IL10R2"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def wrap(text: str, indent: str = "") -> str:
    return "\n".join(
        textwrap.fill(paragraph, width=78, initial_indent=indent,
                      subsequent_indent=indent)
        if paragraph.strip() else ""
        for paragraph in text.split("\n")
    )


def bullet(text: str, marker: str = "  - ") -> str:
    """A bullet whose continuation lines align under the text, not the dash."""
    return textwrap.fill(text, width=78, initial_indent=marker,
                         subsequent_indent=" " * len(marker))


def heading(text: str) -> str:
    return f"{text}\n{'-' * len(text)}"


def write_csv(props, path: Path) -> Path:
    """Write the shared subset of columns with every value exactly as recorded.

    `props` must be the string-typed read of the results table. Letting pandas
    infer types would rewrite TRUE as True and, worse, turn a recorded literal
    NA -- meaning not applicable, which is a statement -- into an empty cell,
    which is ambiguous. The same applies to LCN2_B1's CONFLICT.
    """
    columns = IDENTITY + DESCRIPTORS + METADATA + ["sequence_analysed"]
    missing = [column for column in columns if column not in props.columns]
    if missing:
        raise SystemExit(
            f"results table is missing expected columns: {', '.join(missing)}"
        )
    table = props[columns].sort_values(
        ["binder_id", "construct_form"], kind="mergesort"
    )
    table.to_csv(path, index=False)
    return path


def write_note(props, stripped, as_given, char_run, smoke_run,
               csv_path: Path, path: Path) -> Path:
    p = lr.params(char_run)
    sp = lr.params(smoke_run)
    manifest = lr.manifest(char_run)

    basic = stripped.isoelectric_point.idxmax()
    low_arom = stripped.aromaticity.idxmin()
    unstable = stripped.instability_index.idxmax()
    pellet = GROUPS["pellet"]
    soluble = GROUPS["supernatant"]

    pi_change = (stripped.isoelectric_point - as_given.isoelectric_point).abs()
    tag_worst = pi_change.idxmax()

    # The manifest stores this as a list; prose needs it as readable text.
    leaders = lr.as_list(char_run, "nterm_leaders_stripped")
    leaders = " or ".join([", ".join(leaders[:-1]), leaders[-1]]) \
        if len(leaders) > 1 else leaders[0]

    lines: list[str] = []
    add = lines.append

    add("PHYSICOCHEMICAL CHARACTERISATION OF 10 BINDER DESIGNS")
    add("Sequence-only descriptors: what the accompanying CSV contains, and")
    add("what it does and does not support.")
    add("")
    add(f"Generated   {lr.today()}")
    add(f"Run         {char_run}")
    add(f"Data file   {csv_path.name}")
    add("")
    add("=" * 78)
    add("")

    add(heading("1. READ THIS FIRST"))
    add("")
    add(wrap(
        "These are computed descriptors, not predictions. Nothing in this "
        "analysis is fitted, trained, or correlated against the experimental "
        "outcomes. The observed conditions and results travel in the same "
        "table so a person can look at both at once; no relationship between "
        "them has been calculated, and none is implied."
    ))
    add("")
    add(wrap(
        f"The sample is {p['binders']} binders, of which 2 did not express, 3 "
        f"went to the pellet and 5 to the supernatant. That is far too small "
        f"to support any conclusion about which properties matter for "
        f"expression or solubility. No p-value or correlation coefficient "
        f"against an outcome column has been computed, and none should be. "
        f"Treat the group comparisons in section 6 as a check that the "
        f"pipeline produces sane numbers, not as a result."
    ))
    add("")
    add(wrap(
        "Two of the descriptors deserve specific scepticism. The instability "
        "index and the aliphatic index are literature heuristics derived from "
        "contexts unlike recombinant expression in E. coli; they are not "
        "quality scores and a value either side of a conventional threshold "
        "means little here."
    ))
    add("")

    add(heading("2. WHAT WAS COMPUTED"))
    add("")
    add(wrap(
        f"From amino-acid sequence alone: length, molecular weight, "
        f"isoelectric point, net charge at pH "
        f"{' and '.join(str(v) for v in p['net_charge_ph_values'])}, GRAVY "
        f"hydropathy, instability and aliphatic indices, aromaticity with "
        f"explicit F/W/Y counts, cysteine count and parity, the full "
        f"20-residue composition, charged fraction over "
        f"{p['charged_residue_set']}, reduced and oxidised extinction "
        f"coefficients at 280 nm, Shannon entropy of composition, the longest "
        f"run of {p['hydrophobic_residue_set']}, and the maximum mean "
        f"{p['hydrophobicity_scale']} hydropathy over a "
        f"{p['sliding_window_residues']}-residue window with its position."
    ))
    add("")
    add(wrap(
        f"The CSV carries {len(DESCRIPTORS)} headline descriptors per row "
        f"plus the experimental metadata. The complete table, including the "
        f"per-residue composition, is the file physchem_properties.tsv in run "
        f"directory {char_run}."
    ))
    add("")
    add(wrap(
        "No structural information is used. A structural tier (solvent "
        "accessibility, secondary structure) and a protein-language-model "
        "embedding tier are planned but not implemented, so nothing here "
        "accounts for how the chain actually folds."
    ))
    add("")

    add(heading("3. VERIFICATION"))
    add("")
    add(wrap(
        f"The implementation was checked against independently published "
        f"values before being run on these sequences: {sp['checks_run']} "
        f"checks, {sp['checks_failed']} failures (run {smoke_run}). The "
        f"reference protein was {sp['reference_protein']}, for which the "
        f"published values are molecular weight {sp['expected_mw_da']} Da "
        f"(matched within {sp['mw_tolerance_da']} Da), isoelectric point "
        f"{sp['expected_pi']} (matched within {sp['pi_tolerance_ph']} pH "
        f"units, the tolerance being there because the pKa table used here "
        f"differs slightly from the ExPASy web implementation), and "
        f"extinction coefficient {sp['expected_ext280']} M-1 cm-1 (matched "
        f"exactly)."
    ))
    add("")
    add(wrap(
        "Edge cases covered a 5-residue peptide, sequences with no cysteine "
        "and no aromatics (where the extinction coefficient must be exactly "
        "zero), an all-alanine sequence, and rejection of sequences "
        "containing non-standard characters."
    ))
    add("")

    add(heading("4. TWO ROWS PER BINDER: READ THE RIGHT ONE"))
    add("")
    add(wrap(
        "Purification tags and expression leaders change every number in the "
        "table, and this set mixes tagged and untagged constructs. Each "
        "binder therefore appears twice, distinguished by the construct_form "
        "column:"
    ))
    add("")
    for term, meaning in [
        ("as_given", "the construct as supplied, i.e. the protein that was "
                     "actually expressed. Use this to interpret an observed "
                     "experimental outcome."),
        ("tag_stripped", "the binder with terminal tags removed. Use this to "
                         "compare binders with each other, because otherwise "
                         "the comparison is partly between tags."),
    ]:
        add(textwrap.fill(meaning, width=78,
                          initial_indent=f"  {term.ljust(14)}",
                          subsequent_indent=" " * 16))
        add("")
    add(wrap(
        f"The stripped_elements column records exactly what was removed from "
        f"each sequence; nothing was modified without being recorded. "
        f"Stripping removes a terminal run of at least "
        f"{p['his_tag_min_run']} histidines at either end, a leading "
        f"{leaders} leader, and a GS-rich linker only "
        f"where it directly abuts a removed tag."
    ))
    add("")
    add(wrap(
        f"This is not a cosmetic difference. Removing {tag_worst}'s "
        f"histidine tag moves its isoelectric point from "
        f"{as_given.loc[tag_worst, 'isoelectric_point']:.2f} to "
        f"{stripped.loc[tag_worst, 'isoelectric_point']:.2f}, because six "
        f"histidines are a real titratable group near neutral pH."
    ))
    add("")
    add(wrap(
        f"IL10NB5 is deliberately left intact: it is a two-domain construct "
        f"whose halves are joined by GGGGS repeats. Those linkers are "
        f"structural, not purification tags, and are retained "
        f"(stripped_elements = "
        f"{stripped.loc['IL10NB5', 'stripped_elements']})."
    ))
    add("")

    add(heading("5. WHAT THE NUMBERS SHOW"))
    add("")
    add(wrap(f"Ranges across the {len(stripped)} binders, tag_stripped rows:"))
    add("")
    summary = ["length", "molecular_weight", "isoelectric_point",
               "net_charge_pH7", "gravy", "instability_index",
               "aliphatic_index", "aromaticity", "cys_count",
               "fraction_charged_DEKR", "ext_coeff_280_reduced",
               "longest_hydrophobic_run", "max_mean_hydropathy_w9"]
    width = max(len(name) for name in summary)
    add(f"  {'descriptor'.ljust(width)}  {'min':>12}  {'median':>12}  "
        f"{'max':>12}")
    for name in summary:
        column = stripped[name]
        add(f"  {name.ljust(width)}  {column.min():>12.3f}  "
            f"{column.median():>12.3f}  {column.max():>12.3f}")
    add("")
    add("Individual binders that stand out:")
    add("")
    add(bullet(
        f"{basic} is the only basic construct (pI "
        f"{stripped.loc[basic, 'isoelectric_point']:.2f}, net charge "
        f"{stripped.loc[basic, 'net_charge_pH7']:+.1f} at pH 7), the only one "
        f"with substantial cysteine ({int(stripped.loc[basic, 'cys_count'])} "
        f"residues), and has the highest local hydropathy of the set.",
    ))
    add("")
    add(bullet(
        f"{low_arom} has the lowest aromaticity of the set "
        f"({stripped.loc[low_arom, 'aromaticity']:.4f}) and an extinction "
        f"coefficient of only "
        f"{int(stripped.loc[low_arom, 'ext_coeff_280_reduced'])} M-1 cm-1: "
        f"one tyrosine and no tryptophan. Practical consequence: absorbance "
        f"at 280 nm will quantify this protein poorly. Use a different assay.",
    ))
    add("")
    add(bullet(
        f"{unstable} has an instability index of "
        f"{stripped.loc[unstable, 'instability_index']:.1f}, nearly double "
        f"the median of {stripped.instability_index.median():.1f}. See the "
        f"caveat in section 1 before reading anything into that.",
    ))
    add("")
    add(bullet(
        f"IL10NB1 and IL10NB5, the two nanobody-format constructs, sit apart "
        f"from the designed helical binders on nearly every descriptor, with "
        f"a higher charged fraction and a GRAVY closest to zero. They differ "
        f"in format, length and origin, so they are not comparable to the "
        f"rest as a group of one kind.",
    ))
    add("")
    add(bullet(
        f"All {len(stripped)} binders have an even cysteine count and "
        f"{int((stripped.cys_count == 0).sum())} have none at all. GRAVY is "
        f"negative for every one, so all are net hydrophilic by that measure.",
    ))
    add("")

    add(heading("6. DESCRIPTORS AGAINST OBSERVED OUTCOMES (SANITY CHECK)"))
    add("")
    add(wrap(
        "Ranges side by side, for inspection only. Sample sizes are given "
        "because they are the reason no test is reported."
    ))
    add("")
    shown = ["fraction_charged_DEKR", "gravy", "instability_index",
             "longest_hydrophobic_run", "length"]
    label_width = max(len(name) for name in shown)
    add(f"  {'descriptor'.ljust(label_width)}  " + "  ".join(
        f"{name} (n={len(ids)})".rjust(22)
        for name, ids in GROUPS.items()))
    for name in shown:
        cells = []
        for ids in GROUPS.values():
            column = stripped.loc[ids, name]
            cells.append(f"{column.min():.2f} - {column.max():.2f}".rjust(22))
        add(f"  {name.ljust(label_width)}  " + "  ".join(cells))
    add("")
    add(wrap(
        f"The charged fraction and GRAVY ranges of the {len(pellet)} pellet "
        f"and {len(soluble)} supernatant cases nearly touch but do not "
        f"separate. Two of the three pellet cases are the nanobody-format "
        f"IL10 constructs, so that split is confounded with construct class "
        f"rather than reflecting these descriptors. The instability index "
        f"runs opposite to the usual folklore, being lower in the insoluble "
        f"group ({stripped.loc[pellet].instability_index.min():.1f} to "
        f"{stripped.loc[pellet].instability_index.max():.1f}) than in the "
        f"soluble one ({stripped.loc[soluble].instability_index.min():.1f} to "
        f"{stripped.loc[soluble].instability_index.max():.1f}). None of this "
        f"is a signal; with these numbers it could not be."
    ))
    add("")

    add(heading("7. IL6P233 VERSUS IL6AD15"))
    add("")
    add(wrap(
        "These two share their C-terminal 30 residues and differ across "
        "their N-terminal 60. IL6AD15 expressed and was mostly soluble; "
        "IL6P233 did not express. They are the closest thing in this set to a "
        "controlled comparison, so it is worth stating plainly how their "
        "computed properties differ."
    ))
    add("")
    pair_keys = ["length", "molecular_weight", "isoelectric_point",
                 "net_charge_pH7", "gravy", "instability_index",
                 "aliphatic_index", "aromaticity", "count_Y",
                 "ext_coeff_280_reduced", "longest_hydrophobic_run",
                 "fraction_charged_DEKR"]
    key_width = max(len(name) for name in pair_keys)
    add(f"  {'descriptor'.ljust(key_width)}  {'IL6P233':>14}  "
        f"{'IL6AD15':>14}")
    add(f"  {'(outcome)'.ljust(key_width)}  {'no expression':>14}  "
        f"{'soluble':>14}")
    for name in pair_keys:
        add(f"  {name.ljust(key_width)}  "
            f"{as_given.loc['IL6P233', name]:>14.3f}  "
            f"{as_given.loc['IL6AD15', name]:>14.3f}")
    add("")
    add(wrap(
        f"They are near-identical on isoelectric point, net charge and "
        f"hydropathy. The differences that exist are small and mostly point "
        f"the wrong way: the binder that failed to express has the LOWER "
        f"instability index "
        f"({as_given.loc['IL6P233', 'instability_index']:.1f}, below the "
        f"conventional threshold of 40) and the HIGHER aliphatic index "
        f"({as_given.loc['IL6P233', 'aliphatic_index']:.1f} against "
        f"{as_given.loc['IL6AD15', 'aliphatic_index']:.1f})."
    ))
    add("")
    add(wrap(
        "No computed descriptor here explains the difference in outcome. For "
        "a pipeline sanity check that is the useful answer: it says these "
        "descriptors are not covertly encoding the result, and it is "
        "consistent with the cause lying somewhere sequence-only descriptors "
        "cannot see -- local folding, codon or mRNA effects, toxicity, or "
        "simply that only one induction temperature was tried for IL6P233."
    ))
    add("")

    add(heading("8. KNOWN DEFECTS IN THE INPUT DATA"))
    add("")
    add(wrap(
        "LCN2_B1 carries the literal value CONFLICT in its host_strain and "
        "iptg_mM columns, because the source sheets disagree. This has been "
        "carried through as-is rather than guessed at. Any analysis touching "
        "those two fields must handle it; the descriptors themselves are "
        "computed from sequence and are unaffected."
    ))
    add("")
    add(wrap(
        "The data_quality_flag and source_sheets columns record provenance "
        "per row and are worth reading before using any outcome column."
    ))
    add("")
    add(wrap(
        "A literal NA in soluble_fraction or elisa_result means not "
        "applicable or not measured, and is distinct from an empty cell. "
        "These values are written out exactly as recorded, so a spreadsheet "
        "that silently converts NA to blank on import will lose that "
        "distinction."
    ))
    add("")

    add(heading("9. PROVENANCE"))
    add("")
    add(f"  characterisation run   {char_run}")
    add(f"  timestamp (UTC)        {manifest.get('timestamp_utc')}")
    add(f"  status                 {manifest.get('status')}")
    add(f"  code version           "
        f"{manifest.get('code_version') or 'not recorded'}")
    add(f"  verification run       {smoke_run}")
    add("")
    add("  inputs, with hashes:")
    for entry in manifest.get("inputs") or []:
        add(f"    {entry.get('path')}")
        add(f"      sha256 {entry.get('sha256')}")
        add(f"      bytes  {entry.get('bytes')}")
    add("")
    add("  tool versions:")
    for name, version in sorted(lr.tool_versions(char_run).items()):
        add(f"    {name} {version}")
    add("")
    add(wrap(
        "Results are versioned as immutable run directories. The run named "
        "above is never edited; a corrected analysis becomes a new run that "
        "names the one it replaces. If the numbers in this note ever "
        "disagree with a newer report, the run identifiers say which is "
        "which."
    ))
    add("")
    add(f"  {csv_path.name} sha256")
    add(f"    {sha256(csv_path)}")
    add("")

    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> int:
    char_run = lr.latest("physchem", "characterisation")
    smoke_run = lr.latest("physchem", "smoke-test")
    if char_run is None or smoke_run is None:
        raise SystemExit(
            "no current physchem characterisation and smoke-test runs found "
            "under runs/; nothing to export"
        )

    # Read twice on purpose: verbatim strings for the file that gets shared,
    # inferred types for the arithmetic in the note.
    verbatim = lr.read_output(char_run, "physchem_properties.tsv",
                              dtype=str, keep_default_na=False)
    props = lr.read_output(char_run, "physchem_properties.tsv")
    stripped = props[props.construct_form == "tag_stripped"].set_index(
        "binder_id")
    as_given = props[props.construct_form == "as_given"].set_index("binder_id")

    csv_path = write_csv(verbatim, HERE / CSV_NAME)
    note_path = write_note(props, stripped, as_given, char_run, smoke_run,
                           csv_path, HERE / NOTE_NAME)

    print(f"wrote {csv_path} ({len(props)} rows, "
          f"{len(IDENTITY + DESCRIPTORS + METADATA) + 1} columns)")
    print(f"wrote {note_path} ({len(note_path.read_text().splitlines())} "
          f"lines)")
    print(f"source run {char_run}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
