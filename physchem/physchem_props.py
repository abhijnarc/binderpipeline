"""Sequence-only physicochemical descriptors for protein binder sequences.

This module is descriptive: it computes properties and returns them. It does
not fit, train, or correlate anything against experimental outcomes.

Structure
---------
Properties are computed by small single-purpose functions (one per property)
so that any one of them can be verified or swapped independently. They are
assembled by a *stage*: a named callable that takes a sequence and returns an
ordered mapping of column name to value. Stages are registered in `STAGES`.

Only the sequence-only stage exists today. A structural stage (SASA,
secondary structure, structure-corrected solubility) and an ESM2 embedding
stage are planned; each will register itself as an additional stage and
contribute its own columns, without changes to the core here.
"""

from __future__ import annotations

import logging
import math
import re
from collections import OrderedDict
from dataclasses import dataclass
from typing import Callable, Mapping

import peptides
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis

LOG = logging.getLogger("physchem")

# --- Fixed parameters. Changing any of these changes the output; they are
# --- written verbatim into the provenance file.
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
CHARGE_PH_VALUES = (7.0, 8.0)
HYDROPHOBICITY_WINDOW = 9
HYDROPHOBICITY_SCALE_NAME = "KyteDoolitle"  # Biopython's spelling
HYDROPHOBIC_RESIDUES = "AVILMFWYC"
CHARGED_RESIDUES = "DEKR"
AROMATIC_RESIDUES = "FWY"
# Biopython's charge_at_pH and isoelectric_point share one pKa table, so net
# charge and pI below are mutually consistent. peptides.py uses a different
# default scale and would not agree; it is used only for aliphatic index.
PKA_SOURCE = "Bio.SeqUtils.IsoelectricPoint (Bjellqvist/EMBOSS table)"
FLOAT_PRECISION = 4

MIN_LENGTH_AFTER_STRIPPING = 5


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------


class SequenceValidationError(ValueError):
    """A sequence failed validation; carries the binder id and the reason."""

    def __init__(self, binder_id: str, reason: str):
        self.binder_id = binder_id
        self.reason = reason
        super().__init__(f"{binder_id}: {reason}")


def validate_sequence(binder_id: str, sequence: str, declared_length) -> str:
    """Return the sequence if valid, else raise SequenceValidationError.

    Checks that the sequence is non-empty, contains only the 20 standard
    amino acids, and that its actual length matches the declared
    `sequence_length`. The offending character is named in the reason.
    """
    if sequence is None or sequence == "":
        raise SequenceValidationError(binder_id, "sequence is empty")

    if sequence != sequence.strip():
        raise SequenceValidationError(
            binder_id, "sequence has leading or trailing whitespace"
        )

    for position, char in enumerate(sequence, start=1):
        if char not in STANDARD_AA:
            raise SequenceValidationError(
                binder_id,
                f"non-standard character {char!r} at position {position} "
                f"(allowed: the 20 standard amino acids, {STANDARD_AA})",
            )

    if declared_length is not None and str(declared_length).strip() != "":
        try:
            declared = int(str(declared_length).strip())
        except ValueError:
            raise SequenceValidationError(
                binder_id,
                f"sequence_length {declared_length!r} is not an integer",
            ) from None
        if declared != len(sequence):
            raise SequenceValidationError(
                binder_id,
                f"sequence_length {declared} does not match actual "
                f"length {len(sequence)}",
            )

    return sequence


# --------------------------------------------------------------------------
# Tag and leader stripping
# --------------------------------------------------------------------------

# Longest leader first, so MGSS is preferred over MGS over M.
N_TERMINAL_LEADERS = ("MGSS", "MGS", "M")
HIS_RUN = re.compile(r"H{6,}")
# A GS-rich linker, for stripping purposes, is a run of at least three
# residues drawn only from G and S. The threshold keeps incidental GG/GP
# dipeptides at a newly exposed terminus from being treated as linker.
GS_LINKER_MIN_LENGTH = 3
GS_LINKER_RESIDUES = set("GS")


def _leading_gs_linker(sequence: str) -> str:
    run = ""
    for char in sequence:
        if char in GS_LINKER_RESIDUES:
            run += char
        else:
            break
    return run if len(run) >= GS_LINKER_MIN_LENGTH else ""


def _trailing_gs_linker(sequence: str) -> str:
    run = ""
    for char in reversed(sequence):
        if char in GS_LINKER_RESIDUES:
            run = char + run
        else:
            break
    return run if len(run) >= GS_LINKER_MIN_LENGTH else ""


@dataclass(frozen=True)
class StripResult:
    """The outcome of tag stripping: the sequence plus what was removed."""

    sequence: str
    stripped_elements: str  # ";"-joined "<what>=<substring>", or "none"
    changed: bool


def strip_tags(sequence: str) -> StripResult:
    """Remove terminal purification tags and expression leaders.

    Removes, recording each removal exactly:

    * a leading MGSS, MGS, or M expression leader;
    * a terminal run of six or more H (a His tag) at either terminus;
    * a GS-rich linker *directly adjacent to* something already removed.

    Internal GS-repeat linkers are structural (they join domains in
    single-chain multi-domain constructs) and are never touched: a linker is
    only removed when it sits at a terminus exposed by a tag that was just
    stripped.

    Stripping is abandoned if it would leave fewer than
    MIN_LENGTH_AFTER_STRIPPING residues.
    """
    original = sequence
    working = sequence
    removed: list[str] = []

    # N-terminus: His tag, then expression leader, then any linker the
    # removal exposed.
    match = HIS_RUN.match(working)
    if match:
        removed.append(f"nterm_his_tag={match.group(0)}")
        working = working[match.end() :]

    for leader in N_TERMINAL_LEADERS:
        if working.startswith(leader):
            removed.append(f"nterm_leader={leader}")
            working = working[len(leader) :]
            break

    if removed:
        linker = _leading_gs_linker(working)
        if linker:
            removed.append(f"nterm_adjacent_gs_linker={linker}")
            working = working[len(linker) :]

    # C-terminus: His tag, then any linker the removal exposed.
    c_removed_before = len(removed)
    # Only a run anchored at the very end is a tag; an internal HHHHHH is
    # part of the sequence's own composition.
    trailing_match = None
    for candidate in HIS_RUN.finditer(working):
        if candidate.end() == len(working):
            trailing_match = candidate
    if trailing_match:
        removed.append(f"cterm_his_tag={trailing_match.group(0)}")
        working = working[: trailing_match.start()]

    if len(removed) > c_removed_before:
        linker = _trailing_gs_linker(working)
        if linker:
            removed.append(f"cterm_adjacent_gs_linker={linker}")
            working = working[: len(working) - len(linker)]

    if not removed:
        return StripResult(original, "none", False)

    if len(working) < MIN_LENGTH_AFTER_STRIPPING:
        LOG.warning(
            "stripping abandoned: would leave %d residues (< %d); "
            "keeping sequence as given",
            len(working),
            MIN_LENGTH_AFTER_STRIPPING,
        )
        return StripResult(
            original,
            f"none (stripping abandoned: {';'.join(removed)} would leave "
            f"{len(working)} residues)",
            False,
        )

    return StripResult(working, ";".join(removed), True)


# --------------------------------------------------------------------------
# Individual properties. One function each, so any one can be checked or
# replaced on its own.
# --------------------------------------------------------------------------


def sequence_length(sequence: str) -> int:
    return len(sequence)


def molecular_weight(sequence: str) -> float:
    return ProteinAnalysis(sequence).molecular_weight()


def isoelectric_point(sequence: str) -> float:
    return ProteinAnalysis(sequence).isoelectric_point()


def net_charge_at_ph(sequence: str, ph: float) -> float:
    return ProteinAnalysis(sequence).charge_at_pH(ph)


def gravy(sequence: str) -> float:
    return ProteinAnalysis(sequence).gravy(scale=HYDROPHOBICITY_SCALE_NAME)


def instability_index(sequence: str) -> float:
    return ProteinAnalysis(sequence).instability_index()


def aliphatic_index(sequence: str) -> float:
    """Ikai aliphatic index. ProtParam has no equivalent; peptides.py does."""
    return peptides.Peptide(sequence).aliphatic_index()


def aromaticity(sequence: str) -> float:
    """Lobry aromaticity: the F+W+Y fraction of residues."""
    return ProteinAnalysis(sequence).aromaticity()


def residue_count(sequence: str, residue: str) -> int:
    return sequence.count(residue)


def cysteine_count(sequence: str) -> int:
    return sequence.count("C")


def cysteine_parity(sequence: str) -> str:
    """'even' or 'odd'. Odd counts leave an unpaired thiol."""
    return "even" if cysteine_count(sequence) % 2 == 0 else "odd"


def composition_percent(sequence: str) -> "OrderedDict[str, float]":
    """Percentage of each of the 20 standard amino acids."""
    counts = ProteinAnalysis(sequence).count_amino_acids()
    total = len(sequence)
    return OrderedDict(
        (aa, 100.0 * counts.get(aa, 0) / total) for aa in STANDARD_AA
    )


def fraction_charged(sequence: str) -> float:
    """Fraction of residues that are D, E, K, or R."""
    charged = sum(sequence.count(aa) for aa in CHARGED_RESIDUES)
    return charged / len(sequence)


def extinction_coefficients(sequence: str) -> tuple[int, int]:
    """(reduced, oxidised) molar extinction coefficient at 280 nm, M-1 cm-1.

    Reduced assumes all cysteines are free thiols; oxidised assumes every
    pair forms a cystine. Zero for a sequence with no W, Y, or cystine.
    """
    reduced, oxidised = ProteinAnalysis(sequence).molar_extinction_coefficient()
    return int(reduced), int(oxidised)


def shannon_entropy(sequence: str) -> float:
    """Shannon entropy of the residue composition, in bits.

    A sequence-complexity measure over composition only (it ignores residue
    order). Ranges from 0 for a homopolymer to log2(20) = 4.3219 for an even
    mix of all 20 residues.
    """
    total = len(sequence)
    entropy = 0.0
    for aa in STANDARD_AA:
        count = sequence.count(aa)
        if count:
            p = count / total
            entropy -= p * math.log2(p)
    return entropy


def longest_hydrophobic_run(sequence: str) -> int:
    """Length of the longest run of consecutive AVILMFWYC residues."""
    longest = 0
    current = 0
    for char in sequence:
        if char in HYDROPHOBIC_RESIDUES:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def max_window_hydrophobicity(
    sequence: str, window: int = HYDROPHOBICITY_WINDOW
) -> tuple[float, int, int]:
    """Maximum mean Kyte-Doolittle hydropathy over a sliding window.

    Returns (max_mean, window_start, window_centre), both positions 1-based
    and referring to the sequence as passed in. On a tie the earliest window
    wins, which keeps the output deterministic. For a sequence shorter than
    the window, the whole sequence is scored as a single window.
    """
    scale = ProtParamData.kd
    effective = min(window, len(sequence))
    values = [scale[char] for char in sequence]

    best_mean = None
    best_start = 1
    for start in range(0, len(values) - effective + 1):
        mean = sum(values[start : start + effective]) / effective
        if best_mean is None or mean > best_mean:
            best_mean = mean
            best_start = start + 1

    centre = best_start + (effective - 1) // 2
    return best_mean, best_start, centre


# --------------------------------------------------------------------------
# Stages
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Stage:
    """A named group of columns computed from a sequence.

    `compute` returns an ordered mapping of column name to value. `requires`
    names the inputs the stage needs; the sequence-only stage needs nothing
    but the sequence. Later tiers (structure, embeddings) will declare their
    own requirements and be skipped when those inputs are absent.
    """

    name: str
    compute: Callable[[str], "Mapping[str, object]"]
    requires: tuple[str, ...] = ("sequence",)


def _round(value: float) -> float:
    return round(float(value), FLOAT_PRECISION)


def compute_sequence_properties(sequence: str) -> "OrderedDict[str, object]":
    """All sequence-only descriptors for one sequence, in column order."""
    row: "OrderedDict[str, object]" = OrderedDict()

    row["length"] = sequence_length(sequence)
    row["molecular_weight"] = _round(molecular_weight(sequence))
    row["isoelectric_point"] = _round(isoelectric_point(sequence))
    for ph in CHARGE_PH_VALUES:
        label = f"net_charge_pH{ph:g}".replace(".", "_")
        row[label] = _round(net_charge_at_ph(sequence, ph))
    row["gravy"] = _round(gravy(sequence))
    row["instability_index"] = _round(instability_index(sequence))
    row["aliphatic_index"] = _round(aliphatic_index(sequence))

    row["aromaticity"] = _round(aromaticity(sequence))
    for aa in AROMATIC_RESIDUES:
        row[f"count_{aa}"] = residue_count(sequence, aa)

    row["cys_count"] = cysteine_count(sequence)
    row["cys_parity"] = cysteine_parity(sequence)

    row["fraction_charged_DEKR"] = _round(fraction_charged(sequence))

    reduced, oxidised = extinction_coefficients(sequence)
    row["ext_coeff_280_reduced"] = reduced
    row["ext_coeff_280_oxidised"] = oxidised

    row["shannon_entropy_bits"] = _round(shannon_entropy(sequence))
    row["longest_hydrophobic_run"] = longest_hydrophobic_run(sequence)

    best_mean, start, centre = max_window_hydrophobicity(sequence)
    row[f"max_mean_hydropathy_w{HYDROPHOBICITY_WINDOW}"] = _round(best_mean)
    row[f"max_mean_hydropathy_w{HYDROPHOBICITY_WINDOW}_start"] = start
    row[f"max_mean_hydropathy_w{HYDROPHOBICITY_WINDOW}_centre"] = centre

    for aa, percent in composition_percent(sequence).items():
        row[f"pct_{aa}"] = _round(percent)

    return row


SEQUENCE_STAGE = Stage(name="sequence", compute=compute_sequence_properties)

# Ordered list of active stages. Later tiers append here.
STAGES: list[Stage] = [SEQUENCE_STAGE]


def compute_all_stages(sequence: str) -> "OrderedDict[str, object]":
    """Run every registered stage over one sequence and merge the columns."""
    row: "OrderedDict[str, object]" = OrderedDict()
    for stage in STAGES:
        for key, value in stage.compute(sequence).items():
            if key in row:
                raise ValueError(
                    f"stage {stage.name!r} redefines existing column {key!r}"
                )
            row[key] = value
    return row


# --------------------------------------------------------------------------
# Per-binder assembly
# --------------------------------------------------------------------------

CONSTRUCT_AS_GIVEN = "as_given"
CONSTRUCT_TAG_STRIPPED = "tag_stripped"

IDENTITY_COLUMNS = (
    "binder_id",
    "target",
    "construct_form",
    "stripped_elements",
    "sequence_analysed",
)


def binder_rows(
    binder_id: str, target: str, sequence: str
) -> list["OrderedDict[str, object]"]:
    """Two rows for one binder: the sequence as given, and tag-stripped.

    Both rows are always emitted. When nothing was stripped the two rows are
    identical apart from `construct_form` and `stripped_elements`, which is
    itself informative: it says the construct carried no recognised tag.
    """
    rows = []
    stripped = strip_tags(sequence)

    LOG.info(
        "%s: as_given length %d; tag_stripped length %d; removed %s",
        binder_id,
        len(sequence),
        len(stripped.sequence),
        stripped.stripped_elements,
    )

    for form, seq, elements in (
        (CONSTRUCT_AS_GIVEN, sequence, "none"),
        (CONSTRUCT_TAG_STRIPPED, stripped.sequence, stripped.stripped_elements),
    ):
        row: "OrderedDict[str, object]" = OrderedDict()
        row["binder_id"] = binder_id
        row["target"] = target
        row["construct_form"] = form
        row["stripped_elements"] = elements
        row["sequence_analysed"] = seq
        row.update(compute_all_stages(seq))
        rows.append(row)

    return rows
