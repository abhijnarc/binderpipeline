"""Smoke test for the sequence-only physicochemical pipeline.

Checks the implementation against independently known values before it is
trusted on real data. Run directly:

    python3 physchem/test_physchem.py

Expected values are stated literally here, with their source, so that a
future failure is obvious and attributable rather than mysterious.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import physchem_props as pc  # noqa: E402

FAILURES: list[str] = []
CHECKS = 0


def check(label: str, got, expected, tolerance=None) -> None:
    """Record one comparison. Numeric when tolerance is given, else exact."""
    global CHECKS
    CHECKS += 1
    if tolerance is None:
        ok = got == expected
        detail = f"expected {expected!r}, got {got!r}"
    else:
        ok = abs(float(got) - float(expected)) <= tolerance
        detail = f"expected {expected} +/- {tolerance}, got {got}"
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {label}: {detail}")
    if not ok:
        FAILURES.append(label)


# --------------------------------------------------------------------------
# 1. Reference protein: bovine ubiquitin (UniProt P0CG53, residues 1-76).
#
# Published ExPASy ProtParam values for this 76-residue sequence:
#   molecular weight      8564.84 Da
#   theoretical pI        6.56
#   extinction (reduced)  1490 M-1 cm-1  (1 Tyr, 0 Trp, 0 Cys)
# Ubiquitin has no Trp and no Cys, so reduced and oxidised coefficients are
# equal and come from the single tyrosine: 1 x 1490.
# --------------------------------------------------------------------------
UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)


def test_reference_protein() -> None:
    print("\n1. Reference protein: bovine ubiquitin (UniProt P0CG53, 1-76)")
    check("ubiquitin length", pc.sequence_length(UBIQUITIN), 76)
    # ProtParam reports 8564.84; allow a little rounding slack.
    check(
        "ubiquitin molecular weight (ProtParam 8564.84 Da)",
        pc.molecular_weight(UBIQUITIN),
        8564.84,
        tolerance=0.5,
    )
    # ProtParam reports 6.56. Biopython's pKa table differs slightly from the
    # ExPASy web implementation, so allow 0.2 of a pH unit.
    check(
        "ubiquitin pI (ProtParam 6.56)",
        pc.isoelectric_point(UBIQUITIN),
        6.56,
        tolerance=0.2,
    )
    reduced, oxidised = pc.extinction_coefficients(UBIQUITIN)
    check("ubiquitin extinction reduced (ProtParam 1490)", reduced, 1490)
    check("ubiquitin extinction oxidised (no Cys, so 1490)", oxidised, 1490)
    check("ubiquitin Trp count (known 0)", pc.residue_count(UBIQUITIN, "W"), 0)
    check("ubiquitin Tyr count (known 1)", pc.residue_count(UBIQUITIN, "Y"), 1)
    check("ubiquitin Cys count (known 0)", pc.cysteine_count(UBIQUITIN), 0)
    check("ubiquitin Cys parity", pc.cysteine_parity(UBIQUITIN), "even")
    # Net charge near the pI must be near zero; this ties charge_at_pH and
    # isoelectric_point together and would catch a mismatched pKa table.
    at_pi = pc.net_charge_at_ph(UBIQUITIN, pc.isoelectric_point(UBIQUITIN))
    check("ubiquitin net charge at its own pI ~ 0", at_pi, 0.0, tolerance=0.05)


# --------------------------------------------------------------------------
# 2. Edge cases.
# --------------------------------------------------------------------------


def test_edge_cases() -> None:
    print("\n2. Edge cases")

    # A 5-residue peptide. MW of free Gly is 75.07; a pentaglycine loses four
    # waters: 5*75.07 - 4*18.015 = 375.35 - 72.06 = 303.29 Da.
    penta = "GGGGG"
    check("pentaglycine length", pc.sequence_length(penta), 5)
    check(
        "pentaglycine molecular weight (303.29 Da)",
        pc.molecular_weight(penta),
        303.29,
        tolerance=0.1,
    )
    check("pentaglycine entropy (one residue type, 0 bits)",
          pc.shannon_entropy(penta), 0.0, tolerance=1e-9)
    # Shorter than the window: scored as a single window of its own length.
    best, start, _ = pc.max_window_hydrophobicity(penta)
    check("pentaglycine window hydropathy (Gly = -0.4)", best, -0.4,
          tolerance=1e-9)
    check("pentaglycine window start", start, 1)

    # Zero cysteine.
    no_cys = "AEKLYWDEAKRLEEAFK"
    check("zero-Cys count", pc.cysteine_count(no_cys), 0)
    check("zero-Cys parity reported even", pc.cysteine_parity(no_cys), "even")

    # Zero aromatics: extinction coefficient must be exactly 0.
    no_aromatics = "AEKLDEAKRLEEAK"
    check("no-aromatics F count", pc.residue_count(no_aromatics, "F"), 0)
    check("no-aromatics W count", pc.residue_count(no_aromatics, "W"), 0)
    check("no-aromatics Y count", pc.residue_count(no_aromatics, "Y"), 0)
    check("no-aromatics aromaticity", pc.aromaticity(no_aromatics), 0.0,
          tolerance=1e-9)
    reduced, oxidised = pc.extinction_coefficients(no_aromatics)
    check("no-aromatics extinction reduced is 0", reduced, 0)
    check("no-aromatics extinction oxidised is 0", oxidised, 0)

    # All-alanine: GRAVY must equal the Kyte-Doolittle value for Ala (1.8),
    # entropy must be 0, and the hydrophobic run must span the sequence.
    all_ala = "AAAAAAAAAAAAAAAAAAAA"
    check("all-Ala GRAVY equals KD Ala value 1.8", pc.gravy(all_ala), 1.8,
          tolerance=1e-9)
    check("all-Ala entropy 0 bits", pc.shannon_entropy(all_ala), 0.0,
          tolerance=1e-9)
    check("all-Ala longest hydrophobic run equals length",
          pc.longest_hydrophobic_run(all_ala), 20)
    check("all-Ala composition pct_A is 100",
          pc.composition_percent(all_ala)["A"], 100.0, tolerance=1e-9)
    check("all-Ala fraction charged is 0", pc.fraction_charged(all_ala), 0.0,
          tolerance=1e-9)

    # Composition percentages must sum to 100 for any valid sequence.
    total = sum(pc.composition_percent(UBIQUITIN).values())
    check("ubiquitin composition sums to 100%", total, 100.0, tolerance=1e-6)

    # A sequence containing X must be rejected, naming the character.
    try:
        pc.validate_sequence("TEST_X", "AEKLXDEAK", 9)
        check("sequence with X rejected", "not rejected", "rejected")
    except pc.SequenceValidationError as exc:
        check("sequence with X rejected", True, True)
        check("X rejection names the character", "'X'" in exc.reason, True)
        print(f"         reason: {exc.reason}")

    # A sequence containing a hyphen must be rejected.
    try:
        pc.validate_sequence("TEST_HYPHEN", "AEKL-DEAK", 9)
        check("sequence with hyphen rejected", "not rejected", "rejected")
    except pc.SequenceValidationError as exc:
        check("sequence with hyphen rejected", True, True)
        check("hyphen rejection names the character", "'-'" in exc.reason, True)
        print(f"         reason: {exc.reason}")

    # A declared length that disagrees with the real length must be rejected.
    try:
        pc.validate_sequence("TEST_LEN", "AEKLDEAK", 99)
        check("length mismatch rejected", "not rejected", "rejected")
    except pc.SequenceValidationError as exc:
        check("length mismatch rejected", True, True)
        print(f"         reason: {exc.reason}")

    # A valid sequence must pass untouched.
    check("valid sequence passes validation",
          pc.validate_sequence("TEST_OK", "AEKLDEAK", 8), "AEKLDEAK")


# --------------------------------------------------------------------------
# 3. Tag stripping on constructed examples.
# --------------------------------------------------------------------------


def test_tag_stripping() -> None:
    print("\n3. Tag stripping")
    core = "AEKLYWDEAKRLEEAFKRLEEAK"

    # N-terminal His tag.
    result = pc.strip_tags("HHHHHH" + core)
    check("N-term His tag removed", result.sequence, core)
    check("N-term His tag recorded", result.stripped_elements,
          "nterm_his_tag=HHHHHH")

    # C-terminal His tag.
    result = pc.strip_tags(core + "HHHHHH")
    check("C-term His tag removed", result.sequence, core)
    check("C-term His tag recorded", result.stripped_elements,
          "cterm_his_tag=HHHHHH")

    # MGS leader.
    result = pc.strip_tags("MGS" + core)
    check("MGS leader removed", result.sequence, core)
    check("MGS leader recorded", result.stripped_elements, "nterm_leader=MGS")

    # MGSS leader: the longer leader must win over MGS.
    result = pc.strip_tags("MGSS" + core)
    check("MGSS leader removed", result.sequence, core)
    check("MGSS leader recorded", result.stripped_elements,
          "nterm_leader=MGSS")

    # Bare initiator Met.
    result = pc.strip_tags("M" + core)
    check("initiator M removed", result.sequence, core)
    check("initiator M recorded", result.stripped_elements, "nterm_leader=M")

    # Both termini tagged at once.
    result = pc.strip_tags("MGS" + core + "HHHHHH")
    check("both termini stripped", result.sequence, core)
    check("both removals recorded", result.stripped_elements,
          "nterm_leader=MGS;cterm_his_tag=HHHHHH")

    # A sequence whose own composition contains HH but is not a tag. Also
    # contains an internal HHHHHH-like stretch nowhere near a terminus.
    own_hh = "AEKLHHDEAKRLHHEEAFK"
    result = pc.strip_tags(own_hh)
    check("internal HH not treated as tag", result.sequence, own_hh)
    check("internal HH records no stripping", result.stripped_elements, "none")
    check("internal HH marked unchanged", result.changed, False)

    # Six H in the middle of a sequence is still not a terminal tag.
    internal_six = "AEKLDE" + "HHHHHH" + "AKRLEEAFK"
    result = pc.strip_tags(internal_six)
    check("internal HHHHHH not stripped", result.sequence, internal_six)
    check("internal HHHHHH records no stripping", result.stripped_elements,
          "none")

    # A His tag that ends the sequence but is longer than six.
    result = pc.strip_tags(core + "HHHHHHHH")
    check("8xHis tag removed whole", result.sequence, core)
    check("8xHis tag recorded in full", result.stripped_elements,
          "cterm_his_tag=HHHHHHHH")

    # Internal GGGGS linkers must survive stripping. This mirrors the real
    # VH-linker-VH construct: two domains joined by GS repeats, no tag.
    vh = "QVQLVESGGGLVQPGGSLRLSCAAS"
    linker = "GGGGSGGGGSGGGGSGGGGS"
    two_domain = vh + linker + vh
    result = pc.strip_tags(two_domain)
    check("internal GGGGS linker survives (no tag present)",
          result.sequence, two_domain)
    check("internal linker records no stripping", result.stripped_elements,
          "none")

    # Same construct with a C-terminal His tag: the tag goes, the internal
    # linker stays.
    result = pc.strip_tags(two_domain + "HHHHHH")
    check("tag removed from linkered construct", result.sequence, two_domain)
    check("internal linker still present after stripping",
          linker in result.sequence, True)
    check("only the tag recorded", result.stripped_elements,
          "cterm_his_tag=HHHHHH")

    # A GS linker *adjacent to* a removed tag is stripped, unlike an
    # internal one.
    result = pc.strip_tags("MGS" + "GGGGS" + core)
    check("GS linker adjacent to removed leader is stripped",
          result.sequence, core)
    check("adjacent linker recorded separately", result.stripped_elements,
          "nterm_leader=MGS;nterm_adjacent_gs_linker=GGGGS")

    # A short GG at a newly exposed terminus is not a linker: MGP-style
    # starts must keep their G.
    result = pc.strip_tags("MGGPFDYFWEIMWVVAQ")
    check("short GG after M not treated as linker",
          result.sequence, "GGPFDYFWEIMWVVAQ")
    check("short GG records only the leader", result.stripped_elements,
          "nterm_leader=M")

    # Stripping that would leave almost nothing is abandoned, not applied.
    result = pc.strip_tags("MGSHHHHHH")
    check("over-stripping abandoned", result.sequence, "MGSHHHHHH")
    check("abandonment marked unchanged", result.changed, False)
    check("abandonment explained",
          "abandoned" in result.stripped_elements, True)


def main() -> int:
    print("Sequence-only physicochemical pipeline: smoke test")
    print(f"peptides {pc.peptides.__version__}; Biopython in use for ProtParam")

    test_reference_protein()
    test_edge_cases()
    test_tag_stripping()

    print(f"\n{CHECKS} checks run; {len(FAILURES)} failed")
    if FAILURES:
        print("FAILED checks:")
        for name in FAILURES:
            print(f"  - {name}")
        print("\nSMOKE TEST FAILED - do not run on real data.")
        return 1
    print("SMOKE TEST PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
