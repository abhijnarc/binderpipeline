# Sequence-only physicochemical characterisation of binder sequences

## What this is, and what it is not

This directory holds a **descriptive** table of physicochemical properties
computed from amino-acid sequences alone.

**These columns are descriptors. They are not a prediction of expression,
solubility, yield, or binding.** Nothing in this pipeline is fitted, trained,
or correlated against the experimental outcome columns. The experimental
metadata (`expression`, `soluble_fraction`, `expressed_binary`,
`soluble_binary`, `elisa_result`) is carried along so that a human can look at
descriptors and outcomes side by side; no relationship between the two is
computed, asserted, or implied. Several of these descriptors — the instability
index and the aliphatic index in particular — are heuristics from the
literature with weak and context-dependent relevance to recombinant
expression in *E. coli*, and should not be read as quality scores.

## How to run

```bash
pip install --break-system-packages peptides biopython pandas

python3 physchem/test_physchem.py          # smoke test; must pass first
python3 physchem/run_physchem.py           # writes this directory
```

Everything runs locally; no network access is used at any point. Output is
deterministic — the same input gives a byte-identical table.

Files: [`physchem_props.py`](../physchem_props.py) (module),
[`run_physchem.py`](../run_physchem.py) (CLI),
[`test_physchem.py`](../test_physchem.py) (smoke test).

## Outputs

| File | Contents |
|---|---|
| `physchem_properties.tsv` | Main table, one row per (binder, construct form) |
| `physchem_provenance.txt` | Tool versions, run date, exact parameters used |
| `rejected_sequences.tsv` | Validation failures with reasons (empty if all passed) |
| `README.md` | This file |

## Column groups

The main table has three contiguous blocks of columns, in this order.

### 1. Identity

| Column | Meaning |
|---|---|
| `binder_id` | Binder identifier, from the input |
| `target` | Intended target, from the input |
| `construct_form` | `as_given` or `tag_stripped` (see below) |
| `stripped_elements` | Exactly what was removed, or `none` |
| `sequence_analysed` | The sequence every computed column below refers to |

### 2. Computed properties

All computed from `sequence_analysed`, rounded to 4 decimal places.

| Column | Meaning |
|---|---|
| `length` | Number of residues |
| `molecular_weight` | Average molecular weight, Da |
| `isoelectric_point` | Theoretical pI |
| `net_charge_pH7` | Net charge at pH 7.0 |
| `net_charge_pH8` | Net charge at pH 8.0 (typical lysis buffer) |
| `gravy` | Grand average of hydropathicity: mean Kyte-Doolittle value per residue. Positive is more hydrophobic |
| `instability_index` | Guruprasad dipeptide-composition heuristic. A literature convention calls > 40 "unstable"; it is a sequence heuristic, not a measurement of folding stability or half-life |
| `aliphatic_index` | Ikai aliphatic index: relative volume occupied by A, V, I, L. Historically associated with thermostability in globular proteins |
| `aromaticity` | Lobry aromaticity: the F + W + Y fraction of residues |
| `count_F`, `count_W`, `count_Y` | Explicit counts of each aromatic residue |
| `cys_count` | Number of cysteines |
| `cys_parity` | `even` or `odd`. An odd count leaves at least one unpaired thiol; the *E. coli* cytoplasm is reducing, so this is about disulfide potential, not about what forms in vivo |
| `fraction_charged_DEKR` | Fraction of residues that are D, E, K, or R |
| `ext_coeff_280_reduced` | Molar extinction coefficient at 280 nm, M⁻¹cm⁻¹, with all cysteines reduced |
| `ext_coeff_280_oxidised` | Same, assuming every cysteine pair forms a cystine. Equal to the reduced value when there is no cysteine; `0` for a sequence with no W, Y, or cystine |
| `shannon_entropy_bits` | Shannon entropy of the residue composition, in bits. 0 for a homopolymer, log₂20 = 4.3219 for an even mix of all 20 residues. A composition-only complexity measure: it ignores residue order, so it flags compositional bias but not periodicity or repeats |
| `longest_hydrophobic_run` | Length of the longest run of consecutive residues from AVILMFWYC |
| `max_mean_hydropathy_w9` | Highest mean Kyte-Doolittle hydropathy over any 9-residue window |
| `max_mean_hydropathy_w9_start` | 1-based position of that window's first residue |
| `max_mean_hydropathy_w9_centre` | 1-based position of that window's central residue |
| `pct_A` … `pct_Y` | Full amino-acid composition as percentages, 20 columns, alphabetical by one-letter code. Sums to 100 |

Positions in the hydropathy columns refer to `sequence_analysed`, so they
differ between the two construct forms of the same binder whenever an
N-terminal element was stripped.

### 3. Experimental metadata

`sequence`, `sequence_length`, `n_term_first3`, `has_his_tag`, `host_strain`,
`vector`, `temperature_C`, `iptg_mM`, `duration_h`,
`multiple_conditions_tested`, `expression`, `soluble_fraction`,
`expressed_binary`, `soluble_binary`, `elisa_result`, `source_sheets`,
`data_quality_flag`, `notes`.

These are **carried through from the input unchanged** and are read as text,
never coerced to numbers. This matters for two reasons:

- `LCN2_B1` carries `CONFLICT` in `host_strain` and `iptg_mM` because the
  source sheets disagree (Decision Tree: BL21, 0.5 mM; Marker tracker:
  Rosetta, 0.3 mM). That sentinel is preserved verbatim. No value is imputed,
  and it cannot disturb any other row, because no numeric parsing is applied
  to these columns at all.
- Multi-valued fields such as `iptg_mM` = `0.3; 0.5; 1.0` and `temperature_C`
  = `37; 30` record several tested conditions and likewise survive as written.

Note that `sequence` (metadata, the input sequence) and `sequence_analysed`
(identity, what was actually computed on) differ in `tag_stripped` rows.

## Construct forms and tag handling

Purification tags and expression leaders change every number in the table, and
this dataset mixes tagged and untagged constructs. So each binder appears
**twice**:

- `as_given` — the sequence exactly as supplied.
- `tag_stripped` — with terminal tags and leaders removed.

Both rows are always emitted. When nothing was stripped the two rows are
identical except for `construct_form` and `stripped_elements`; that is itself
informative, as it says the construct carried no recognised terminal tag.
Compare like with like: use `tag_stripped` when comparing binders to each
other, and `as_given` when reasoning about the protein that was actually
expressed and purified.

Stripping rules, in order:

1. A terminal run of **six or more H** at either terminus (a His tag).
2. A leading **MGSS**, **MGS**, or **M** (longest match wins).
3. A **GS-rich linker** — a run of ≥ 3 residues drawn only from G and S —
   but *only* when directly adjacent to something just stripped.

A sequence is never modified without the removal being recorded in
`stripped_elements`. Rule 3 is deliberately narrow: **internal GS-repeat
linkers are structural and are never stripped.** `IL10NB5` is a VH-linker-VH
construct whose two domains are joined by `GGGGS` repeats; those repeats are
part of the molecule and are retained. The smoke test covers this case
explicitly, along with a sequence whose own composition contains `HH` and an
internal `HHHHHH` that is not a tag.

If stripping would leave fewer than 5 residues it is abandoned, the sequence
is kept as given, and the attempted removal is recorded in
`stripped_elements`.

## Validation

The input is pre-cleaned, but it is validated on read rather than trusted.
Each sequence must contain only the 20 standard amino acids, and its actual
length must match the declared `sequence_length`. A failure names the binder
and the offending character and position. Failures are written to
`rejected_sequences.tsv` with their reason, valid rows are still processed,
and the run then exits non-zero so the failure cannot pass unnoticed.
Duplicate `binder_id` values are rejected the same way.

## Tiers

This is the sequence-only tier. Properties are computed by small
single-purpose functions assembled by a named *stage*; stages are registered
in `physchem_props.STAGES`. A structural tier (SASA, secondary structure,
structure-corrected solubility) and an ESM2 embedding tier are planned and
will register as additional stages contributing their own columns. Neither is
implemented, and there are no stubs for them — the stage list reflects only
what actually ran, and `physchem_provenance.txt` records which stages those
were.
