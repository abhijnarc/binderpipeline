# Changes from the original scripts

| File | What changed | Why |
|---|---|---|
| `scripts/get_int_score_entropy.py` | Reads/writes prefix-based filenames; avoids division by zero and computes entropy only for positive fractions. | Supports separate pilot outputs and handles empty interaction totals without creating invalid entropy values. |
| `scripts/ifanalyze_loop.sh` | Uses the installed static InterfaceAnalyzer path, accepts IDs or a clone list, loops over models, and checks inputs/outputs. | Makes invocation configurable for this installation and reports missing inputs or outputs. |
| `scripts/ifmetric_addl.py` | Reads/writes prefix-based filenames; clarifies which derived metrics are included in aggregation. | Supports pilot-specific tables and documents the existing aggregation columns. |
| `scripts/ifmetric_parser_all.py` | Adds `PREFIX` and `COMP_LIST` environment overrides for filenames and input-list path. | Allows the parser to use selected inputs and separate output names. |
| `scripts/merge_corr.py` | Reworks table merging and correlation output; calculates pairwise-complete Spearman and Kendall statistics, sample counts, p-values, and an absolute Spearman threshold flag; uses prefix-based outputs. | Produces explicit per-feature correlation results while tolerating missing/non-numeric feature values. |
| `scripts/plip_loop.sh` | Uses the installed PLIP path, accepts IDs or a clone list, loops over models, quotes arguments, and checks inputs. | Makes the wrapper fit the current environment and process selected models with clearer failure handling. |
| `scripts/plip_parser_all.py` | Adds fixed interaction columns, zero-count rows for empty interaction sets, and `PREFIX`/`COMP_LIST` overrides. | Keeps empty poses representable downstream and allows pilot-specific inputs and outputs. |
| `scripts/prodigy_all.py` | Uses the installed Prodigy path; runs one argument-list subprocess per model; handles command/parse failures as missing values; reports `logKd` as log10 molar Kd; adds input/output overrides. | Avoids shell-string execution and makes the output units and failure values explicit. |
| `scripts/relax_parallel.sh` | Replaces hardcoded complex IDs with command-line/list input; points to the installed static FastRelax binary; adds configurable input directory/model count and file checks; renames Rosetta `_0001.pdb` output to the `_0002.pdb` name expected downstream. | Makes the batch inputs configurable and aligns Rosetta output names with current consumers. |

## Additional files outside the original scripts commit

The following files are present in the working tree but were not part of the original scripts commit.

- Root documentation/configuration: `.gitignore`, `ALPHAFOLD_SERVER.md`, `ENVIRONMENT.md`, `README.md`, `CHANGES.md`.
- Added helpers under `scripts/`: `dockq_all.sh`, `prepare_alphafold_server_inputs.py`.
- Locally present, git-ignored root inputs: `comp_list`, `clone_list`, `expt_kd.csv` (the two lists each have 10 lines).
- `alphafold3_server/`: `pilot10_sequences.fasta`, `pilot10_chain_mapping.tsv`, `metadata/af3_confidence.tsv`, ten files under `jobs/`, plus `.gitkeep` files under `inputs/`, `logs/`, `metadata/`, and `predictions/`.
- `pilot10/`: `PHASE_B.md`, `RUNLOG.md`, `selection.tsv`, `prodigy_crystal.csv`, and `prodigy_crystal.log`.
- `pdbbind/` is locally present but ignored by Git; it is external dataset content, not part of the original scripts commit.

## Known issues not yet fixed

- **Chain names A/B are assumed.** `scripts/interface_AB.xml:14`, `scripts/plip_loop.sh:27`, and `scripts/prodigy_all.py:25` select chains using A/B. AlphaFold Server output chain IDs and the crystal chain mapping still need review before analysis.
- **Rosetta output naming needs the wrapper.** Rosetta writes `<complex>_<model>_0001.pdb` at `scripts/relax_parallel.sh:31`; the wrapper renames it to `_0002.pdb` at lines 28 and 48 because the current downstream wrappers read `_0002.pdb` (`scripts/ifanalyze_loop.sh:24`, `scripts/plip_loop.sh:24`, `scripts/prodigy_all.py:19`). Running Rosetta outside this wrapper will not perform that rename.
- **Relative input defaults depend on the caller's working directory.** Examples: `scripts/relax_parallel.sh:8-9`, `scripts/ifanalyze_loop.sh:6`, `scripts/plip_loop.sh:6`, `scripts/dockq_all.sh:6,12`, `scripts/prodigy_all.py:55`, `scripts/plip_parser_all.py:180`, and `scripts/ifmetric_parser_all.py:47`. Environment overrides or explicit IDs are available for some of these inputs.
- **Hardcoded complexes in `relax_parallel.sh` were removed.** The original committed script hardcoded `1a22` and `3mca` at original line 8. The current script accepts command-line IDs or reads `CLONE_LIST` at `scripts/relax_parallel.sh:12-18`.
