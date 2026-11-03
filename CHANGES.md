# Changes from the original scripts

| File | What changed | Why |
|---|---|---|
| `scripts/get_int_score_entropy.py` | Reads/writes prefix-based filenames (default `bench20`); avoids division by zero and computes entropy only for positive fractions. | Supports separate pilot outputs and handles empty interaction totals without creating invalid entropy values. |
| `scripts/ifanalyze_loop.sh` | Uses the installed static InterfaceAnalyzer path, accepts IDs or a clone list, loops over models, and checks inputs/outputs. | Makes invocation configurable for this installation and reports missing inputs or outputs. |
| `scripts/ifmetric_addl.py` | Reads/writes prefix-based filenames; clarifies which derived metrics are included in aggregation. | Supports pilot-specific tables and documents the existing aggregation columns. |
| `scripts/ifmetric_parser_all.py` | Adds `PREFIX`, `COMP_LIST`, and `MODEL_COUNT` (1-5, default 5) environment overrides; `COMP_LIST` defaults to the root `comp_list`. | Allows the parser to use selected inputs, separate output names, and a one-pose smoke test. |
| `scripts/merge_corr.py` | Reworks table merging and correlation output; calculates pairwise-complete Spearman and Kendall statistics, sample counts, p-values, and an absolute Spearman threshold flag; uses prefix-based outputs. | Produces explicit per-feature correlation results while tolerating missing/non-numeric feature values. |
| `scripts/plip_loop.sh` | Uses the installed PLIP path, accepts IDs or a clone list, loops over models, quotes arguments, and checks inputs. | Makes the wrapper fit the current environment and process selected models with clearer failure handling. |
| `scripts/plip_parser_all.py` | Adds fixed interaction columns, zero-count rows for empty interaction sets, A<->B inter-chain filtering, and `PREFIX`/`COMP_LIST`/`MODEL_COUNT` overrides. | Keeps empty poses representable downstream, excludes intra-chain contacts, and allows selected inputs and a one-pose smoke test. |
| `scripts/prodigy_all.py` | Uses the installed Prodigy path; runs one argument-list subprocess per model; handles command/parse failures as missing values; reports `logKd` as log10 molar Kd; adds input/output and `MODEL_COUNT` overrides; adds a `--crystal` mode that runs once per native structure using the chains in `bench20/selection.tsv` and writes `<PREFIX>/prodigy_crystal.csv`. | Avoids shell-string execution and makes the output units and failure values explicit. |
| `scripts/relax_parallel.sh` | Replaces hardcoded complex IDs with command-line/list input; points to the installed static FastRelax binary; adds configurable input directory/model count and file checks; renames Rosetta `_0001.pdb` output to the `_0002.pdb` name expected downstream. | Makes the batch inputs configurable and aligns Rosetta output names with current consumers. |

## Additional files outside the original scripts commit

The following files are present in the working tree but were not part of the original scripts commit.

- Root documentation/configuration: `.gitignore`, `ALPHAFOLD_SERVER.md`, `ENVIRONMENT.md`, `README.md`, `CHANGES.md`.
- Added helpers under `scripts/`: `dockq_all.sh`, `prepare_alphafold_server_inputs.py` (pilot10 job generator), `select_complexes.py` (bench20 selection), `prepare_bench20_af3.py` (bench20 jobs, FASTA, and chain mapping), and `physiochemical.tsv`.
- Locally present, git-ignored root inputs: `comp_list`, `clone_list` (20 bench20 IDs each), and `expt_kd.csv` (`complex,logKd`, 20 rows).
- `bench20/`: `PHASE_B.md` and `RUNLOG.md` (tracked); `selection.tsv`, `experimental_methods.tsv`, and `prodigy_crystal.csv` (git-ignored).
- `alphafold3_server/`: `metadata/af3_confidence.tsv` (header only) and `.gitkeep` files under `inputs/`, `logs/`, `metadata/`, and `predictions/` (tracked); `bench20_sequences.fasta`, `bench20_chain_mapping.tsv`, and the 20 authoritative jobs in `jobs/bench20_ready/` (git-ignored).
- `archive/`: the superseded pilot10 run (`pilot10/` notes, tables, FASTA, mapping, and its 10 never-submitted jobs under `af3_jobs/`) and earlier Phase A job drafts plus a superseded-selection PRODIGY table under `superseded_phaseA/`. Only the Markdown notes are tracked.
- `pdbbind/` is locally present but ignored by Git; it is external dataset content, not part of the original scripts commit.

## Known issues not yet fixed

- **Chain names A/B are assumed.** `scripts/interface_AB.xml:14`, `scripts/plip_loop.sh:28`, and `scripts/prodigy_all.py:31` select chains using A/B. This matches AlphaFold Server output; `scripts/dockq_all.sh` maps A/B to the native chains via `alphafold3_server/bench20_chain_mapping.tsv`. Check the downloaded mmCIF chain IDs before analysis.
- **Rosetta output naming needs the wrapper.** Rosetta writes `<complex>_<model>_0001.pdb` (`scripts/relax_parallel.sh:32`); the wrapper renames it to `_0002.pdb` (lines 29 and 49) because the downstream consumers read `_0002.pdb` (`scripts/ifanalyze_loop.sh:25`, `scripts/plip_loop.sh:25`, `scripts/prodigy_all.py:25`). Running Rosetta outside this wrapper will not perform that rename.
- **Per-pose files are read and written in the current working directory.** The ID lists and `expt_kd.csv` default to repository-root paths, but relaxed PDBs, score files, PLIP XMLs (`PLIPDIR`, default `./plip`), and all `<PREFIX>_*.csv` tables are relative to the caller's directory. Run every Phase B stage from the same working directory (`bench20/work/`).
- **Converted-model folder default differs from the instructions.** `scripts/relax_parallel.sh:9` and `scripts/dockq_all.sh:6` default `AF3_MODEL_DIR` to `alphafold3_server/predictions/AF3_models`, while `bench20/PHASE_B.md` converts into `alphafold3_server/predictions/raw/AF3_models`. Export `AF3_MODEL_DIR` as PHASE_B shows.
- **Hardcoded complexes in `relax_parallel.sh` were removed.** The original committed script hardcoded `1a22` and `3mca` at original line 8. The current script accepts command-line IDs or reads `CLONE_LIST` at `scripts/relax_parallel.sh:13-20`.
