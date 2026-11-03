# AlphaFold 3 Server workflow

## Decision and supported workflow

This project uses the official [AlphaFold Server](https://alphafoldserver.com/) for predictions. It avoids installing AlphaFold 3, model parameters, or genetic databases locally, which also avoids the workstation's 12 GB GPU limit. This repository does not submit jobs to a service endpoint or automate browser interaction. The supported reproducible flow is to prepare JSON jobs locally, import/submit them in the Server interface, then download each result archive from the job page.

Google DeepMind documents an AlphaFold Server JSON dialect and importing jobs through the Server interface in its [server guide](https://github.com/google-deepmind/alphafold/blob/main/server/README.md). Protein entities use a sequence of the 20 standard amino acids and a copy count. The server input does not provide a PDB chain-ID field. Entity order is therefore recorded in our manifest and must be checked against chain IDs in downloaded structures. Do not use the separate local AlphaFold 3 JSON dialect (`alphafold3`).

The Server accepts standard one-letter amino acid codes only; B, J, O, U, X, gaps, and nonstandard residues are unsupported. The validator rejects these rather than substituting, trimming, or translating residues. The Server's current use terms and eligibility/usage conditions apply; review the [official FAQ](https://alphafoldserver.com/faq) and [output terms](https://alphafoldserver.com/output-terms) before submission. The output terms restrict certain automated interaction-prediction uses and ML models similar to AlphaFold; confirm the planned affinity-analysis use fits the accepted terms before submitting or using outputs. Do not upload confidential or sensitive sequences without authorization.

## Project findings

The repository has `scripts/select_complexes.py` for the 20-complex selection and `scripts/prepare_bench20_af3.py` to build Server jobs, a FASTA and native-chain mapping from the selected crystal structures. The index contains affinity text and chain labels in comments. `pdbbind/` is excluded by `.gitignore`; retain the selected structure IDs, native chain IDs and sequence-source notes.

`scripts/prepare_alphafold_server_inputs.py` consumes a CSV with one row per input chain:

```csv
complex_id,crystal_structure,chain_id,sequence
1abc,pdbbind/1abc/1abc_protein.cif,A,ACDEFGHIKLMNPQRSTVWY
1abc,pdbbind/1abc/1abc_protein.cif,B,ACDEFGHIKLMNPQRSTVWY
```

The example is schematic only; it is not a real PDBBind sequence. `crystal_structure` records the original file path/name used for the sequence extraction. Keep each chain's full selected sequence as supplied by the chosen extraction procedure. The script does not infer boundaries from the index or silently remove unresolved residues. It checks required fields, minimum two chains per complex, duplicate complex/chain keys, invalid amino-acid characters and inconsistent source-structure references. Identical sequences in distinct chain IDs are retained because they may represent genuine homomers.

Run from the repository root:

```bash
python scripts/prepare_alphafold_server_inputs.py chain_sequences.csv
```

The command fails before writing jobs if validation fails. Existing job JSON or manifest files are never overwritten. Add the resulting directory to the project's normal backup/versioning process as appropriate; generated prediction and structure files may be large and should not be committed by default.

## Output layout and manual procedure

The preparation command creates:

```text
alphafold3_server/
  inputs/       one Server-dialect JSON job per complex (each file is a one-job JSON list); submission_manifest.csv
  predictions/  untouched downloaded ZIP archives and extracted prediction files
  metadata/     submission records and any copied JSON confidence outputs
  logs/         local validation/submission notes
```

The 20-complex jobs are under `alphafold3_server/jobs/bench20_ready/`; use only that folder. Earlier Phase A drafts and the 10 pilot10 job JSONs are archived under `archive/superseded_phaseA/` and `archive/pilot10/af3_jobs/`; do not submit them. The sequence FASTA and native-to-A/B mapping are `alphafold3_server/bench20_sequences.fasta` and `alphafold3_server/bench20_chain_mapping.tsv`. The selection was generated with `python scripts/select_complexes.py`, then the AF3 artifacts with `python scripts/prepare_bench20_af3.py`. The selector uses a fixed seed and records its X-ray method checks in `bench20/experimental_methods.tsv`; the job generator refuses to overwrite existing JSON files.

1. Select the desired PDBBind complexes and extract the intended chain sequences from their crystal structures. Record exact chain identifiers and source structure in the input CSV. Resolve whether a chain should be a complete biological construct or only a crystallized fragment before extraction; the preparation tool cannot decide this.
2. Generate jobs with the command above. Review every sequence and mapping in the manifest. Server entity order follows the CSV row order for that complex. The Server does not preserve input PDB chain IDs as an explicit input setting; confirm the produced mmCIF chain IDs after download.
3. Open the official Server website, import each corresponding JSON job using its documented interface, review the displayed entities, and submit manually. Record the visible job identifier and submission date in the manifest. No public automated submission API is assumed here.
4. When the result is ready, use the job page's Download control to retrieve its ZIP archive. Keep the original archive unchanged under `predictions/`; extract a working copy there, and record its archive filename and any job/result identifiers in the manifest or a separate metadata record. Do not rename the original downloaded files.
5. Preserve all per-sample structures and JSON output files. The Server returns five samples per seed (an empty `modelSeeds` list requests a random seed); download archives include mmCIF structures and confidence/output JSON, and may include MSA/template artifacts. Consult the current FAQ for exact filenames because these can evolve.

## Confidence and reproducibility

The Server documents pLDDT (local confidence), PAE (relative-position error), pTM and ipTM (complex/interface confidence), plus a ranking score for ranking samples. The per-sample summary confidence JSON contains scalar keys `ptm`, `iptm`, and `ranking_score`. Record these values and the source filenames in `metadata/af3_confidence.tsv`; the header-only file is provided until predictions are downloaded. Higher ipTM generally indicates higher interface confidence; these are model confidence estimates, not experimental binding-affinity measurements. Retain the original JSON for every structure, not only the top-ranked sample. See [Server FAQ: confidence and outputs](https://alphafoldserver.com/faq) and the [official output terms](https://alphafoldserver.com/output-terms).

The current official FAQ states a limit of 30 submitted jobs per day and permits uploading JSON files containing up to 100 jobs. Thus 20 individual jobs fit within the published 30-job daily limit if the account has at least 20 jobs of quota remaining. Limits and terms can change; check the [live FAQ](https://alphafoldserver.com/faq) immediately before submission.

## Bench20 job state and working directory

The bench20 JSON files are ready under `alphafold3_server/jobs/bench20_ready/`. They have not been submitted. Follow [bench20/PHASE_B.md](bench20/PHASE_B.md) for manual submission/download instructions. Do not infer that downloaded predictions exist from the prepared jobs.

The analysis scripts read and write many filenames relative to the process working directory. For the bench20 run, use one dedicated run directory and pass absolute paths for list/model inputs. For example, from the repository root:

```bash
mkdir -p bench20/work alphafold3_server/predictions/raw/AF3_models
cd bench20/work
export PREFIX=bench20
export COMP_LIST="$OLDPWD/comp_list"
export CLONE_LIST="$OLDPWD/clone_list"
export AF3_MODEL_DIR="$OLDPWD/alphafold3_server/predictions/raw/AF3_models"
export RELAXED_MODEL_DIR="$PWD"
export CHAIN_MAPPING="$OLDPWD/alphafold3_server/bench20_chain_mapping.tsv"
export EXPT_FILE="$OLDPWD/expt_kd.csv"
```

After extracting each archive into `alphafold3_server/predictions/raw/<id>/`, convert the five model CIFs from `alphafold3_server/predictions/raw/` by running `python <repo>/scripts/cif_to_pdb.py` with that directory as the working directory. The converter writes to `AF3_models/` and requires that output directory to exist; it has not yet been validated against a real Server archive. Check the output names, chain IDs, and sequences before using them. Then use the dedicated `bench20/work` directory as above.

Run the relaxation, InterfaceAnalyzer, DockQ, and PLIP shell wrappers from this working directory (or pass the 20 IDs explicitly). The relaxed PDBs, scorefiles, PLIP files, and parser CSVs are written/read here; DockQ JSON/log output goes to the repository-level `dockq_results/` by default. Run `plip_parser_all.py` and `ifmetric_parser_all.py` here so they see their XML/score files. Run `ifmetric_addl.py`, `get_int_score_entropy.py`, and `merge_corr.py` here after the preceding CSVs exist. `PREFIX` names the generated tables, while `COMP_LIST`, `CLONE_LIST`, `AF3_MODEL_DIR`, `RELAXED_MODEL_DIR`, `CHAIN_MAPPING`, and `EXPT_FILE` avoid relying on the scripts' legacy relative-path defaults. Outputs are not silently moved between directories.

The shell wrappers accept complex IDs as positional arguments, which avoids reading `CLONE_LIST`. `prodigy_all.py`, `plip_parser_all.py`, and `ifmetric_parser_all.py` use `COMP_LIST`; set it explicitly as above. `merge_corr.py` also accepts `INTF_FILE`, `IFM_FILE`, and `PROD_FILE` when input tables are stored outside the current directory. The PLIP and Prodigy executable paths and the supplied GCC 7 Rosetta binary paths are fixed to the verified paths on this workstation.

For reproducibility, retain the exact input CSV, generated JSON, manifest, downloaded archive, submission date/job ID, and relevant output JSON. Server-side model/database versions and defaults can change; include retrieval date and any visible settings in the metadata. Bench20 jobs specify one fixed seed. Server-side model/database versions may still change, so exact prediction reproduction is not guaranteed.

## Downstream compatibility assessment

The downloaded structural predictions are mmCIF. The repository's `scripts/cif_to_pdb.py` reads local `????/*.cif` and writes PDBs under `AF3_models/`; it expects that directory to exist and assumes a directory/file naming layout that differs from the per-job download layout. It should not be run directly on the raw archive layout without organizing a working copy. It uses Biopython's `MMCIFParser`/`PDBIO`; conversion may lose mmCIF-only metadata and PDB has chain/residue limits, so retain and cite the original mmCIF. The bench20 wrappers expect converted filenames `<id>_0.pdb` through `<id>_4.pdb`; actual Server archives still need a representative conversion check.

The current wrappers reveal additional integration requirements:

| Tool/script | Observed input assumption | Integration action |
|---|---|---|
| DockQ (`dockq_all.sh`) | Reads original `<id>_<0..4>.pdb` and relaxed `<id>_<0..4>_0002.pdb`; native from `pdbbind/P-P/<id>_complex.pdb`. | Writes one JSON result per original and relaxed sample under `dockq_results/`. Requires two structures with a valid chain/residue correspondence. Tested on a PDBBind control, not an AF3 prediction. |
| `relax_parallel.sh` | Reads five PDBs from `AF3_MODEL_DIR`; defaults to `../AF3/AF3_models`. | Uses the verified GCC 7 FastRelax binary and renames Rosetta `_0001.pdb` output to the `_0002.pdb` convention consumed downstream. One-model run passed. |
| `ifanalyze_loop.sh` | Reads relaxed `<id>_<0..4>_0002.pdb`, assumes the two chains form the interface. | Uses the verified InterfaceAnalyzer binary. Writes `<id>_<sample>_ifmetric.sc`, including the five fields parsed by `ifmetric_parser_all.py`. One-model run and parser passed. |
| `ifmetric_parser_all.py` | Reads `<id>_<sample>_ifmetric.sc`. | Output filename and SCORE header were verified with the installed InterfaceAnalyzer. Set `COMP_LIST` to the absolute bench20 list path; score files are read from the working directory. |
| PLIP scripts | PLIP wrapper reads each relaxed `_0002.pdb`, selects A/B, writes `<id>_<sample>.xml` and `.txt`. | Updated executable path; one-model wrapper and five-sample parser smoke checks passed. Confirm prediction chain IDs are A/B. |
| `prodigy_all.py` | Reads five relaxed `_0002.pdb` files with selection A B. | Updated executable path; five-input script run and CSV output passed using the same representative structure. Set `COMP_LIST` explicitly; PDB files are read from the working directory. |

For DockQ comparisons, the prediction and crystal reference must represent corresponding chains and residue numbering; server-assigned IDs are not guaranteed to equal PDBBind IDs. Keep a reviewed mapping and make any chain renumbering on copies only. The wrappers were updated and smoke-tested on a small PDBBind complex. The CIF converter was not tested on an AlphaFold Server download, and no prediction was submitted or downloaded, so actual AF3 chain mapping, confidence extraction, and prediction DockQ remain pending.

## Remaining manual steps / open points

- Select complexes, decide biologically appropriate sequence boundaries, and create `chain_sequences.csv` from the actual PDBBind structures.
- Review the generated jobs and chain mappings, then submit each job through the official Server interface.
- Download and archive the results; fill submission metadata and verify actual output chain IDs.
- Verify CIF-to-PDB conversion and chain mapping on a representative real Server download before a benchmark run.
- Populate the confidence TSV with per-sample pTM/ipTM values from the downloaded summary JSON files.
- Before a bench20 run, follow the dedicated working-directory and environment-variable recipe above; verify the chosen model directory and chain mapping against each downloaded archive.
