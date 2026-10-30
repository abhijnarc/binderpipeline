# Pilot 10 run log

Task tracker: the requested `TaskCreate` tool is not available in this session; this file tracks the phases and decisions.

## Phase A status

- [x] Read `ENVIRONMENT.md`, `ALPHAFOLD_SERVER.md`, and scripts.
- [x] Select and structurally validate ten PDBBind PP complexes.
- [x] Write root-level ID lists and experimental logKd table.
- [x] Prepare AF3 Server jobs, FASTA, and native-to-AF3 chain mapping.
- [x] Apply requested script fixes and verify their behavior.
- [x] Run crystal PRODIGY baseline and report correlations.
- [x] Write manual Phase B instructions; stop before any Server submission.

## Findings and selection

Selected: `2ptc`, `1kac`, `1ay7`, `2b0z`, `2c0l`, `2hle`, `5gjk`, `5dob`, `6fub`, `1op9`.

- All ten exist in both the PDBBind PP index and `scripts/INDEX_general_PP.2020R1.txt`; selected index lines match exactly.
- All have exactly two observed protein chains, each >50 residues, total length <=600, standard amino-acid residues only, and exact SEQRES/ATOM sequence matches.
- RCSB entry metadata confirms X-ray diffraction and resolution <=3.0 A for all ten. PDBBind complex files omit `EXPDTA`, so this method check used RCSB's official entry records.
- Parent PDF stems are unique. Selected types include enzyme-inhibitor, receptor-ligand signaling, viral complexes, electron transfer, plant immunity, chromatin remodeling, and a VHH-antigen complex.
- Experimental Kd spans 60 fM to 100 uM (9.22 log10 orders of magnitude).
- Borderline index comments: `2c0l` lists Kd=109nM with +/-34nM and `2hle` Kd=40nM with +/-20nM. Their index affinity fields themselves are exact `Kd=` values with no `<`, `>`, or `~`, so they pass the stated rule; the uncertainty is retained in `selection.tsv`. `2b0z` is the F82I cytochrome c complex, but no second member of that family/parent was selected.
- Atom/SEQRES discrepancies: none in selected rows.

## Inputs and baseline

- Root `comp_list` and `clone_list` each contain the ten lowercase IDs, one per line; `expt_kd.csv` is `complex,logKd`, with log10(Kd in M) in column 2.
- Created ten one-job AlphaFold Server JSON files under `alphafold3_server/jobs/`, a two-chain-per-complex FASTA, and `alphafold3_server/pilot10_chain_mapping.tsv`. Every job has `modelSeeds: []`, protein chains only, and order native chain 1 -> AF3 A, native chain 2 -> AF3 B.
- Sequence derivation used actual ATOM residues; SEQRES matched exactly for both selected chains. All residues are standard amino acids and no modified amino-acid residues were found.
- Crystal PRODIGY baseline: 10/10 completed using each native chain pair and the 4.5 A cutoff. Spearman rho=0.449850 (p=0.192076), Kendall tau=0.404520 (p=0.105998), n=10. This is a noisy pilot reference, not a scientific conclusion. Results and tool output are in ignored `pilot10/prodigy_crystal.csv` and `pilot10/prodigy_crystal.log`.

## Script changes and verification

- `prodigy_all.py`: one subprocess call per pose; Kd is kept in molar units and transformed with `np.log10` to match experimental log10(Kd in M); failures become NaN so means skip them instead of using 9999.
- The six requested analysis scripts use `PREFIX` from the environment, default `pilot10`, for prefixed input/output CSV names. `COMP_LIST` and `CLONE_LIST` can be overridden without changing their previous relative-path fallbacks.
- `plip_parser_all.py`: empty-interaction XMLs now emit zero rows for A/B across the four existing interaction categories. PLIP's `--chains "[['A'],['B']]"` divides A as receptor and B as ligand; the smoke XML contained only `reschain=A`, `reschain_lig=B`, so the wrapper counts cross-chain A-to-B contacts, not within-chain contacts.
- `merge_corr.py`: computes Spearman and Kendall with pairwise finite rows, includes n and p-values, flags |Spearman rho| >=0.3, and writes prefixed correlation/master CSV files.
- `ifmetric_addl.py`: aggregation definition retained; `columns[2:]` includes `dG_norm` and `unsat_density`, confirmed in smoke output.
- Smoke checks in `/tmp/pilot_script_checks`: empty PLIP XML -> eight zero-count A/B rows -> zero weighted score and entropy; prefix outputs; Prodigy failure -> NaN; InterfaceAnalyzer parser aggregation; merge correlation NaN handling. Interaction weights and entropy definition were not changed.
- `.gitignore` excludes pilot tables, root lists, AF3 job JSONs, FASTA, and mapping data. `pdbbind/` was not modified.

## Baseline PRODIGY statistics

- n=10 paired structures.
- Spearman rho=0.449850, p=0.192076.
- Kendall tau=0.404520, p=0.105998.
- Do not draw scientific conclusions from this pilot-sized sample.

## Remaining gate

The ten AF3 jobs have not been submitted. No model structures or confidence JSONs exist locally. Wait for the user to review the job inputs and terms, submit/download the archives, populate `af3_confidence.tsv`, and explicitly authorize Phase B.
