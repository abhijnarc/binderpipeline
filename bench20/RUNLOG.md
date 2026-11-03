# Bench20 run log

TaskCreate is unavailable in this session; tracking the requested task checklist here.

- [x] Select 20 complexes and write selection/list/Kd files.
- [x] Update downstream scripts and run directory/path handling.
- [x] Prepare AF3 Server jobs, FASTA, chain mapping, and manual instructions.
- [x] Run crystal PRODIGY baseline and report rank correlations.
- [x] Stop at Phase A pending user confirmation for AF3 downloads.

Do not edit physicochemical-property work or install software. Do not modify `pdbbind/`.

## Phase A completion (2026-10-30)

- [x] Selected exactly 20 structures from the local PDBBind index and atom records; fixed seed `20261030`.
- [x] Wrote `bench20/selection.tsv`, root `comp_list`/`clone_list`, and root `expt_kd.csv`.
- [x] Updated requested analysis scripts, PLIP inter-chain filtering/empty handling, PREFIX defaults, and root-anchored list paths.
- [x] Prepared 20 one-seed protein-only Server jobs, FASTA, and native/A-B chain mapping. Authoritative job folder: `alphafold3_server/jobs/bench20_ready/`.
- [x] Ran crystal PRODIGY once per complex at 4.5 A and wrote `bench20/prodigy_crystal.csv`.
- [x] Crystal baseline: n=20, Spearman rho=0.0925, Kendall tau=0.0686.
- [ ] Phase B not started; waiting for user confirmation after Server downloads.

The AlphaFold job draft folders from earlier selection passes were retained (not deleted) and are now archived under `archive/superseded_phaseA/`, together with the PRODIGY baseline from a superseded selection (`prodigy_crystal_superseded_selection.csv`, formerly root `bench20_prodigy_crystal.csv`); the 10 pilot10 JSONs are in `archive/pilot10/af3_jobs/`. Submit only `bench20_ready/`. The final selection includes 31 chain rows with a SEQRES/ATOM length difference and 9 with matching lengths; details are in `alphafold3_server/bench20_chain_mapping.tsv`. The selector checks for close homologs at >=70% identity over >=70% coverage, unique index parent-PDF references, and X-ray method metadata from the RCSB entry API. No selected chain pair exceeded the homology cutoff.
