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

## Phase B (2026-11-03)

- [x] AF3 Server downloads present for all 20 complexes in `alphafold3_server/predictions/pdbbind_<id>_bench20/` (5 model CIFs + 5 summary confidences each). Job requests match submitted sequences and seed; every model has chains A/B in submitted order; model_0 is top-ranked.
- [x] Working copies staged as `predictions/raw/<id>/<id>_0..4.cif`; `af3_confidence.tsv` populated (100 rows; Server job IDs not available in the downloads); 100 PDBs converted to `predictions/raw/AF3_models/`.
- Decision (option a): keep all 20 complexes. Flag low-confidence interfaces 1zvh (ipTM 0.24-0.27), 3qhy (0.15-0.18), and 5szj (~0.26); report correlations with and without them.
- Note: 3w8h's own PDB entry appeared among AF3 template hits (single-chain templates).
- [x] Smoke test, 1avx model 0, `MODEL_COUNT=1` (outputs kept in `bench20/work/smoke/`): FastRelax 1239 s, total_score -1288.3, chains A/B retained; InterfaceAnalyzer dG_separated -59.3, dSASA 2016, sc 0.70; PLIP 18 HB / 8 hydrophobic / 2 salt bridges; PRODIGY logKd -7.46 (expt -13.22); DockQ 0.762 original, 0.704 relaxed. Parsers, `ifmetric_addl.py`, and `get_int_score_entropy.py` all ran.
- Note: `plip_parser_all.py` must run inside `PLIPDIR` (default `./plip`) to see the XMLs; copy its `<PREFIX>_intsummary.csv` up to the work directory before `get_int_score_entropy.py`.
- [x] Full FastRelax: 100/100 models, 20 complexes in parallel (`xargs -P 20`, 5 models sequential each), logs in `bench20/work/relax_logs/`. 2026-11-03 18:5x to 2026-11-04 01:36; 4bkx (658 res) was rate-limiting at ~50 min/model.
- [x] InterfaceAnalyzer 100/100, PLIP 100/100, DockQ 200/200 (original + relaxed), parsers, `ifmetric_addl.py`, `get_int_score_entropy.py`, `prodigy_all.py`, and `merge_corr.py`; no stage failures. Stage logs in `bench20/work/stage_logs/`.
- DockQ (`bench20/work/bench20_dockq.csv`): mean 0.688 original, 0.657 relaxed. 16 complexes medium/high; incorrect: 1zvh, 3qhy, 5szj (low ipTM) and 5bnq (ipTM 0.85 but DockQ 0.03, fnat 0 in all five models; confirmed by an independent alignment-based contact check, so not a mapping artefact). Relaxation dropped 4bkx DockQ from ~0.92 to 0.56-0.69 in four of five models.
- Correlations written for three subsets: `bench20_correlations.csv` (n=20), `bench20_excl_lowiptm_correlations.csv` (n=17, option a), `bench20_excl_incorrect_correlations.csv` (n=16, also excludes 5bnq). Only HB passes |rho|>=0.3 in all three (rho -0.38/-0.45/-0.46, Spearman p 0.07-0.10); no feature reaches p<0.05 by Spearman, and none would survive multiple-testing correction.
- PRODIGY on relaxed AF3 models vs crystal: rho 0.70 between the two; vs experiment rho 0.162 (relaxed) and 0.093 (crystal), n=20 (`bench20_prodigy_crystal_vs_relaxed.csv`).
- [x] DockQ comparison plot: `bench20/work/bench20_dockq_comparison.png` / `.pdf` (per-complex dumbbell, AF3 vs relaxed, mean of 5 models with model range, CAPRI class boundaries, ipTM column).
