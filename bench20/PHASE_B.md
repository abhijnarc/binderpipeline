# Phase B: manual AlphaFold Server submission and download

Do not start Phase B analysis until the user confirms that the downloads are complete. The AF3 jobs below are prepared but have not been submitted.

## Before submitting

- The current AlphaFold Server FAQ publishes a limit of 30 submissions per day. These 20 jobs fit in one day if at least 20 submissions remain in the account quota. Check the live limit on the submission day.
- Review the Server Output Terms and obtain the supervisor/institution's interpretation before using predictions in this affinity-prediction workflow. The terms restrict output use in automated systems that predict protein binding or interaction, and separately restrict training structure-prediction-like ML models. Whether the proposed PPI Kd model is permitted should be resolved before submission.
- Use only the authoritative jobs in `alphafold3_server/jobs/bench20_ready/`; earlier Phase A drafts are archived under `archive/superseded_phaseA/` and the 10 pilot10 JSONs under `archive/pilot10/af3_jobs/`; do not submit either set. The jobs contain protein sequences only, one seed (`20261030`), and one copy of each partner. The first selected native chain is AF3 entity A; the second is entity B. The job JSONs are in `alphafold3_server/jobs/bench20_ready/`.

## Submit and download

1. Upload each JSON from `alphafold3_server/jobs/bench20_ready/` through the Server UI, review the two displayed sequences and submit. Record the Server job ID and submission date for every complex.
2. Download each completed job ZIP without renaming it. Keep originals in `alphafold3_server/predictions/archives/`.
3. Extract each archive to a separate staging directory under `alphafold3_server/predictions/extracted/<id>/`. Keep all five CIFs and all five per-sample confidence JSONs. Preserve other output files as well.
4. In a working copy, place the five structure files in `alphafold3_server/predictions/raw/<id>/` and name them `<id>_0.cif` through `<id>_4.cif` (copy/rename only the working copies; retain the archive and original extracted filenames). Identify each sample by its paired confidence JSON before assigning the index.
5. For each sample's summary confidence JSON, record `ptm`, `iptm`, and `ranking_score` in `alphafold3_server/metadata/af3_confidence.tsv`, alongside complex ID, model index, CIF filename, confidence filename, Server job ID, and download date. Keep the original JSONs.
6. Check the mmCIF chain IDs and sequence order against `alphafold3_server/bench20_chain_mapping.tsv`. Confirm the input partner mapped to A is the first chain in that table and B is the second. Record any discrepancy before running conversion or DockQ.

The converter expects four-character complex directories and CIF basenames whose final underscore-delimited portion is the sample index. From the raw staging directory, prepare its output folder and run it this way:

```bash
cd alphafold3_server/predictions/raw
mkdir -p AF3_models
python <repo>/scripts/cif_to_pdb.py
```

Check that it wrote 100 files named `<id>_0.pdb` through `<id>_4.pdb` for each of the 20 complexes, and inspect chain IDs and sequences. Keep mmCIF originals because conversion to PDB can lose metadata.

## Resume command after user confirmation

From the repository root, set paths for the dedicated working directory:

```bash
mkdir -p bench20/work
cd bench20/work
export PREFIX=bench20
export COMP_LIST="<repo>/comp_list"
export CLONE_LIST="<repo>/clone_list"
export EXPT_FILE="<repo>/expt_kd.csv"
export CHAIN_MAPPING="<repo>/alphafold3_server/bench20_chain_mapping.tsv"
export AF3_MODEL_DIR="<repo>/alphafold3_server/predictions/raw/AF3_models"
export RELAXED_MODEL_DIR="$PWD"
```

Replace `<repo>` with the absolute repository path. Smoke-test FastRelax on one complex and one pose first:

```bash
MODEL_COUNT=1 <repo>/scripts/relax_parallel.sh <complex_id>
```

Keep `MODEL_COUNT=1` exported (with `COMP_LIST` pointing at a one-ID list) to run the rest of the stages on that single pose; the shell wrappers and `plip_parser_all.py`, `ifmetric_parser_all.py`, and `prodigy_all.py` all honour it. Unset it (default 5) for the full run.

After inspecting its output, run all five models for all 20 complexes. The current wrapper runs sequentially (one process), which stays below the requested `cores - 2` concurrency cap. Then run DockQ, InterfaceAnalyzer, and PLIP wrappers; parse their outputs; run `prodigy_all.py`, `ifmetric_addl.py`, `get_int_score_entropy.py`, and `merge_corr.py`. Compare crystal and relaxed-model PRODIGY, create the requested DockQ comparison plot, and report stage completion per complex. Do not replace complexes that fail; stop and ask before proceeding if more than two fail any stage.

This completes only the manual handoff. Phase B analysis remains gated on the user's confirmation after the 20 archives and 100 CIF/confidence pairs have been downloaded and checked.
