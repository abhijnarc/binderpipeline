# Phase B handoff: AlphaFold Server downloads

Phase B is paused until the ten result archives are downloaded and the user explicitly asks to resume. No AlphaFold Server jobs have been submitted by this workflow.

## Before submission

Review Google's [AlphaFold Server FAQ](https://alphafoldserver.com/faq), [terms of service](https://alphafoldserver.com/terms-of-service), and [output terms](https://alphafoldserver.com/output-terms). The FAQ currently lists 30 submitted jobs per day, so these ten jobs fit in one day if at least ten jobs remain in the account's daily quota. The service supports importing up to 100 draft jobs per JSON file; the provided files each contain one job.

## Submit these jobs

Upload each JSON file using the Server's **Upload JSON** control, check the displayed chain sequences/order, and submit:

```text
alphafold3_server/jobs/2ptc.json
alphafold3_server/jobs/1kac.json
alphafold3_server/jobs/1ay7.json
alphafold3_server/jobs/2b0z.json
alphafold3_server/jobs/2c0l.json
alphafold3_server/jobs/2hle.json
alphafold3_server/jobs/5gjk.json
alphafold3_server/jobs/5dob.json
alphafold3_server/jobs/6fub.json
alphafold3_server/jobs/1op9.json
```

Each job has two protein entities, one random seed, and no ligand/ion entities. The first sequence maps to AF3 chain A and the second to chain B. Compare the displayed sequences with `alphafold3_server/pilot10_sequences.fasta` and the native chain IDs with `alphafold3_server/pilot10_chain_mapping.tsv` before submitting.

## Download and organize

For each completed job, use the result page's **Download** control. Preserve the downloaded ZIP's original filename and all filenames inside it. Copy each archive into a complex-specific folder and extract into a matching raw folder. Example for 2ptc (repeat with each ID):

```bash
mkdir -p alphafold3_server/predictions/archives/2ptc
mkdir -p alphafold3_server/predictions/raw/2ptc
cp -n /path/to/the/downloaded/2ptc-result.zip alphafold3_server/predictions/archives/2ptc/
unzip -n alphafold3_server/predictions/archives/2ptc/*.zip -d alphafold3_server/predictions/raw/2ptc/
```

Replace the example source path with the downloaded archive path; `cp -n` and `unzip -n` avoid overwriting existing files. Repeat for all ten complexes. Record the archive's actual filename, Server job ID, and download date in the confidence table or run log.

Before asking to resume, confirm for every complex that the raw directory contains five model mmCIF files and five matching `*_summary_confidences.json` files. Template mmCIFs may also be present in a nested templates directory; count only model structures. Check that each model has chains A and B in the submitted sequence order. Do not rename the archive members.

## Record pTM/ipTM

For each of the five model summary-confidence JSON files, copy its scalar `ptm`, `iptm`, and `ranking_score` values into `alphafold3_server/metadata/af3_confidence.tsv`. Use one row per prediction model:

```text
complex_id  model_index  structure_file  confidence_summary_file  ptm  iptm  ranking_score  server_job_id  download_date
```

The file already contains these tab-separated column headers. Keep `model_index` aligned with the model number in the original archive filename, and use paths relative to the repository root. Retain each original summary JSON unchanged.

## Resume Phase B

After the ten archives are extracted, the output chain IDs are checked, and the confidence TSV has 50 model rows, reply in this chat with:

> Phase B ready: all 10 AlphaFold Server archives are downloaded and organized, chains A/B are verified, and `af3_confidence.tsv` is populated. Proceed with Phase B.

Do not start Phase B before that confirmation. The next work will first validate archive counts and one complex's conversion, then continue through the five-model DockQ, relaxation, interface, PLIP, Prodigy, and correlation stages.
