# Environment Setup and Verification Report

This report records environment checks, data-presence checks, and the Phase A crystal PRODIGY baseline. No AlphaFold Server or relaxed-model pipeline stage has been run.

## Environment checks

| Tool or check | Version / path | Result |
|---|---|---|
| Conda environment `plip` | Python 3.9.13 | PASS |
| PLIP | 3.0.1; `/home/wetlab/miniforge3/envs/plip/bin/plip` | PASS; `plip --help` exits 0 |
| Open Babel | 3.1.1; `from openbabel import pybel` | PASS; import succeeds |
| Conda environment `ppi-affinity` | Python 3.11.18 | PASS |
| prodigy-prot | 2.2.4; `/home/wetlab/miniforge3/envs/ppi-affinity/bin/prodigy` | PASS; `prodigy --help` exits 0 |
| DockQ | 2.1.3; `/home/wetlab/miniforge3/envs/ppi-affinity/bin/DockQ` | PASS; `DockQ --help` exits 0 |
| Biopython | 1.88 | PASS; import succeeds |
| scikit-learn | 1.9.1 | PASS; import succeeds |
| pandas | 3.0.6 | PASS; import succeeds |
| scipy | 1.17.1 | PASS; import succeeds |
| numpy | 1.26.4 | PASS; import succeeds |
| Rosetta `relax.default.linuxgccrelease` | Requested `.default` binary name | **FAIL**; not present in the installation |
| Rosetta `rosetta_scripts.default.linuxgccrelease` | Requested `.default` binary name | **FAIL**; not present in the installation |
| Installed Rosetta FastRelax | `/home/wetlab/tools/rosetta.binary.ubuntu.release-408/main/source/build/src/release/linux/5.4/64/x86/gcc/7/static/relax.static.linuxgccrelease` | PASS; present, executable, `-help` exits 0 |
| Installed RosettaScripts | `/home/wetlab/tools/rosetta.binary.ubuntu.release-408/main/source/build/src/release/linux/5.4/64/x86/gcc/7/static/rosetta_scripts.static.linuxgccrelease` | PASS; present, executable, `-help` exits 0 |
| Installed Rosetta InterfaceAnalyzer | `/home/wetlab/tools/rosetta.binary.ubuntu.release-408/main/source/build/src/release/linux/5.4/64/x86/gcc/7/static/InterfaceAnalyzer.static.linuxgccrelease` | PASS; present, executable, `-help` exits 0 |

The supplied Rosetta build uses `.static.linuxgccrelease` names. The wrappers currently point to the installed static FastRelax and InterfaceAnalyzer binaries. The `.default` aliases requested above are unavailable; no installation or compilation was performed.

## Physchem track dependencies

The property source is [`scripts/physiochemical.tsv`](scripts/physiochemical.tsv), which contains 19 parameter rows although the initial request refers to 18. The README table follows the supplied names and tool suggestions. Cysteine number/arrangement and aggregation propensity combine sequence-only and structure-dependent work, so they are marked Tier 1 + 2.

| Tool / dependency | ppi-affinity status | Version | Notes |
|---|---|---:|---|
| Biopython (`Bio`) | Installed | 1.88 | Available for sequence descriptors and structure parsing |
| `peptides` | **Missing** | — | Python module not found in `ppi-affinity` |
| FreeSASA (`freesasa`) | Installed | 2.2.1 | Python module available for SASA calculations |
| `mkdssp` | **Missing** | — | Executable not on PATH; needed for the planned DSSP-based secondary-structure calculation |
| CamSol Intrinsic | Not installed / external service | — | Sequence-based web service is available after registration at the [official CamSol site](https://www-cohsoftware.ch.cam.ac.uk/). The official material says most algorithms can run locally and licensing is available, but a public standalone CLI for CamSol Intrinsic was not confirmed; check access/licensing before use. See the [CamSol method page](https://www-vendruscolo.ch.cam.ac.uk/camsolmethod.html) and [Cambridge Enterprise listing](https://licensing.enterprise.cam.ac.uk/product/camsol-method/). |
| AGGRESCAN3D | Missing in `ppi-affinity` and not on PATH | — | A standalone package with a command-line interface is available from the [official AGGRESCAN3D repository](https://bitbucket.org/lcbio/aggrescan3d/). The A3D 2.0 server also supports REST-based command-line operation, documented in the [A3D 2.0 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6602499/). No installation was performed. |

No analysis code or dependencies were installed for this track.

## Data and ignored files

- `pdbbind/` exists with `index/` and `P-P/`; `P-P/` contains 2,798 `*_complex.pdb` files.
- `pdbbind/index/INDEX_general_PP.2020R1.lst` has 2,805 lines including comments/header. The PDBBind directory is git-ignored and has no tracked files.
- `scripts/INDEX_general_PP.2020R1.txt` is present.
- Root `comp_list` and `clone_list` each contain 20 lines; `expt_kd.csv` is present. These are ignored by `.gitignore`.
- Structure, score, log, and CSV outputs are also git-ignored. PDBBind and AlphaFold Server prediction archives must be downloaded separately; generated outputs must be produced locally. The `bench20` selection and crystal PRODIGY baseline are now present; AF3 jobs are prepared but not submitted. Baseline rank correlations on 20 complexes: Spearman rho = 0.0925, Kendall tau = 0.0686.

## Setup status

| Component | Status | Current evidence / remaining work |
|---|---|---|
| `plip` environment | SET UP | Python, PLIP and Open Babel checks pass |
| `ppi-affinity` environment | SET UP | Prodigy, DockQ and scientific-package imports pass |
| Rosetta binaries | SET UP WITH NAME DIFFERENCE | Installed static FastRelax and RosettaScripts binaries start with `-help`; requested `.default` names are absent |
| PDBBind files | PRESENT, IGNORED | 2,798 PP structures and index file are present locally |
| AlphaFold Server stage | SET UP, NOT RUN | Input/job files and instructions are present; predictions and confidence outputs have not been downloaded |
| Analysis scripts | PHASE A COMPLETE | 20-complex selection and crystal PRODIGY baseline run; AF3/relaxed-model analysis remains pending |
| DockQ comparison | SET UP, NOT RUN | `scripts/dockq_all.sh` and DockQ are present; real prediction/crystal comparisons remain to be run |
| Machine-learning and extra validation stages | PLANNED | Not implemented in the current script set |

## Next steps

1. Review `bench20/selection.tsv`, `alphafold3_server/bench20_chain_mapping.tsv`, and the experimental Kd table.
2. Submit the 20 prepared sequences through AlphaFold Server and download all five CIF models and confidence JSON files per complex, after reviewing output-use terms.
3. Review chain IDs, convert the downloaded models, and then begin a controlled real-data pipeline run.
