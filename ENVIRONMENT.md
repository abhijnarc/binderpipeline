# Environment Setup and Verification Report

This report records environment and file-presence checks only. The analysis pipeline has not been run on real data, and no pipeline component is claimed as tested on real data.

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

## Data and ignored files

- `pdbbind/` exists with `index/` and `P-P/`; `P-P/` contains 2,798 `*_complex.pdb` files.
- `pdbbind/index/INDEX_general_PP.2020R1.lst` has 2,805 lines including comments/header. The PDBBind directory is git-ignored and has no tracked files.
- `scripts/INDEX_general_PP.2020R1.txt` is present.
- Root `comp_list` and `clone_list` each contain 10 lines; `expt_kd.csv` is present. These are ignored by `.gitignore`.
- Structure, score, log, and CSV outputs are also git-ignored. PDBBind and AlphaFold Server prediction archives must be downloaded separately; generated outputs must be produced locally.

## Setup status

| Component | Status | Current evidence / remaining work |
|---|---|---|
| `plip` environment | SET UP | Python, PLIP and Open Babel checks pass |
| `ppi-affinity` environment | SET UP | Prodigy, DockQ and scientific-package imports pass |
| Rosetta binaries | SET UP WITH NAME DIFFERENCE | Installed static FastRelax and RosettaScripts binaries start with `-help`; requested `.default` names are absent |
| PDBBind files | PRESENT, IGNORED | 2,798 PP structures and index file are present locally |
| AlphaFold Server stage | SET UP, NOT RUN | Input/job files and instructions are present; predictions and confidence outputs have not been downloaded |
| Analysis scripts | PRESENT, NOT RUN ON REAL DATA | No pipeline execution is claimed in this report |
| DockQ comparison | SET UP, NOT RUN | `scripts/dockq_all.sh` and DockQ are present; real prediction/crystal comparisons remain to be run |
| Machine-learning and extra validation stages | PLANNED | Not implemented in the current script set |

## Next steps

1. Review the selected complex inputs, chain mapping, and experimental Kd table.
2. Submit the selected sequences through AlphaFold Server and download all five CIF models and confidence JSON files per complex.
3. Review chain IDs, convert the downloaded models, and then begin a controlled real-data pipeline run.
