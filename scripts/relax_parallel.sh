#!/usr/bin/env bash
set -euo pipefail

# Inputs are the five PDBs converted from original AlphaFold Server CIFs:
# <AF3_MODEL_DIR>/<complex>_<0..4>.pdb. Outputs retain *_0002.pdb names.
ROSETTA_RELAX=/home/wetlab/tools/rosetta.binary.ubuntu.release-408/main/source/build/src/release/linux/5.4/64/x86/gcc/7/static/relax.static.linuxgccrelease
ROSETTA_DB=/home/wetlab/tools/rosetta.binary.ubuntu.release-408/main/database
AF3_MODEL_DIR=${AF3_MODEL_DIR:-../AF3/AF3_models}
CLONE_LIST=${CLONE_LIST:-../clone_list}
MODEL_COUNT=${MODEL_COUNT:-5}

if (($#)); then
    complexes=("$@")
elif [[ -r "$CLONE_LIST" ]]; then
    mapfile -t complexes < <(sed '/^[[:space:]]*$/d' "$CLONE_LIST")
else
    echo "Missing complex IDs: pass them as arguments or provide $CLONE_LIST" >&2
    exit 2
fi

[[ -x "$ROSETTA_RELAX" ]] || { echo "FastRelax executable not found: $ROSETTA_RELAX" >&2; exit 2; }
[[ -d "$ROSETTA_DB" ]] || { echo "Rosetta database not found: $ROSETTA_DB" >&2; exit 2; }
[[ "$MODEL_COUNT" =~ ^[1-5]$ ]] || { echo "MODEL_COUNT must be between 1 and 5" >&2; exit 2; }

for pref in "${complexes[@]}"; do
    for ((i=0; i<MODEL_COUNT; i++)); do
        input="$AF3_MODEL_DIR/${pref}_${i}.pdb"
        output="${pref}_${i}_0002.pdb"
        score="${pref}_${i}_relax.sc"
        log="${pref}_${i}_relax.log"
        generated="${pref}_${i}_0001.pdb"
        [[ -s "$input" ]] || { echo "Missing prediction model: $input" >&2; exit 1; }
        [[ ! -e "$output" ]] || { echo "Refusing to overwrite relaxed model: $output" >&2; exit 1; }
        [[ ! -e "$score" ]] || { echo "Refusing to overwrite scorefile: $score" >&2; exit 1; }

        echo "FastRelax ${input} -> ${output}"
        "$ROSETTA_RELAX" \
            -database "$ROSETTA_DB" \
            -s "$input" \
            -relax:fast \
            -use_input_sc \
            -constrain_relax_to_start_coords \
            -relax:coord_constrain_sidechains \
            -ex1 -ex2 \
            -nstruct 1 \
            -out:file:scorefile "$score" >"$log" 2>&1
        [[ -s "$generated" ]] || { echo "FastRelax did not produce $generated; inspect $log" >&2; exit 1; }
        mv -- "$generated" "$output"
    done
done
