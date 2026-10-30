#!/usr/bin/env bash
set -euo pipefail

PLIP=/home/wetlab/miniforge3/envs/plip/bin/plip
PLIPDIR=${PLIPDIR:-./plip}
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

[[ -x "$PLIP" ]] || { echo "PLIP executable not found: $PLIP" >&2; exit 2; }
[[ "$MODEL_COUNT" =~ ^[1-5]$ ]] || { echo "MODEL_COUNT must be between 1 and 5" >&2; exit 2; }
mkdir -p "$PLIPDIR"

for pref in "${complexes[@]}"; do
    for ((i=0; i<MODEL_COUNT; i++)); do
        pdb="${pref}_${i}_0002.pdb"
        name="${pref}_${i}"
        [[ -s "$pdb" ]] || { echo "Missing relaxed model: $pdb" >&2; exit 1; }
        "$PLIP" -f "$pdb" --out "$PLIPDIR" --txt --xml --name "$name" --chains "[['A'],['B']]"
    done
done
