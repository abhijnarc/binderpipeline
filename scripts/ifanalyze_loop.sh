#!/usr/bin/env bash
set -euo pipefail

ROSETTA_IA=/home/wetlab/tools/rosetta.binary.ubuntu.release-408/main/source/build/src/release/linux/5.4/64/x86/gcc/7/static/InterfaceAnalyzer.static.linuxgccrelease
ROSETTA_DB=/home/wetlab/tools/rosetta.binary.ubuntu.release-408/main/database
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

[[ -x "$ROSETTA_IA" ]] || { echo "InterfaceAnalyzer executable not found: $ROSETTA_IA" >&2; exit 2; }
[[ -d "$ROSETTA_DB" ]] || { echo "Rosetta database not found: $ROSETTA_DB" >&2; exit 2; }
[[ "$MODEL_COUNT" =~ ^[1-5]$ ]] || { echo "MODEL_COUNT must be between 1 and 5" >&2; exit 2; }

for pref in "${complexes[@]}"; do
    for ((i=0; i<MODEL_COUNT; i++)); do
        pdb="${pref}_${i}_0002.pdb"
        score="${pref}_${i}_ifmetric.sc"
        log="${pref}_${i}_ifmetric.log"
        [[ -s "$pdb" ]] || { echo "Missing relaxed model: $pdb" >&2; exit 1; }
        [[ ! -e "$score" ]] || { echo "Refusing to overwrite scorefile: $score" >&2; exit 1; }
        echo "InterfaceAnalyzer ${pdb} -> ${score}"
        "$ROSETTA_IA" \
            -database "$ROSETTA_DB" \
            -s "$pdb" \
            -pack_separated true \
            -compute_packstat true \
            -compute_interface_sc true \
            -out:file:score_only "$score" >"$log" 2>&1
        [[ -s "$score" ]] || { echo "InterfaceAnalyzer did not produce $score; inspect $log" >&2; exit 1; }
    done
done
