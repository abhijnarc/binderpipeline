#!/usr/bin/env bash
set -euo pipefail

DOCKQ=/home/wetlab/miniforge3/envs/ppi-affinity/bin/DockQ
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
AF3_MODEL_DIR=${AF3_MODEL_DIR:-../AF3/AF3_models}
RELAXED_MODEL_DIR=${RELAXED_MODEL_DIR:-.}
NATIVE_DIR=${NATIVE_DIR:-$PROJECT_ROOT/pdbbind/P-P}
RESULT_DIR=${DOCKQ_RESULT_DIR:-$PROJECT_ROOT/dockq_results}
CHAIN_MAPPING=${CHAIN_MAPPING:-$PROJECT_ROOT/alphafold3_server/pilot10_chain_mapping.tsv}
MODEL_COUNT=${MODEL_COUNT:-5}
CLONE_LIST=${CLONE_LIST:-../clone_list}

if (($#)); then
    complexes=("$@")
elif [[ -r "$CLONE_LIST" ]]; then
    mapfile -t complexes < <(sed '/^[[:space:]]*$/d' "$CLONE_LIST")
else
    echo "Missing complex IDs: pass them as arguments or provide $CLONE_LIST" >&2
    exit 2
fi

[[ -x "$DOCKQ" ]] || { echo "DockQ executable not found: $DOCKQ" >&2; exit 2; }
[[ -r "$CHAIN_MAPPING" ]] || { echo "Chain mapping TSV not found: $CHAIN_MAPPING" >&2; exit 2; }
[[ "$MODEL_COUNT" =~ ^[1-5]$ ]] || { echo "MODEL_COUNT must be between 1 and 5" >&2; exit 2; }
mkdir -p "$RESULT_DIR/original" "$RESULT_DIR/relaxed"

for pref in "${complexes[@]}"; do
    native="$NATIVE_DIR/${pref}_complex.pdb"
    [[ -s "$native" ]] || { echo "Missing native crystal structure: $native" >&2; exit 1; }
    mapfile -t native_chains < <(awk -F '\t' -v id="$pref" 'NR > 1 && $1 == id {print $3}' "$CHAIN_MAPPING")
    [[ "${#native_chains[@]}" -eq 2 ]] || { echo "Expected two native chain mappings for $pref in $CHAIN_MAPPING" >&2; exit 1; }
    chain_mapping="AB:${native_chains[0]}${native_chains[1]}"
    for ((i=0; i<MODEL_COUNT; i++)); do
        original="$AF3_MODEL_DIR/${pref}_${i}.pdb"
        relaxed="$RELAXED_MODEL_DIR/${pref}_${i}_0002.pdb"
        original_result="$RESULT_DIR/original/${pref}_${i}.json"
        relaxed_result="$RESULT_DIR/relaxed/${pref}_${i}.json"
        [[ -s "$original" ]] || { echo "Missing original prediction: $original" >&2; exit 1; }
        [[ -s "$relaxed" ]] || { echo "Missing relaxed prediction: $relaxed" >&2; exit 1; }
        [[ ! -e "$original_result" && ! -e "$relaxed_result" ]] || {
            echo "Refusing to overwrite DockQ results for ${pref} model ${i}" >&2
            exit 1
        }
        "$DOCKQ" "$original" "$native" --mapping "$chain_mapping" --json "$original_result" >"${original_result%.json}.log" 2>&1
        "$DOCKQ" "$relaxed" "$native" --mapping "$chain_mapping" --json "$relaxed_result" >"${relaxed_result%.json}.log" 2>&1
    done
done
