#!/usr/bin/env bash
#SBATCH --job-name=PETKLS_GraphProps
#SBATCH --output=PETKLS_GraphProps_%j.out
#SBATCH --error=PETKLS_GraphProps_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  sbatch run_graph_properties.sh INPUT_ROOT OUTPUT_ROOT [LABELS_FILE] [EXTRA_ARGS...]

Example:
  sbatch run_graph_properties.sh \
    /data/pet_connectomes_kls \
    /data/pet_kls_graph_properties \
    /data/Schaefer2018_100Parcels_17Networks_labels.csv

Optional environment variables:
  PETKLS_PYTHON      Python executable (default: python3)
  PETKLS_MODULE      Environment module to load, e.g. Anaconda3
  PETKLS_CONDA_ENV   Existing conda environment to activate

The environment must contain numpy, pandas, and bctpy. The default run keeps
all positive KLS edges. Extra Python options can be appended after LABELS_FILE,
for example: --threshold-mode quantile --threshold 0.3 --overwrite
EOF
}

if (( $# < 2 )); then
    usage >&2
    exit 64
fi

INPUT_ROOT=$1
OUTPUT_ROOT=$2
LABELS_FILE=${3:-}
if (( $# >= 3 )); then
    shift 3
else
    shift 2
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHON_BIN=${PETKLS_PYTHON:-python3}

if [[ -n ${PETKLS_MODULE:-} ]]; then
    if ! command -v module >/dev/null 2>&1; then
        echo "ERROR: PETKLS_MODULE is set but the module command is unavailable." >&2
        exit 69
    fi
    module load "$PETKLS_MODULE"
fi

if [[ -n ${PETKLS_CONDA_ENV:-} ]]; then
    if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: PETKLS_CONDA_ENV is set but conda is unavailable." >&2
        exit 69
    fi
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$PETKLS_CONDA_ENV"
    PYTHON_BIN=${PETKLS_PYTHON:-python}
fi

if ! "$PYTHON_BIN" -c 'import numpy, pandas, bct' >/dev/null 2>&1; then
    echo "ERROR: Python dependencies are missing. Install numpy, pandas, and bctpy." >&2
    exit 69
fi

COMMAND=(
    "$PYTHON_BIN"
    "$SCRIPT_DIR/compute_graph_properties.py"
    --input-root "$INPUT_ROOT"
    --output-root "$OUTPUT_ROOT"
)
if [[ -n "$LABELS_FILE" ]]; then
    COMMAND+=(--labels-file "$LABELS_FILE")
fi
COMMAND+=("$@")

printf 'Running:'
printf ' %q' "${COMMAND[@]}"
printf '\n'
"${COMMAND[@]}"
