#!/usr/bin/env bash
#SBATCH --job-name=petkls-graphprops
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=petkls-graphprops_%j.out
#SBATCH --error=petkls-graphprops_%j.err

# SLURM wrapper for recompute_graph_properties.py.
#
# Usage:
#   sbatch run_graph_properties.sh INPUT_CSV OUTPUT_DIRECTORY [PYTHON_OPTIONS...]
#
# Optional environment variables:
#   PETKLS_PYTHON    Python executable (default: python3)
#   PETKLS_MODULE    Environment module to load before execution
#   PETKLS_CONDA_ENV Existing Conda environment to activate

set -euo pipefail

usage() {
    echo "Usage: sbatch $0 INPUT_CSV OUTPUT_DIRECTORY [PYTHON_OPTIONS...]" >&2
    echo "Example: sbatch $0 data/kls_edges.csv outputs --overwrite" >&2
}

if [[ $# -lt 2 ]]; then
    usage
    exit 2
fi

input_csv=$1
output_directory=$2
shift 2

# SLURM normally starts the job in the submission directory. Make this
# explicit so relative input and output paths behave consistently.
cd "${SLURM_SUBMIT_DIR:-$PWD}"

script_directory=$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1
    pwd -P
)
python_script="${script_directory}/recompute_graph_properties.py"

if [[ ! -f "$python_script" ]]; then
    echo "ERROR: Python script not found: $python_script" >&2
    exit 2
fi

if [[ -n "${PETKLS_MODULE:-}" ]]; then
    if ! command -v module >/dev/null 2>&1; then
        echo "ERROR: PETKLS_MODULE is set but the module command is unavailable." >&2
        exit 2
    fi
    module load "$PETKLS_MODULE"
fi

if [[ -n "${PETKLS_CONDA_ENV:-}" ]]; then
    if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: PETKLS_CONDA_ENV is set but conda is unavailable." >&2
        exit 2
    fi
    conda_base=$(conda info --base)
    # shellcheck source=/dev/null
    source "${conda_base}/etc/profile.d/conda.sh"
    conda activate "$PETKLS_CONDA_ENV"
fi

python_executable=${PETKLS_PYTHON:-python3}
if ! command -v "$python_executable" >/dev/null 2>&1; then
    echo "ERROR: Python executable not found: $python_executable" >&2
    exit 2
fi

echo "[JOB] ID:       ${SLURM_JOB_ID:-local}"
echo "[JOB] Host:     $(hostname)"
echo "[JOB] Started:  $(date --iso-8601=seconds)"
echo "[JOB] Input:    $input_csv"
echo "[JOB] Output:   $output_directory"
echo "[JOB] Python:   $python_executable"

"$python_executable" "$python_script" \
    --input "$input_csv" \
    --output-dir "$output_directory" \
    "$@"

echo "[JOB] Finished: $(date --iso-8601=seconds)"

