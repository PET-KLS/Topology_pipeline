# Schaefer-100 KLS graph-property recomputation

Standalone Python workflow for reconstructing subject-level Schaefer-100
Kullback–Leibler similarity (KLS) networks and recomputing weighted regional
and global graph properties. The workflow evaluates the complete positive-
weight network and proportional-threshold sensitivity networks retaining the
strongest 30%, 20%, and 10% of all 4,950 possible undirected edges.

No participant-level data are included in this repository.

## Repository contents

| File | Purpose |
| --- | --- |
| `recompute_graph_properties.py` | Standalone command-line analysis |
| `run_graph_properties.sh` | Optional SLURM submission wrapper |
| `requirements.txt` | Python dependencies |
| `CITATION.cff` | Software citation metadata |
| `LICENSE` | BSD 3-Clause software license |

## What the script does

- detects Schaefer-100 edge columns and verifies all 4,950 undirected pairs;
- reconciles duplicate edge directions after a symmetry check;
- validates finite KLS weights in the interval `[0, 1]`;
- computes weighted graph measures using KLS as similarity and `1 / KLS` as
  shortest-path distance;
- repeats global and regional calculations at full, 30%, 20%, and 10% density;
- records failures and quality-control information without dropping rows from
  the master output; and
- never overwrites the input CSV.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Input format

The input is one CSV row per scan or subject. It must contain:

- a `subject_id` column (customizable with `--subject-column`);
- an optional `session_id` column (customizable with `--session-column`); and
- all 4,950 unique Schaefer-100 ROI pairs as edge columns.

ROI labels may optionally begin with `17Networks_`. Edge columns must follow
`<ROI_A>-<ROI_B>`. If both directions are present, their values must agree to
within the configured symmetry tolerance. Non-edge columns are preserved.

Because participant-level neuroimaging tables may be sensitive, keep input and
output CSV files outside Git or under the ignored `data/` and `outputs/`
directories.

## Usage

```bash
python recompute_graph_properties.py \
  --input data/clinical_kls_edges.csv \
  --output-dir outputs
```

On Windows PowerShell:

```powershell
python .\recompute_graph_properties.py `
  --input .\data\clinical_kls_edges.csv `
  --output-dir .\outputs
```

Run `python recompute_graph_properties.py --help` for all options. In
particular, `--compute-nodal-local-efficiency` enables the computationally
expensive regional local-efficiency calculation, and `--overwrite` permits
replacement of existing output files.

## SLURM execution

Submit the same analysis to a SLURM cluster with:

```bash
sbatch run_graph_properties.sh \
  data/clinical_kls_edges.csv \
  outputs
```

Additional Python options are forwarded unchanged. For example:

```bash
sbatch run_graph_properties.sh \
  data/clinical_kls_edges.csv \
  outputs \
  --compute-nodal-local-efficiency \
  --overwrite
```

The wrapper supports three optional environment variables:

- `PETKLS_PYTHON`: Python executable, default `python3`;
- `PETKLS_MODULE`: environment module to load; and
- `PETKLS_CONDA_ENV`: existing Conda environment to activate.

For example:

```bash
PETKLS_MODULE=Python/3.11 sbatch run_graph_properties.sh data/kls.csv outputs
```

Cluster resource directives near the top of `run_graph_properties.sh` use
portable defaults (one CPU, 8 GB memory, and four hours). Adjust the time,
memory, partition, account, or quality-of-service directives for the target
cluster. No institution-specific partition or account is hardcoded.

## Outputs

For an input named `clinical_kls_edges.csv`, the output directory contains:

| File | Contents |
| --- | --- |
| `clinical_kls_edges_graphprops_recomputed.csv` | Original columns, corrected full-network regional measures, global measures for all constructions, and graph status fields |
| `graph_properties_regional_all_densities.csv` | Tidy ROI-level measures for every successful construction |
| `graph_processing_manifest_all_densities.csv` | One quality-control/status row per input row and construction |

The input is read only and is never modified.

## Methodological notes

- KLS values are treated as connection strengths/similarities.
- Weighted shortest-path measures use inverse KLS (`1 / KLS`) as edge length.
- Proportional thresholds retain an exact number of the strongest positive
  edges. Stable sorting makes tied-edge selection deterministic.
- Participation coefficient and within-module degree z-score use the fixed
  Schaefer/Yeo-17 partition.
- Louvain modularity is evaluated repeatedly with deterministic seeds; the
  maximum finite modularity is retained.
- Characteristic path length excludes infinite paths. The output therefore
  also reports `connected_pair_fraction`, which should be inspected for sparse
  constructions.
- Local nodal efficiency is disabled by default because it is substantially
  more expensive than the other measures. When disabled, its output columns
  contain missing values and the manifest records that it was not computed.

## Reproducibility

Record the Python version and installed dependency versions alongside each
analysis. For a frozen environment, run:

```bash
python -m pip freeze > environment-lock.txt
```

Do not commit participant-level data or generated tables unless their public
release has been explicitly approved and appropriately de-identified.

## License and citation

Code is released under the BSD 3-Clause License. Citation metadata are in
[`CITATION.cff`](CITATION.cff). Add the associated article DOI and repository
URL to that file once they are available.
