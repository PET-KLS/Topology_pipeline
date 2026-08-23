# PET-KLS graph properties for Schaefer-100 connectomes

This repository computes regional and global graph properties from symmetric
100 × 100 PET-KLS connectomes parcellated with the
Schaefer-100/Yeo-17 atlas. It contains one Python program and one optional
SLURM wrapper:

- `compute_graph_properties.py`: validates the matrices and computes all graph
  properties.
- `run_graph_properties.sh`: submits the Python program as a SLURM job.

KLS values are interpreted as non-negative connection strengths between zero
and one. Larger KLS therefore represents greater regional similarity. Measures
defined directly on connection strengths use KLS as supplied. Measures based
on shortest paths instead use the connection-length transformation

\[
L_{ij}=\frac{1}{W_{ij}},
\]

where \(W_{ij}>0\) is the KLS weight. Thus, a stronger KLS connection becomes a
shorter path. Zero-weight pairs are treated as absent edges.

## Requirements

- Python 3.10 or later
- NumPy
- pandas
- NetworkX
- bctpy

Install the Python dependencies with:

```bash
python -m pip install numpy pandas networkx bctpy
```

## Input data

Each input connectome must represent one subject and session and must:

- contain 100 rows and 100 columns;
- be symmetric within the requested tolerance;
- contain only finite numeric KLS values in the interval \([0,1]\); and
- use a zero diagonal.

The reader accepts a headerless matrix and safely recognizes the common pandas
numeric header and/or index. It does not silently crop matrices with unexpected
dimensions.

The labels file must contain exactly 100 ordered Schaefer-100/Yeo-17 region
labels. Standard labels such as `17Networks_LH_VisCent_ExStr_1` are used to
derive the 17 fixed atlas systems. By default, labels are read from the second
column (`--label-column 1`, using zero-based indexing).

Batch input files are discovered using this default pattern:

```text
sub-*/ses-*/*KLS_connectome.csv
```

## Graph-property definitions

Let \(W\) denote the processed KLS strength matrix, \(N=100\), \(L=1/W\) for
positive edges, and \(D\) the weighted shortest-path matrix derived from
\(L\). Thresholding, if requested, is applied before every metric is computed.

### Regional metrics

| Output column | Definition | Basis and interpretation |
| --- | --- | --- |
| `strength` | \(s_i=\sum_{j\ne i}W_{ij}\) | Sum of the KLS strengths incident on region \(i\). Higher values indicate greater overall similarity with the rest of the brain. |
| `mean_incident_kls` | \(s_i/(N-1)\) | Mean KLS between region \(i\) and all other atlas regions. For an unthresholded 100-node matrix, this is `strength / 99`. |
| `degree` | \(k_i=\sum_{j\ne i}\mathbf{1}(W_{ij}>0)\) | Number of retained nonzero edges. Without thresholding, this will usually equal 99 for every region and is therefore mainly a density/QC descriptor. |
| `clustering_coeff` | NetworkX weighted clustering coefficient | Geometric intensity of triangles around a node, computed from KLS strengths. NetworkX normalizes weights by the largest network weight. No distance conversion is applied. |
| `betweenness` | Normalized fraction of shortest paths passing through region \(i\) | Computed by NetworkX using \(L_{ij}=1/W_{ij}\) as edge distance. The default is scientifically interpretable because high KLS becomes short distance. |
| `eigenvector_centrality` | \(Wx=\lambda x\) | A region is central when it has strong connections to other central regions. Computed with the Brain Connectivity Toolbox implementation for undirected weighted networks. |
| `nodal_efficiency` | Weighted local efficiency | Efficiency of communication among the neighbors of region \(i\), computed with `bct.efficiency_wei(W, local=True)`. KLS is supplied as strength because the BCT routine performs its required path conversion internally. |
| `atlas_module` | Schaefer/Yeo-17 system label | Fixed atlas system derived from each ordered ROI label. This is not estimated from the subject's connectome. |
| `louvain_community` | Best subject-specific Louvain assignment | Data-driven community label from repeated seeded Louvain optimization on KLS strengths. It is reported separately and is not used for the atlas-based participation metrics. |
| `participation_coeff` | \(P_i=1-\sum_m(s_{im}/s_i)^2\) | Distribution of regional strength across the fixed Schaefer/Yeo-17 systems. Higher values indicate more evenly distributed between-system connectivity. |
| `within_module_degree_z` | Within-system strength z-score | Standardized strength of region \(i\) relative to other regions in the same fixed Schaefer/Yeo-17 system. |

### Global metrics

| Output column | Definition | Basis and interpretation |
| --- | --- | --- |
| `global_efficiency` | Mean \(1/D_{ij}\) over ordered off-diagonal pairs | Weighted integration computed with `bct.efficiency_wei(W)`. The function receives KLS strengths and performs the required internal path conversion. Disconnected pairs contribute zero efficiency. |
| `char_path_len` | Mean finite \(D_{ij}\), excluding the diagonal | Characteristic weighted path length after converting KLS strengths to \(1/W\). Infinite paths are excluded and network connectedness is reported separately. |
| `avg_clustering_coef` | Mean regional weighted clustering | Whole-network segregation summary based directly on KLS strengths. |
| `density` | Retained edges divided by 4,950 possible undirected edges | Proportion of nonzero off-diagonal connections after optional thresholding. |
| `n_vertices` | Number of graph nodes | Expected to equal 100. |
| `n_edges` | Number of retained undirected edges | Maximum 4,950 for 100 nodes. |
| `avg_strength` | Mean regional strength | Average sum of incident KLS weights. For an undirected network, this equals twice the total edge weight divided by 100. |
| `avg_mean_incident_kls` | Mean regional `mean_incident_kls` | Whole-connectome mean KLS expressed per possible incident edge. |
| `avg_betweenness` | Mean normalized regional betweenness | Whole-network average of inverse-KLS shortest-path betweenness. |
| `modularity_Q` | Highest Louvain modularity across seeded runs | Strength of the best data-driven community division found across `--modularity-runs` attempts. Default: 10 runs with seeds beginning at 42. |
| `n_louvain_communities` | Number of communities in the retained Louvain solution | Describes the selected subject-specific modular partition. |
| `connected_pair_fraction` | Finite ordered off-diagonal paths divided by \(N(N-1)\) | Equals one for a connected graph and identifies disconnection after thresholding. |
| `betweenness_mode` | Path convention recorded in the output | `inverse` is the corrected default. `notebook-legacy` exists only to reproduce the historical calculation that incorrectly treated larger KLS as longer distance. |
| `participation_partition` | Partition used for participation metrics | Always recorded as `Schaefer_Yeo17` in this implementation. |

The `shortest_paths.csv` output contains the complete \(D\) matrix computed
from inverse-KLS connection lengths.

## Thresholding

The default is `--threshold-mode none`, which retains every positive KLS edge.
This preserves weighted information and matches the unthresholded analysis.
Two optional modes are available for sensitivity analyses:

- `absolute`: removes weights below the supplied KLS value;
- `quantile`: removes the stated lower fraction of positive weights within each
  connectome.

For example, `--threshold-mode quantile --threshold 0.3` removes the lowest 30%
of positive weights. It does **not** apply an absolute KLS cutoff of 0.3.

Threshold choice changes graph density and can affect every downstream metric.
Analyses should state the mode and value and, where relevant, evaluate the
stability of findings across proportional thresholds. Zero should only denote
an absent edge; missing or invalid measurements should not be encoded as zero.

## Usage

### One connectome

```bash
python compute_graph_properties.py \
  --input-csv INPUT_CONNECTOME.csv \
  --output-root OUTPUT_DIRECTORY \
  --labels-file SCHAEFER100_YEO17_LABELS.csv
```

### Batch processing

```bash
python compute_graph_properties.py \
  --input-root INPUT_DIRECTORY \
  --output-root OUTPUT_DIRECTORY \
  --labels-file SCHAEFER100_YEO17_LABELS.csv
```

### SLURM

```bash
sbatch run_graph_properties.sh \
  INPUT_DIRECTORY \
  OUTPUT_DIRECTORY \
  SCHAEFER100_YEO17_LABELS.csv
```

Extra Python arguments may follow the labels file:

```bash
sbatch run_graph_properties.sh \
  INPUT_DIRECTORY \
  OUTPUT_DIRECTORY \
  SCHAEFER100_YEO17_LABELS.csv \
  --threshold-mode quantile \
  --threshold 0.3
```

The wrapper supports these optional environment variables:

- `PETKLS_PYTHON`: Python executable, default `python3`;
- `PETKLS_MODULE`: environment module to load; and
- `PETKLS_CONDA_ENV`: existing conda environment to activate.

## Outputs

For each `sub-*`/`ses-*` pair, the program writes:

```text
OUTPUT_DIRECTORY/
├── processing_manifest.csv
└── sub-<ID>/
    └── ses-<ID>/
        ├── sub-<ID>_ses-<ID>_Schaefer2018_100Parcels_17Networks_global_metrics.csv
        ├── sub-<ID>_ses-<ID>_Schaefer2018_100Parcels_17Networks_regional_metrics.csv
        └── sub-<ID>_ses-<ID>_Schaefer2018_100Parcels_17Networks_shortest_paths.csv
```

Outputs are written atomically. Complete subject/session outputs are skipped
unless `--overwrite` is supplied. The manifest records completed, skipped, and
failed inputs without stopping the remaining batch.

The global output also records matrix QC, processing parameters, thresholding,
random seeds, and the metric partitions used, making each result auditable.

## Methodological notes

- Similarity and distance are not interchangeable. NetworkX weighted
  betweenness treats its edge-weight attribute as a distance, so inverse KLS is
  required for the corrected analysis.
- `--betweenness-mode notebook-legacy` is provided only for comparison with
  historical notebook results and should not be used for scientific
  interpretation.
- Fixed Schaefer/Yeo-17 participation and subject-specific Louvain modularity
  answer different questions. The former quantifies integration across known
  functional systems; the latter describes a data-driven community structure.
- Louvain solutions can vary because modularity optimization is degenerate.
  The program therefore performs repeated deterministic seeded runs and retains
  the solution with the highest modularity value.
- Characteristic path length excludes infinite paths. Always inspect
  `connected_pair_fraction`, particularly after thresholding.

## Software references

- [NetworkX weighted clustering](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.cluster.clustering.html)
- [NetworkX betweenness centrality](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.betweenness_centrality.html)
- [Brain Connectivity Toolbox measures](https://sites.google.com/site/bctnet/list-of-measures)
- [bctpy](https://github.com/aestrivex/bctpy)

