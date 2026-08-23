#!/usr/bin/env python3
"""Compute graph properties for Schaefer-100 PET-KLS connectomes.

Required packages:
    numpy, pandas, bctpy

The program accepts either one connectome (``--input-csv``) or a directory
containing BIDS-like ``sub-*/ses-*`` folders (``--input-root``). KLS values
are interpreted as non-negative connection strengths. For path-based
measures, non-zero strengths are converted to lengths as ``1 / KLS``.

Examples
--------
Single connectome::

    python compute_graph_properties.py \
      --input-csv sub-001_ses-001_KLS_connectome.csv \
      --output-root derivatives/graph_properties

Batch::

    python compute_graph_properties.py \
      --input-root derivatives/pet_connectomes_kls \
      --output-root derivatives/graph_properties \
      --labels-file Schaefer2018_100Parcels_17Networks_labels.csv
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

try:
    import bct
except ImportError as exc:
    raise SystemExit(
        "Missing dependency 'bct' (provided by bctpy). Install numpy, "
        "pandas, and bctpy in the execution environment."
    ) from exc


N_ROIS = 100
DEFAULT_PATTERN = "sub-*/ses-*/*KLS_connectome.csv"
ATLAS_NAME = "Schaefer2018_100Parcels_17Networks"


@dataclass(frozen=True)
class MatrixQC:
    raw_rows: int
    raw_columns: int
    removed_header_row: bool
    removed_index_column: bool
    max_asymmetry_before_fix: float
    zero_weight_pairs: int


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute Schaefer-100 graph properties from PET-KLS matrices."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input-csv", type=Path, help="One KLS connectome CSV.")
    source.add_argument(
        "--input-root",
        type=Path,
        help="Root searched recursively using --pattern.",
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--pattern", default=DEFAULT_PATTERN)
    parser.add_argument(
        "--labels-file",
        type=Path,
        help="Optional CSV/TXT containing 100 ROI labels.",
    )
    parser.add_argument(
        "--label-column",
        type=int,
        default=1,
        help="Zero-based label column (default: 1).",
    )
    parser.add_argument(
        "--threshold-mode",
        choices=("none", "absolute", "quantile"),
        default="none",
        help=(
            "none keeps all positive edges; absolute removes weights below "
            "--threshold; quantile removes the stated lower fraction."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Threshold value used by absolute or quantile mode.",
    )
    parser.add_argument(
        "--modularity-runs",
        type=int,
        default=10,
        help="Seeded Louvain runs; the highest-Q solution is retained.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--symmetry-tolerance",
        type=float,
        default=1e-8,
        help="Maximum accepted |W-W.T| before averaging tiny differences.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace complete outputs instead of skipping them.",
    )
    return parser


def _is_counter(values: np.ndarray, n: int) -> bool:
    """Return True only for an exact 0..n-1 or 1..n numeric counter."""
    if values.shape != (n,) or not np.isfinite(values).all():
        return False
    return bool(
        np.allclose(values, np.arange(n), rtol=0, atol=1e-12)
        or np.allclose(values, np.arange(1, n + 1), rtol=0, atol=1e-12)
    )


def read_kls_matrix(path: Path, symmetry_tolerance: float) -> tuple[np.ndarray, MatrixQC]:
    """Read a headerless matrix or safely remove a recognized CSV header/index."""
    if not path.is_file():
        raise FileNotFoundError(path)

    raw_frame = pd.read_csv(path, header=None, dtype=str)
    raw_shape = raw_frame.shape
    numeric = raw_frame.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    removed_header = False
    removed_index = False

    # PET-KLS matrices written by pandas with index=False have a 0..99 header
    # and therefore appear as 101 x 100 when deliberately read with header=None.
    if numeric.shape == (N_ROIS + 1, N_ROIS) and _is_counter(numeric[0], N_ROIS):
        numeric = numeric[1:, :]
        removed_header = True
    elif numeric.shape == (N_ROIS, N_ROIS + 1) and _is_counter(numeric[:, 0], N_ROIS):
        numeric = numeric[:, 1:]
        removed_index = True
    elif numeric.shape == (N_ROIS + 1, N_ROIS + 1):
        if _is_counter(numeric[0, 1:], N_ROIS) and _is_counter(
            numeric[1:, 0], N_ROIS
        ):
            numeric = numeric[1:, 1:]
            removed_header = True
            removed_index = True

    if numeric.shape != (N_ROIS, N_ROIS):
        raise ValueError(
            f"expected a {N_ROIS}x{N_ROIS} matrix (optionally with a recognized "
            f"pandas header/index), found {raw_shape}"
        )
    if not np.isfinite(numeric).all():
        locations = np.argwhere(~np.isfinite(numeric))[:5].tolist()
        raise ValueError(f"matrix contains non-finite/non-numeric values at {locations}")
    if np.any(numeric < -1e-12):
        raise ValueError("KLS matrix contains negative weights")
    if np.any(numeric > 1 + 1e-8):
        raise ValueError("KLS matrix contains weights above 1")

    numeric = np.clip(numeric, 0.0, 1.0)
    asymmetry = float(np.max(np.abs(numeric - numeric.T)))
    if asymmetry > symmetry_tolerance:
        raise ValueError(
            f"matrix is not symmetric (maximum |W-W.T|={asymmetry:.6g}; "
            f"tolerance={symmetry_tolerance:.6g})"
        )
    numeric = (numeric + numeric.T) / 2.0
    np.fill_diagonal(numeric, 0.0)
    upper = numeric[np.triu_indices(N_ROIS, 1)]

    qc = MatrixQC(
        raw_rows=raw_shape[0],
        raw_columns=raw_shape[1],
        removed_header_row=removed_header,
        removed_index_column=removed_index,
        max_asymmetry_before_fix=asymmetry,
        zero_weight_pairs=int(np.count_nonzero(upper == 0)),
    )
    return numeric, qc


def load_labels(path: Path | None, column: int) -> list[str]:
    if path is None:
        return [f"ROI_{i:03d}" for i in range(1, N_ROIS + 1)]
    if not path.is_file():
        raise FileNotFoundError(path)

    # sep=None allows comma-, semicolon-, or tab-delimited label files.
    frame = pd.read_csv(path, sep=None, engine="python", header=None, dtype=str)
    if column < 0 or column >= frame.shape[1]:
        raise ValueError(
            f"label column {column} is unavailable; file has {frame.shape[1]} columns"
        )
    labels = frame.iloc[:, column].fillna("").str.strip().tolist()
    if len(labels) == N_ROIS + 1:
        labels = labels[1:]
    if len(labels) != N_ROIS or any(not value for value in labels):
        raise ValueError(
            f"labels file must yield exactly {N_ROIS} non-empty labels; found {len(labels)}"
        )
    if len(set(labels)) != N_ROIS:
        raise ValueError("ROI labels must be unique")
    return labels


def apply_threshold(matrix: np.ndarray, mode: str, value: float) -> np.ndarray:
    weights = matrix.copy()
    positive = weights[np.triu_indices_from(weights, 1)]
    positive = positive[positive > 0]

    if mode == "none":
        if value != 0:
            raise ValueError("--threshold must be 0 when --threshold-mode=none")
        return weights
    if positive.size == 0:
        raise ValueError("matrix has no positive edges")
    if mode == "absolute":
        if not 0 <= value <= 1:
            raise ValueError("absolute threshold must be between 0 and 1")
        cutoff = value
    else:
        if not 0 <= value < 1:
            raise ValueError("quantile threshold must be in [0, 1)")
        cutoff = float(np.quantile(positive, value))

    weights[weights < cutoff] = 0.0
    np.fill_diagonal(weights, 0.0)
    return weights


def strengths_to_lengths(weights: np.ndarray) -> np.ndarray:
    lengths = np.zeros_like(weights)
    mask = weights > 0
    lengths[mask] = 1.0 / weights[mask]
    return lengths


def best_louvain(
    weights: np.ndarray, runs: int, seed: int
) -> tuple[np.ndarray, float]:
    if runs < 1:
        raise ValueError("--modularity-runs must be at least 1")
    best_community: np.ndarray | None = None
    best_q = -np.inf
    failures: list[str] = []
    for offset in range(runs):
        try:
            community, q_value = bct.modularity_louvain_und(
                weights, seed=seed + offset
            )
            q_float = float(q_value)
            if np.isfinite(q_float) and q_float > best_q:
                best_community = np.asarray(community, dtype=int)
                best_q = q_float
        except Exception as exc:  # retain other successful seeded runs
            failures.append(str(exc))
    if best_community is None:
        detail = failures[0] if failures else "no finite modularity result"
        raise RuntimeError(f"all Louvain runs failed: {detail}")
    return best_community, best_q


def compute_metrics(
    weights: np.ndarray, labels: Sequence[str], modularity_runs: int, seed: int
) -> tuple[dict[str, float | int], pd.DataFrame, pd.DataFrame]:
    lengths = strengths_to_lengths(weights)
    distances, _ = bct.distance_wei(lengths)
    characteristic_path, _, _, _, _ = bct.charpath(
        distances, include_diagonal=False, include_infinite=False
    )
    community, modularity_q = best_louvain(weights, modularity_runs, seed)

    strength = np.sum(weights, axis=0)
    degree = np.count_nonzero(weights, axis=0)
    clustering = np.asarray(bct.clustering_coef_wu(weights), dtype=float)
    betweenness = np.asarray(bct.betweenness_wei(lengths), dtype=float)
    eigenvector = np.asarray(bct.eigenvector_centrality_und(weights), dtype=float)
    nodal_efficiency = np.asarray(bct.efficiency_wei(weights, local=True), dtype=float)
    participation = np.asarray(bct.participation_coef(weights, community), dtype=float)
    within_module_z = np.asarray(
        bct.module_degree_zscore(weights, community), dtype=float
    )
    density, n_vertices, n_edges = bct.density_und(weights)

    finite_off_diagonal = np.isfinite(distances) & ~np.eye(N_ROIS, dtype=bool)
    connected_fraction = float(np.count_nonzero(finite_off_diagonal) / (N_ROIS * (N_ROIS - 1)))
    global_metrics: dict[str, float | int] = {
        "global_efficiency": float(bct.efficiency_wei(weights)),
        "characteristic_path_length": float(characteristic_path),
        "mean_clustering_coefficient": float(np.mean(clustering)),
        "density": float(density),
        "n_vertices": int(n_vertices),
        "n_edges": int(n_edges),
        "mean_strength": float(np.mean(strength)),
        "mean_betweenness": float(np.mean(betweenness)),
        "modularity_q": modularity_q,
        "n_communities": int(np.unique(community).size),
        "connected_pair_fraction": connected_fraction,
    }
    regional = pd.DataFrame(
        {
            "roi": labels,
            "strength": strength,
            "degree": degree,
            "betweenness": betweenness,
            "clustering": clustering,
            "eigenvector_centrality": eigenvector,
            "nodal_efficiency": nodal_efficiency,
            "community": community,
            "participation_coefficient": participation,
            "within_module_degree_z": within_module_z,
        }
    )
    shortest_paths = pd.DataFrame(distances, index=labels, columns=labels)
    shortest_paths.index.name = "roi"
    return global_metrics, regional, shortest_paths


def subject_session(path: Path) -> tuple[str, str]:
    subject = next((part for part in path.parts if part.startswith("sub-")), None)
    session = next((part for part in path.parts if part.startswith("ses-")), None)
    if subject and session:
        return subject, session
    match = re.search(r"(sub-[^_]+)_(ses-[^_]+)", path.name)
    if not match:
        raise ValueError("could not derive sub-* and ses-* identifiers from path")
    return match.group(1), match.group(2)


def atomic_to_csv(frame: pd.DataFrame, destination: Path, *, index: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    os.close(handle)
    temporary = Path(temporary_name)
    try:
        frame.to_csv(temporary, index=index)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def output_paths(output_root: Path, subject: str, session: str) -> dict[str, Path]:
    folder = output_root / subject / session
    prefix = f"{subject}_{session}_{ATLAS_NAME}"
    return {
        "global": folder / f"{prefix}_global_metrics.csv",
        "regional": folder / f"{prefix}_regional_metrics.csv",
        "paths": folder / f"{prefix}_shortest_paths.csv",
    }


def discover_inputs(args: argparse.Namespace) -> list[Path]:
    if args.input_csv:
        return [args.input_csv.resolve()]
    if not args.input_root.is_dir():
        raise NotADirectoryError(args.input_root)
    paths = sorted(path.resolve() for path in args.input_root.glob(args.pattern))
    if not paths:
        raise FileNotFoundError(
            f"no files matched {args.pattern!r} beneath {args.input_root}"
        )
    return paths


def reject_output_collisions(paths: Sequence[Path]) -> None:
    seen: dict[tuple[str, str], Path] = {}
    for path in paths:
        key = subject_session(path)
        if key in seen:
            raise ValueError(
                f"multiple connectomes resolve to {key[0]} {key[1]}: "
                f"{seen[key]} and {path}"
            )
        seen[key] = path


def process_one(
    path: Path, args: argparse.Namespace, labels: Sequence[str]
) -> dict[str, object]:
    subject, session = subject_session(path)
    outputs = output_paths(args.output_root, subject, session)
    if not args.overwrite and all(output.is_file() for output in outputs.values()):
        return {
            "input_csv": str(path),
            "subject": subject,
            "session": session,
            "status": "skipped_complete",
            "message": "all outputs already exist",
        }

    matrix, qc = read_kls_matrix(path, args.symmetry_tolerance)
    matrix = apply_threshold(matrix, args.threshold_mode, args.threshold)
    global_metrics, regional, shortest_paths = compute_metrics(
        matrix, labels, args.modularity_runs, args.seed
    )
    global_row: dict[str, object] = {
        "subject": subject,
        "session": session,
        "atlas": ATLAS_NAME,
        "threshold_mode": args.threshold_mode,
        "threshold": args.threshold,
        "modularity_runs": args.modularity_runs,
        "seed": args.seed,
        **global_metrics,
        **qc.__dict__,
    }
    regional.insert(0, "session", session)
    regional.insert(0, "subject", subject)
    atomic_to_csv(pd.DataFrame([global_row]), outputs["global"], index=False)
    atomic_to_csv(regional, outputs["regional"], index=False)
    atomic_to_csv(shortest_paths, outputs["paths"], index=True)
    return {
        "input_csv": str(path),
        "subject": subject,
        "session": session,
        "status": "completed",
        "message": "",
    }


def write_manifest(rows: Iterable[dict[str, object]], output_root: Path) -> None:
    atomic_to_csv(
        pd.DataFrame(rows), output_root / "processing_manifest.csv", index=False
    )


def main() -> int:
    args = build_parser().parse_args()
    args.output_root = args.output_root.resolve()
    try:
        labels = load_labels(args.labels_file, args.label_column)
        inputs = discover_inputs(args)
        reject_output_collisions(inputs)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    records: list[dict[str, object]] = []
    failures = 0
    for index, path in enumerate(inputs, start=1):
        print(f"[{index}/{len(inputs)}] {path}")
        try:
            record = process_one(path, args, labels)
            print(f"  {record['status']}")
        except Exception as exc:
            failures += 1
            try:
                subject, session = subject_session(path)
            except Exception:
                subject, session = "", ""
            record = {
                "input_csv": str(path),
                "subject": subject,
                "session": session,
                "status": "failed",
                "message": f"{type(exc).__name__}: {exc}",
            }
            print(f"  FAILED: {record['message']}", file=sys.stderr)
        records.append(record)
        write_manifest(records, args.output_root)

    completed = sum(record["status"] == "completed" for record in records)
    skipped = sum(record["status"] == "skipped_complete" for record in records)
    print(f"Finished: {completed} completed, {skipped} skipped, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
