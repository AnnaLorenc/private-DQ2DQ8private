# diversity/changed_diversity_AL.py

"""
Diversity analysis functions for AIRR-seq repertoires.

Implements diversity_indices, diversity metrics, and repertoire comparison functions
similar to the R LymphoSeq2 package - polarised and extended to handle both Polars and pandas DataFrames.
"""

import numpy as np
import polars as pl
import pandas as pd
from typing import Union, Optional, List


# ---------------------------------------------------------------------------
# Statistical helper functions (inlined from lymphoseq.analysis.statistics)
# ---------------------------------------------------------------------------


def gini_coefficient(values: Union[List[float], np.ndarray]) -> float:
    """Calculate the Gini coefficient of inequality (0 = equal, 1 = maximal)."""
    if len(values) == 0:
        return 0.0
    values = np.array(values, dtype=float)
    values = values[values > 0]
    if len(values) <= 1:
        return 0.0
    sorted_values = np.sort(values)
    n = len(sorted_values)
    index = np.arange(1, n + 1)
    gini = (2 * np.sum(index * sorted_values)) / (n * np.sum(sorted_values)) - (n + 1) / n
    return max(0.0, min(1.0, gini))


def shannon_entropy(frequencies: Union[List[float], np.ndarray], base: float = 2.0) -> float:
    """Calculate Shannon entropy. Higher values indicate more diversity."""
    if len(frequencies) == 0:
        return 0.0
    frequencies = np.array(frequencies, dtype=float)
    frequencies = frequencies[frequencies > 0]
    if len(frequencies) <= 1:
        return 0.0
    if not np.isclose(np.sum(frequencies), 1.0):
        frequencies = frequencies / np.sum(frequencies)
    entropy = -np.sum(frequencies * np.log(frequencies) / np.log(base))
    return max(0.0, entropy)



def inverse_simpson(counts: Union[np.ndarray, pl.Series, List]) -> float:
    """
    Calculate the Inverse Simpson diversity index.

    Accepts numpy arrays, Polars Series, or plain lists.
    Defined as 1 / sum(p_i^2). Ranges from 1 (single dominant clone)
    to N (all clones equally abundant).
    """
    if isinstance(counts, pl.Series):
        counts = counts.to_numpy()
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return 0.0
    frequencies = counts / total
    simpson = np.sum(frequencies ** 2)
    return 1.0 / simpson if simpson > 0 else 0.0


def _normalized_clonality(counts: np.ndarray) -> float:
    """1 - normalized Shannon entropy from a counts array."""
    total = counts.sum()
    if total == 0:
        return 0.0
    frequencies = counts / total
    n_unique = int(np.sum(counts > 0))
    max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
    norm_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0.0
    return 1.0 - norm_entropy


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def diversity_indices(
    data: Union[pl.DataFrame, pd.DataFrame],
    rarefy: bool = False,
    min_count: Optional[int] = None,
    iterations: int = 100,
    group_by: str = "repertoire_id",
    sequence_cols: Optional[Union[str, List[str]]] = None,
    count_col: str = "duplicate_count",
) -> Union[pl.DataFrame, pd.DataFrame]:
    """
    Calculate diversity indices and diversity metrics for each repertoire.

    Args:
        data: Input data frame with AIRR-formatted sequences
        rarefy: Whether to perform rarefaction analysis
        min_count: Minimum sequence count for rarefaction
        iterations: Number of rarefaction iterations
        group_by: Column to group by (default: "repertoire_id")
        sequence_cols: Column(s) used as sequence index. If multiple columns
            are given the data is first aggregated (sum of count_col) per
            unique combination within each repertoire. Defaults to ["junction_aa"].
        count_col: Name of the abundance/count column (default: "duplicate_count")

    Returns:
        Data frame with diversity metrics for each repertoire:
        - diversity_indices: 1 - normalized Shannon entropy (0-1)
        - gini_coefficient: Measure of inequality (0-1)
        - unique_sequences: Number of unique sequences
        - total_count: Total sequence reads
        - top_sequence: Frequency of most abundant clone
        - convergence: Ratio of total sequences to unique sequences
        - inverse_simpson: Inverse Simpson diversity index
    """
    if sequence_cols is None:
        sequence_cols = ["junction_aa"]
    elif isinstance(sequence_cols, str):
        sequence_cols = [sequence_cols]

    sequence_index_name = "_".join(sequence_cols)
    rarefaction_its = iterations if rarefy else None

    is_polars = isinstance(data, pl.DataFrame)
    if is_polars:
        result = _diversity_indices_polars(data, rarefy, min_count, iterations, group_by, sequence_cols, count_col)
        return result.with_columns([
            pl.lit(sequence_index_name).alias("sequence_index"),
            pl.lit(rarefaction_its).alias("rarefaction_iterations"),
        ])
    else:
        result = _diversity_indices_pandas(data, rarefy, min_count, iterations, group_by, sequence_cols, count_col)
        result["sequence_index"] = sequence_index_name
        result["rarefaction_iterations"] = rarefaction_its
        return result


def _aggregate_polars(
    data: pl.DataFrame,
    group_by: str,
    sequence_cols: List[str],
    count_col: str,
) -> pl.DataFrame:
    """Aggregate count_col by (group_by + sequence_cols), summing counts."""
    print(f"[Polars] Input rows: {len(data)}")
    result = (
        data
        .group_by([group_by] + sequence_cols)
        .agg(pl.col(count_col).sum())
    )
    print(f"[Polars] Rows after aggregation by {[group_by] + sequence_cols}: {len(result)}")
    return result


def _aggregate_pandas(
    data: pd.DataFrame,
    group_by: str,
    sequence_cols: List[str],
    count_col: str,
) -> pd.DataFrame:
    """Aggregate count_col by (group_by + sequence_cols), summing counts."""
    print(f"[Pandas] Input rows: {len(data)}")
    result = (
        data
        .groupby([group_by] + sequence_cols, as_index=False)[count_col]
        .sum()
    )
    print(f"[Pandas] Rows after aggregation by {[group_by] + sequence_cols}: {len(result)}")
    return result


def _diversity_indices_polars(
    data: pl.DataFrame,
    rarefy: bool,
    min_count: Optional[int],
    iterations: int,
    group_by: str,
    sequence_cols: List[str],
    count_col: str,
) -> pl.DataFrame:
    """Calculate diversity indices using native Polars aggregations."""

    required_cols = [group_by, count_col] + sequence_cols
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    data = data.with_columns(
        pl.col(count_col).cast(pl.Float64, strict=False).fill_null(1.0)
    )

    # Pre-aggregate
    data = _aggregate_polars(data, group_by, sequence_cols, count_col)

    if not rarefy:
        results = (
            data
            .group_by(group_by)
            .agg([
                pl.len().alias("total_sequences"),
                pl.struct(sequence_cols).n_unique().alias("unique_sequences"),
                pl.col(count_col).sum().alias("total_count"),
                pl.col(count_col).map_elements(
                    lambda s: _normalized_clonality(s.to_numpy()),
                    return_dtype=pl.Float64
                ).first().alias("clonality"),
                pl.col(count_col).map_elements(
                    lambda s: gini_coefficient(s.to_numpy()),
                    return_dtype=pl.Float64
                ).first().alias("gini_coefficient"),
                (pl.col(count_col).max() / pl.col(count_col).sum())
                    .first().alias("top_sequence"),
                pl.col(count_col).map_elements(
                    lambda s: inverse_simpson(s),
                    return_dtype=pl.Float64
                ).first().alias("inverse_simpson"),
                pl.col(count_col).map_elements(
                    lambda s: shannon_entropy(s.to_numpy() / s.to_numpy().sum()),
                    return_dtype=pl.Float64
                ).first().alias("shannon_index"),
            ])
            .with_columns(
                (pl.col("total_sequences") / pl.col("unique_sequences")).alias("convergence")
            )
        )
        return results
    else:
        if min_count is None:
            min_count = (
                data
                .group_by(group_by)
                .agg(pl.col(count_col).sum().alias("total"))
                .select(pl.col("total").min())
                .item()
            )
        return _perform_rarefaction_polars(data, min_count, iterations, group_by, sequence_cols, count_col)


def _diversity_indices_pandas(
    data: pd.DataFrame,
    rarefy: bool,
    min_count: Optional[int],
    iterations: int,
    group_by: str,
    sequence_cols: List[str],
    count_col: str,
) -> pd.DataFrame:
    """Calculate diversity indices using pandas."""

    required_cols = [group_by, count_col] + sequence_cols
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    data = data.copy()
    data[count_col] = pd.to_numeric(data[count_col], errors="coerce").fillna(1.0)

    # Pre-aggregate
    data = _aggregate_pandas(data, group_by, sequence_cols, count_col)

    if not rarefy:
        def calculate_metrics(group):
            counts = group[count_col].values
            frequencies = counts / counts.sum()
            n_unique = group[sequence_cols].drop_duplicates().shape[0]
            max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
            normalized_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0
            return pd.Series({
                "total_sequences": len(group),
                "unique_sequences": n_unique,
                "total_count": counts.sum(),
                "clonality": 1 - normalized_entropy,
                "gini_coefficient": gini_coefficient(counts),
                "top_sequence": counts.max() / counts.sum(),
                "convergence": len(group) / n_unique,
                "inverse_simpson": inverse_simpson(counts),
                "shannon_index": shannon_entropy(frequencies),
            })

        results = data.groupby(group_by).apply(calculate_metrics).reset_index()
    else:
        if min_count is None:
            min_count = data.groupby(group_by)[count_col].sum().min()
        results = _perform_rarefaction_pandas(data, min_count, iterations, group_by, sequence_cols, count_col)

    return results


# ---------------------------------------------------------------------------
# Rarefaction helpers
# ---------------------------------------------------------------------------

def _weighted_sample(sequences: list, counts: list, sample_size: int) -> np.ndarray:
    """Sample sequences with replacement, weighted by their counts."""
    total = sum(counts)
    probabilities = [c / total for c in counts]
    sampled_indices = np.random.choice(
        len(sequences), size=sample_size, replace=True, p=probabilities
    )
    return np.array([sequences[i] for i in sampled_indices])


def _perform_rarefaction_polars(
    data: pl.DataFrame,
    min_count: int,
    iterations: int,
    group_by: str,
    sequence_cols: List[str],
    count_col: str,
) -> pl.DataFrame:
    """Perform rarefaction analysis using Polars."""
    all_rows = []
    for repertoire in data[group_by].unique().to_list():
        repertoire_data = data.filter(pl.col(group_by) == repertoire)
        total_reads = repertoire_data[count_col].sum()
        if total_reads < min_count:
            continue

        if len(sequence_cols) == 1:
            sequences = repertoire_data[sequence_cols[0]].to_list()
        else:
            sequences = [
                "|".join(str(v) for v in row)
                for row in repertoire_data.select(sequence_cols).iter_rows()
            ]
        counts = repertoire_data[count_col].to_list()
        iteration_metrics = []

        for _ in range(iterations):
            sampled_seqs = _weighted_sample(sequences, counts, min_count)
            unique_seqs, seq_counts = np.unique(sampled_seqs, return_counts=True)
            n_unique = len(unique_seqs)
            max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
            frequencies = seq_counts / seq_counts.sum()
            normalized_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0
            iteration_metrics.append({
                "total_sequences": float(len(sampled_seqs)),
                "unique_sequences": float(n_unique),
                "total_count": float(min_count),
                "clonality": 1.0 - normalized_entropy,
                "gini_coefficient": gini_coefficient(seq_counts),
                "top_sequence": float(seq_counts.max() / seq_counts.sum()),
                "convergence": float(len(sampled_seqs)) / n_unique,
                "inverse_simpson": inverse_simpson(seq_counts),
                "shannon_index": shannon_entropy(frequencies),
            })

        # Average across iterations using Polars
        iter_df = pl.DataFrame(iteration_metrics)
        avg = iter_df.mean().with_columns(pl.lit(repertoire).alias(group_by))
        all_rows.append(avg)

    if not all_rows:
        return pl.DataFrame()
    return pl.concat(all_rows)


def _perform_rarefaction_pandas(
    data: pd.DataFrame,
    min_count: int,
    iterations: int,
    group_by: str,
    sequence_cols: List[str],
    count_col: str,
) -> pd.DataFrame:
    """Perform rarefaction analysis using pandas."""
    results = []

    for repertoire, group in data.groupby(group_by):
        total_reads = group[count_col].sum()
        if total_reads < min_count:
            continue

        if len(sequence_cols) == 1:
            sequences = group[sequence_cols[0]].values
        else:
            sequences = group[sequence_cols].apply(
                lambda r: "|".join(str(v) for v in r), axis=1
            ).values
        counts = group[count_col].values
        probabilities = counts / counts.sum()
        iteration_results = []

        for _ in range(iterations):
            sampled_seqs = np.random.choice(sequences, size=int(min_count), replace=True, p=probabilities)
            unique_seqs, seq_counts = np.unique(sampled_seqs, return_counts=True)
            frequencies = seq_counts / seq_counts.sum()
            n_unique = len(unique_seqs)
            max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
            normalized_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0
            iteration_results.append({
                "total_sequences": len(sampled_seqs),
                "unique_sequences": n_unique,
                "total_count": min_count,
                "clonality": 1 - normalized_entropy,
                "gini_coefficient": gini_coefficient(seq_counts),
                "top_sequence": seq_counts.max() / seq_counts.sum(),
                "convergence": len(sampled_seqs) / n_unique,
                "inverse_simpson": inverse_simpson(seq_counts),
                "shannon_index": shannon_entropy(frequencies),
            })

        avg_metrics = pd.DataFrame(iteration_results).mean().to_dict()
        avg_metrics[group_by] = repertoire
        results.append(avg_metrics)

    return pd.DataFrame(results)


def diversity_metrics(
    data: Union[pl.DataFrame, pd.DataFrame],
    group_by: str = "repertoire_id",
    sequence_cols: Optional[Union[str, List[str]]] = None,
    count_col: str = "duplicate_count",
) -> Union[pl.DataFrame, pd.DataFrame]:
    """Calculate comprehensive diversity metrics for repertoires."""
    return diversity_indices(data, rarefy=False, group_by=group_by, sequence_cols=sequence_cols, count_col=count_col)


