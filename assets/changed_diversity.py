"""
Diversity analysis functions for AIRR-seq repertoires.

Implements clonality, diversity metrics, and repertoire comparison functions
similar to the R LymphoSeq2 package.
"""

import numpy as np
import polars as pl
import pandas as pd
from typing import Union, Optional, List, Dict, Any
import lymphoseq as ls
from lymphoseq.analysis.statistics import gini_coefficient, shannon_entropy


def clonality(
    data: Union[pl.DataFrame, pd.DataFrame],
    rarefy: bool = False,
    min_count: Optional[int] = None,
    iterations: int = 100,
    group_by: str = "repertoire_id"
) -> Union[pl.DataFrame, pd.DataFrame]:
    """
    Calculate clonality and diversity metrics for each repertoire.

    This function calculates standard repertoire diversity metrics including
    clonality, Gini coefficient, convergence, and unique productive sequences.

    Clonality is calculated as 1 - (normalized Shannon entropy), where Shannon
    entropy is normalized by dividing by log2(n_unique_sequences). This ensures
    clonality is in the range [0, 1], where:
    - 0 indicates maximum diversity (all clones equally abundant)
    - 1 indicates minimum diversity (single dominant clone)

    Args:
        data: Input data frame with AIRR-formatted sequences
        rarefy: Whether to perform rarefaction analysis
        min_count: Minimum sequence count for rarefaction
        iterations: Number of rarefaction iterations
        group_by: Column to group by (default: repertoire_id)

    Returns:
        Data frame with diversity metrics for each repertoire:
        - clonality: 1 - normalized Shannon entropy (0-1)
        - gini_coefficient: Measure of inequality (0-1)
        - unique_sequences: Number of unique sequences
        - total_count: Total sequence reads
        - top_sequence: Frequency of most abundant clone
        - convergence: Ratio of total sequences to unique sequences

    Examples:
        >>> import lymphoseq as ls
        >>> data = ls.read_immunoseq("data/")
        >>> diversity = ls.clonality(data)
        >>> print(diversity)
    """
    is_polars = isinstance(data, pl.DataFrame)

    if is_polars:
        return _clonality_polars(data, rarefy, min_count, iterations, group_by)
    else:
        return _clonality_pandas(data, rarefy, min_count, iterations, group_by)


def inverse_simpson(counts: np.ndarray) -> float:
    """
    Calculate the Inverse Simpson diversity index.

    Defined as 1 / sum(p_i^2), where p_i are relative frequencies.
    Ranges from 1 (single dominant clone) to N (all clones equally abundant).

    Args:
        counts: Array of sequence counts

    Returns:
        Inverse Simpson index
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return 0.0
    frequencies = counts / total
    simpson = np.sum(frequencies ** 2)
    return 1.0 / simpson if simpson > 0 else 0.0


def inverse_simpson_polars(counts: pl.Series) -> float:
    """
    Calculate the Inverse Simpson diversity index from a Polars Series.

    Args:
        counts: Polars Series of sequence counts

    Returns:
        Inverse Simpson index
    """
    arr = counts.to_numpy().astype(float)
    total = arr.sum()
    if total == 0:
        return 0.0
    frequencies = arr / total
    simpson = np.sum(frequencies ** 2)
    return 1.0 / simpson if simpson > 0 else 0.0


def _clonality_polars(
    data: pl.DataFrame,
    rarefy: bool,
    min_count: Optional[int],
    iterations: int,
    group_by: str
) -> pl.DataFrame:
    """Calculate clonality using Polars."""

    # Ensure we have required columns
    required_cols = [group_by, "junction_aa", "duplicate_count"]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Convert duplicate_count to numeric if it's string
    if data[group_by].dtype == pl.Utf8:
        data = data.with_columns(
            pl.col("duplicate_count").cast(pl.Float64, strict=False).fill_null(1.0)
        )

    if not rarefy:
        # Standard diversity calculation - need to use apply for complex metrics
        # filtered_data = (
        #     data
        #     .filter(pl.col("productive") == True)
        #     .filter(pl.col("junction_aa") != "")
        # )
        filtered_data = data
        # Calculate using Python for complex metrics
        results = []
        for repertoire in filtered_data[group_by].unique():
            rep_data = filtered_data.filter(pl.col(group_by) == repertoire)
            counts = rep_data["duplicate_count"].to_numpy()
            frequencies = counts / counts.sum()

            # Calculate normalized clonality
            n_unique = rep_data["junction_aa"].n_unique()
            max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
            normalized_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0

            results.append({
                group_by: repertoire,
                "total_sequences": len(rep_data),
                "unique_sequences": n_unique,
                "total_count": int(counts.sum()),
                "clonality": 1 - normalized_entropy,
                "gini_coefficient": gini_coefficient(counts),
                "top_sequence": counts.max() / counts.sum(),
                "convergence": len(rep_data) / n_unique,
                "inverse_simpson": inverse_simpson_polars(rep_data["duplicate_count"]),
            })

        results = pl.DataFrame(results)
    else:
        # Rarefaction analysis
        if min_count is None:
            min_count = (
                data
                .group_by(group_by)
                .agg(pl.col("duplicate_count").sum().alias("total"))
                .select(pl.col("total").min())
                .item()
            )

        results = _perform_rarefaction_polars(data, min_count, iterations, group_by)

    return results


def _clonality_pandas(
    data: pd.DataFrame,
    rarefy: bool,
    min_count: Optional[int],
    iterations: int,
    group_by: str
) -> pd.DataFrame:
    """Calculate clonality using pandas."""

    # Ensure we have required columns
    required_cols = [group_by, "junction_aa", "duplicate_count"]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    productive_data = data
    #    # Filter productive sequences
    # productive_data = data[
    #     (data["productive"] == True) &
    #     (data["junction_aa"] != "") &
    #     (data["junction_aa"].notna())
    # ].copy()

    # Convert duplicate_count to numeric
    productive_data["duplicate_count"] = pd.to_numeric(
        productive_data["duplicate_count"], errors="coerce"
    ).fillna(1.0)

    if not rarefy:
        # Standard diversity calculation
        def calculate_metrics(group):
            counts = group["duplicate_count"].values
            frequencies = counts / counts.sum()
            n_unique = group["junction_aa"].nunique()

            # Calculate normalized clonality
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
            })

        results = productive_data.groupby(group_by).apply(calculate_metrics).reset_index()
    else:
        # Rarefaction analysis
        if min_count is None:
            min_count = productive_data.groupby(group_by)["duplicate_count"].sum().min()

        results = _perform_rarefaction_pandas(productive_data, min_count, iterations, group_by)

    return results


def _calculate_clonality_expr():
    """Polars expression for clonality calculation."""
    return (
        1 - (
            -pl.col("duplicate_count").map_elements(
                lambda counts: sum(
                    (c / sum(counts)) * np.log2(c / sum(counts))
                    for c in counts if c > 0
                ) if sum(counts) > 0 else 0,
                return_dtype=pl.Float64
            ) / np.log2(pl.len())
        )
    )


def _calculate_gini_expr():
    """Polars expression for Gini coefficient calculation."""
    return pl.col("duplicate_count").map_elements(
        lambda counts: gini_coefficient(counts),
        return_dtype=pl.Float64
    )


def _calculate_convergence_expr():
    """Polars expression for convergence calculation."""
    return pl.len() / pl.col("junction_aa").n_unique()


def _perform_rarefaction_polars(
    data: pl.DataFrame,
    min_count: int,
    iterations: int,
    group_by: str
) -> pl.DataFrame:
    """Perform rarefaction analysis using Polars."""

    results = []
    repertoires = data[group_by].unique().to_list()

    for repertoire in repertoires:
        repertoire_data = data.filter(pl.col(group_by) == repertoire)
        total_reads = repertoire_data["duplicate_count"].sum()

        # Skip if repertoire doesn't have enough reads
        if total_reads < min_count:
            continue

        # Prepare data for sampling - expand sequences by their counts
        sequences = repertoire_data["junction_aa"].to_list()
        counts = repertoire_data["duplicate_count"].to_list()

        iteration_metrics = []

        for _ in range(iterations):
            # Sample with replacement proportional to counts
            sampled_seqs = _weighted_sample(sequences, counts, min_count)

            # Count unique sequences and their frequencies
            unique_seqs, seq_counts = np.unique(sampled_seqs, return_counts=True)
            frequencies = seq_counts / seq_counts.sum()

            # Calculate metrics for this iteration
            n_unique = len(unique_seqs)
            max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
            normalized_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0

            metrics = {
                "total_sequences": len(sampled_seqs),
                "unique_sequences": n_unique,
                "total_count": min_count,
                "clonality": 1 - normalized_entropy,
                "gini_coefficient": gini_coefficient(seq_counts),
                "top_sequence": seq_counts.max() / seq_counts.sum(),
                "convergence": len(sampled_seqs) / n_unique,
                "inverse_simpson": inverse_simpson_polars(pl.Series(seq_counts)),
            }
            iteration_metrics.append(metrics)

        # Average across iterations
        avg_metrics = {
            group_by: repertoire,
            "total_sequences": np.mean([m["total_sequences"] for m in iteration_metrics]),
            "unique_sequences": np.mean([m["unique_sequences"] for m in iteration_metrics]),
            "total_count": min_count,
            "clonality": np.mean([m["clonality"] for m in iteration_metrics]),
            "gini_coefficient": np.mean([m["gini_coefficient"] for m in iteration_metrics]),
            "top_sequence": np.mean([m["top_sequence"] for m in iteration_metrics]),
            "convergence": np.mean([m["convergence"] for m in iteration_metrics]),
            "inverse_simpson": np.mean([m["inverse_simpson"] for m in iteration_metrics]),
        }
        results.append(avg_metrics)

    return pl.DataFrame(results)


def _weighted_sample(sequences: list, counts: list, sample_size: int) -> np.ndarray:
    """Sample sequences with replacement, weighted by their counts."""
    # Convert counts to probabilities
    total = sum(counts)
    probabilities = [c / total for c in counts]

    # Sample with replacement
    sampled_indices = np.random.choice(
        len(sequences),
        size=sample_size,
        replace=True,
        p=probabilities
    )

    return np.array([sequences[i] for i in sampled_indices])


def _perform_rarefaction_pandas(
    data: pd.DataFrame,
    min_count: int,
    iterations: int,
    group_by: str
) -> pd.DataFrame:
    """Perform rarefaction analysis using pandas."""

    results = []

    for repertoire, group in data.groupby(group_by):
        total_reads = group["duplicate_count"].sum()

        if total_reads >= min_count:
            # Prepare sequences and weights
            sequences = group["junction_aa"].values
            counts = group["duplicate_count"].values
            probabilities = counts / counts.sum()

            iteration_results = []

            for _ in range(iterations):
                # Sample sequences proportional to their counts
                sampled_seqs = np.random.choice(
                    sequences,
                    size=int(min_count),
                    replace=True,
                    p=probabilities
                )

                # Count unique sequences
                unique_seqs, seq_counts = np.unique(sampled_seqs, return_counts=True)
                frequencies = seq_counts / seq_counts.sum()

                # Calculate metrics for this iteration
                n_unique = len(unique_seqs)
                max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
                normalized_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0

                metrics = {
                    "total_sequences": len(sampled_seqs),
                    "unique_sequences": n_unique,
                    "total_count": min_count,
                    "clonality": 1 - normalized_entropy,
                    "gini_coefficient": gini_coefficient(seq_counts),
                    "top_sequence": seq_counts.max() / seq_counts.sum(),
                    "convergence": len(sampled_seqs) / n_unique,
                    "inverse_simpson": inverse_simpson(seq_counts),
                }
                iteration_results.append(metrics)

            # Average across iterations
            avg_metrics = pd.DataFrame(iteration_results).mean().to_dict()
            avg_metrics[group_by] = repertoire
            results.append(avg_metrics)

    return pd.DataFrame(results)


def diversity_metrics(
    data: Union[pl.DataFrame, pd.DataFrame],
    group_by: str = "repertoire_id"
) -> Union[pl.DataFrame, pd.DataFrame]:
    """
    Calculate comprehensive diversity metrics for repertoires.

    Args:
        data: Input data frame with AIRR-formatted sequences
        group_by: Column to group by

    Returns:
        Data frame with diversity metrics
    """
    return clonality(data, rarefy=False, group_by=group_by)


def common_sequences(
    data: Union[pl.DataFrame, pd.DataFrame],
    group_by: str = "repertoire_id",
    min_samples: int = 2
) -> Union[pl.DataFrame, pd.DataFrame]:
    """
    Find sequences that are common across multiple repertoires.

    Args:
        data: Input data frame with AIRR-formatted sequences
        group_by: Column to group by
        min_samples: Minimum number of samples a sequence must appear in

    Returns:
        Data frame with common sequences
    """
    is_polars = isinstance(data, pl.DataFrame)

    if is_polars:
        # Count how many repertoires each sequence appears in
        sequence_counts = (
            data
            .filter(pl.col("junction_aa") != "")
            .group_by("junction_aa")
            .agg([
                pl.col(group_by).n_unique().alias("repertoire_count"),
                pl.col("duplicate_count").sum().alias("total_count")
            ])
            .filter(pl.col("repertoire_count") >= min_samples)
            .sort("repertoire_count", descending=True)
        )
    else:
        # Pandas version
        sequence_counts = (
            data[data["junction_aa"] != ""]
            .groupby("junction_aa")
            .agg({
                group_by: "nunique",
                "duplicate_count": "sum"
            })
            .rename(columns={group_by: "repertoire_count"})
            .reset_index()
        )
        sequence_counts = sequence_counts[
            sequence_counts["repertoire_count"] >= min_samples
        ].sort_values("repertoire_count", ascending=False)

    return sequence_counts


def rarefaction_curve(
    data: Union[pl.DataFrame, pd.DataFrame],
    repertoire_ids: Optional[List[str]] = None,
    step_size: Optional[int] = None,
    iterations: int = 100,
    group_by: str = "repertoire_id"
) -> Union[pl.DataFrame, pd.DataFrame]:
    """
    Generate rarefaction curves for diversity analysis.

    Creates rarefaction curves by subsampling at multiple depths and calculating
    diversity metrics. This helps compare repertoires with different sequencing
    depths and assess sampling completeness.

    Args:
        data: Input DataFrame with AIRR-formatted sequences. Must contain
            'duplicate_count', 'junction_aa', and group_by columns.
        repertoire_ids: List of repertoire IDs to analyze. If None, uses all
            repertoires (default: None)
        step_size: Number of reads between sampling points. If None, uses
            10 steps across the range (default: None)
        iterations: Number of iterations for each sampling depth (default: 100)
        group_by: Column to group by (default: "repertoire_id")

    Returns:
        DataFrame with rarefaction curve data containing:
        - repertoire_id: Repertoire identifier
        - sample_size: Number of reads sampled
        - unique_sequences: Mean number of unique sequences
        - unique_sequences_sd: Standard deviation across iterations
        - clonality: Mean clonality at this depth
        - clonality_sd: Standard deviation of clonality

    Examples:
        >>> import lymphoseq as ls
        >>> data = ls.read_immunoseq("data/immunoseq/")

        >>> # Generate rarefaction curves for all samples
        >>> rarefaction = ls.rarefaction_curve(data, iterations=100)
        >>> print(rarefaction.head())

        >>> # Generate curves for specific samples
        >>> rarefaction = ls.rarefaction_curve(
        ...     data,
        ...     repertoire_ids=["S1", "S2"],
        ...     step_size=1000
        ... )

        >>> # Plot rarefaction curves (if using pandas)
        >>> import matplotlib.pyplot as plt
        >>> for rep in rarefaction["repertoire_id"].unique():
        ...     subset = rarefaction[rarefaction["repertoire_id"] == rep]
        ...     plt.plot(subset["sample_size"], subset["unique_sequences"], label=rep)
        >>> plt.xlabel("Sample size")
        >>> plt.ylabel("Unique sequences")
        >>> plt.legend()
        >>> plt.show()

    Notes:
        - Rarefaction helps compare diversity across samples with different depths
        - Sampling is done with replacement, weighted by sequence abundance
        - Curves that haven't plateaued suggest undersampling
        - Standard deviations indicate sampling variance

    See Also:
        - clonality(): Calculate diversity metrics with optional rarefaction
        - diversity_metrics(): Comprehensive diversity analysis
    """
    is_polars = isinstance(data, pl.DataFrame)

    # Filter for specified repertoires
    if repertoire_ids is not None:
        if is_polars:
            filtered = data.filter(pl.col(group_by).is_in(repertoire_ids))
        else:
            filtered = data[data[group_by].isin(repertoire_ids)].copy()
    else:
        filtered = data

    # Filter productive sequences
    if is_polars:
        filtered = filtered.filter(
            (pl.col("productive") == True) & (pl.col("junction_aa") != "")
        )
    else:
        filtered = filtered[
            (filtered["productive"] == True) &
            (filtered["junction_aa"] != "") &
            (filtered["junction_aa"].notna())
        ].copy()

    results = []

    # Process each repertoire
    if is_polars:
        repertoire_list = filtered[group_by].unique().to_list()
    else:
        repertoire_list = filtered[group_by].unique()

    for repertoire in repertoire_list:
        if is_polars:
            rep_data = filtered.filter(pl.col(group_by) == repertoire)
            sequences = rep_data["junction_aa"].to_list()
            counts = rep_data["duplicate_count"].to_list()
            total_reads = sum(counts)
        else:
            rep_data = filtered[filtered[group_by] == repertoire]
            sequences = rep_data["junction_aa"].values
            counts = rep_data["duplicate_count"].values
            total_reads = counts.sum()

        # Determine sampling points
        if step_size is None:
            # Use 10 evenly spaced points
            max_size = int(total_reads)
            step_size = max(1, max_size // 10)

        sample_sizes = list(range(step_size, int(total_reads) + 1, step_size))
        if total_reads not in sample_sizes:
            sample_sizes.append(int(total_reads))

        # Calculate probabilities for sampling
        probabilities = [c / total_reads for c in counts]

        # For each sampling depth
        for sample_size in sample_sizes:
            iteration_results = {
                "unique_sequences": [],
                "clonality": []
            }

            # Perform multiple iterations
            for _ in range(iterations):
                # Sample sequences
                sampled_seqs = np.random.choice(
                    sequences,
                    size=sample_size,
                    replace=True,
                    p=probabilities
                )

                # Count unique sequences
                unique_seqs, seq_counts = np.unique(sampled_seqs, return_counts=True)
                frequencies = seq_counts / seq_counts.sum()

                # Calculate metrics
                n_unique = len(unique_seqs)
                max_entropy = np.log2(n_unique) if n_unique > 1 else 1.0
                normalized_entropy = shannon_entropy(frequencies) / max_entropy if max_entropy > 0 else 0

                iteration_results["unique_sequences"].append(n_unique)
                iteration_results["clonality"].append(1 - normalized_entropy)

            # Calculate mean and SD across iterations
            results.append({
                group_by: repertoire,
                "sample_size": sample_size,
                "unique_sequences": np.mean(iteration_results["unique_sequences"]),
                "unique_sequences_sd": np.std(iteration_results["unique_sequences"]),
                "clonality": np.mean(iteration_results["clonality"]),
                "clonality_sd": np.std(iteration_results["clonality"])
            })

    # Convert to DataFrame
    if is_polars:
        result_df = pl.DataFrame(results)
        result_df = result_df.sort([group_by, "sample_size"])
    else:
        result_df = pd.DataFrame(results)
        result_df = result_df.sort_values([group_by, "sample_size"])

    return result_df