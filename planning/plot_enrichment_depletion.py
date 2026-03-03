import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgba

# IMGT positions for CDR3 (equivalent to R's imgt factor)
IMGT_POSITIONS = (
    [str(x) for x in range(109, 112)] +
    [f"111.{x}" for x in range(1, 5)] +
    [f"112.{x}" for x in range(5, 0, -1)] +
    [str(x) for x in range(112, 115)]
)

# Amino acid -> colour mapping (equivalent to kidera_hydro_aacids)
KIDERA_HYDRO_AACIDS = {
    'A': 'green', 'R': 'red',   'N': 'orange', 'D': 'orange',
    'C': 'blue',  'Q': 'red',   'E': 'red',    'G': 'green',
    'H': 'orange','I': 'green', 'L': 'blue',   'K': 'red',
    'M': 'green', 'F': 'blue',  'P': 'green',  'S': 'orange',
    'T': 'orange','W': 'blue',  'Y': 'green',  'V': 'green'
}


def _get_avail_aa_length(lengths):
    """Build availability matrix for type='length' canvas (lengths 9-22)."""
    avail = []
    for z in range(1, 8):  # z = 1..7 maps to lengths 9..15 and 16..22
        row_even = [1]*z + [0]*(15 - 2*z) + [1]*z
        row_odd  = [1]*z + [0]*(15 - 2*z - 1) + [1]*(z + 1)
        avail.extend(row_even)
        avail.extend(row_odd)
    return avail  # 14 * 15 = 210 values


def _get_avail_aa_len_bin():
    """Build availability matrix for type='len_bin' canvas."""
    bins = [
        [1]*2  + [0]*11 + [1]*2,
        [1]*3  + [0]*9  + [1]*3,
        [1]*4  + [0]*7  + [1]*4,
        [1]*5  + [0]*5  + [1]*5,
        [1]*15,
    ]
    return [val for row in bins for val in row]


def plot_canvas_shade(ax, length_type="length", shade_col="red"):
    """
    Draw the empty tiled canvas with shading for unavailable positions.
    Equivalent to R's plot_canvas_shade().
    """
    shade_rgba = to_rgba(shade_col, alpha=0.2)

    if length_type == "length":
        lengths = [f"l{x}" for x in range(9, 23)]
        avail_aa = _get_avail_aa_length(lengths)
        canvas = pd.DataFrame({
            'imgt':   IMGT_POSITIONS * 14,
            'length': [l for l in lengths for _ in IMGT_POSITIONS],
            'avail':  avail_aa
        })
    elif length_type == "len_bin":
        bins = ["(0,11]", "(11,13]", "(13,15]", "(15,17]", "(17,28]"]
        avail_aa = _get_avail_aa_len_bin()
        canvas = pd.DataFrame({
            'imgt':   IMGT_POSITIONS * 5,
            'length': [b for b in bins for _ in IMGT_POSITIONS],
            'avail':  avail_aa
        })
    else:
        raise ValueError(f"Unknown length_type: {length_type}")

    imgt_order = {pos: i for i, pos in enumerate(IMGT_POSITIONS)}
    length_order = canvas['length'].unique().tolist()
    len_idx = {l: i for i, l in enumerate(length_order)}

    for _, row in canvas.iterrows():
        x = imgt_order[row['imgt']]
        y = len_idx[row['length']]
        color = 'white' if row['avail'] == 1 else shade_rgba
        rect = mpatches.FancyBboxPatch(
            (x - 0.5, y - 0.5), 1, 1,
            boxstyle="square,pad=0",
            linewidth=1, edgecolor='darkgray', facecolor=color
        )
        ax.add_patch(rect)

    ax.set_xlim(-0.5, len(IMGT_POSITIONS) - 0.5)
    ax.set_ylim(-0.5, len(length_order) - 0.5)
    ax.set_xticks(range(len(IMGT_POSITIONS)))
    ax.set_xticklabels(IMGT_POSITIONS, rotation=90)
    ax.set_yticks(range(len(length_order)))
    ax.set_yticklabels(length_order)

    return ax, length_order, imgt_order


def plot_enrichment_depletion(
    ready_extracted: pd.DataFrame,
    length_col: str = "length",
    estim_col: str = "estim_hom_vs_DQ2DQ8",
    estim_expr_for_plotting: str = "estim_hom_vs_DQ2DQ8 < 0",
    aminoacid_mapping: dict = None,
    inset_colour: str = "red",
    figsize: tuple = (14, 8)
):
    """
    Python equivalent of R's plot_enrichment_depletion2().

    Parameters
    ----------
    ready_extracted : pd.DataFrame
        DataFrame loaded from e.g. all_together.csv.gz.
        Expected columns: 'length', 'cells', 'imgt', 'aa', and the estim_col.
    length_col : str
        'length' or 'len_bin' — controls canvas type.
    estim_col : str
        Column name used for filtering.
    estim_expr_for_plotting : str
        Pandas query string for filtering rows to plot.
    aminoacid_mapping : dict
        Dict mapping single-letter AA -> colour. Defaults to KIDERA_HYDRO_AACIDS.
    inset_colour : str
        Shade colour for unavailable positions.
    figsize : tuple
        Figure size.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    if aminoacid_mapping is None:
        aminoacid_mapping = KIDERA_HYDRO_AACIDS

    fig, ax = plt.subplots(figsize=figsize)
    ax, length_order, imgt_order = plot_canvas_shade(ax, length_type=length_col, shade_col=inset_colour)

    len_idx = {l: i for i, l in enumerate(length_order)}

    # Filter rows matching the expression
    filtered = ready_extracted.query(estim_expr_for_plotting).copy()

    if filtered.empty:
        print("No data matches the filter expression.")
        return fig, ax

    # Get all unique colours present in this data
    all_colours = set()
    for aa_str in filtered['aa'].dropna():
        for aa in aa_str:
            if aa in aminoacid_mapping:
                all_colours.add(aminoacid_mapping[aa])

    # One pass per colour layer — AAs of other colours replaced by spaces
    for colour in all_colours:
        other_aas = {aa for aa, col in aminoacid_mapping.items() if col != colour}

        layer = filtered.copy()
        # Replace AAs of other colours with spaces
        layer['aa'] = layer['aa'].apply(
            lambda s: ''.join(' ' if c in other_aas else c for c in str(s))
        )
        layer['aa_col'] = colour

        grouped = (
            layer.groupby(['length', 'cells', 'imgt'], as_index=False)
            .agg(aa=('aa', lambda x: ''.join(x)), aa_col=('aa_col', 'first'))
        )

        for _, row in grouped.iterrows():
            if row['imgt'] not in imgt_order or row['length'] not in len_idx:
                continue
            x = imgt_order[row['imgt']]
            y = len_idx[row['length']]
            ax.text(
                x, y, row['aa'].strip(),
                color=colour, fontsize=9, fontweight='bold',
                ha='center', va='center'
            )

    # Legend
    legend_handles = [
        mpatches.Patch(color=col, label=col) for col in sorted(all_colours)
    ]
    ax.legend(handles=legend_handles, title="AA property", bbox_to_anchor=(1.01, 1), loc='upper left')

    ax.set_title(f"Enrichment/Depletion: {estim_expr_for_plotting}")
    plt.tight_layout()
    return fig, ax


def compute_inverse_simpson(
    df: pd.DataFrame,
    freq_col: str = "freq",
    group_cols: list = None
) -> pd.DataFrame:
    """
    Compute the Inverse Simpson index (1 / sum(p_i^2)).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing frequency/proportion data.
        If values are counts rather than proportions, they will be normalised
        within each group automatically.
    freq_col : str
        Column name containing frequencies or counts.
    group_cols : list, optional
        Columns to group by before computing (e.g. ['sample', 'length']).
        If None, computes a single value for the whole DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with group columns (if any) and 'inverse_simpson' column.
    """
    def _inv_simpson(series: pd.Series) -> float:
        p = series / series.sum()          # normalise to proportions
        return 1.0 / (p ** 2).sum()

    if group_cols:
        result = (
            df.groupby(group_cols)[freq_col]
            .apply(_inv_simpson)
            .reset_index()
            .rename(columns={freq_col: "inverse_simpson"})
        )
    else:
        result = pd.DataFrame({"inverse_simpson": [_inv_simpson(df[freq_col])]})

    return result


# Example usage:
# df = pd.read_csv("all_together.csv.gz")
# fig, ax = plot_enrichment_depletion(df)
# plt.show()
