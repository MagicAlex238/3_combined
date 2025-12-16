import ast
import re
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    import plotly.express as px
    _HAS_PLOTLY = True
except Exception:
    _HAS_PLOTLY = False


def _ensure_list(x):
    """Convert x to a python list if it is a string representation of a list."""
    if pd.isna(x):
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, (set, tuple)):
        return list(x)
    # try to parse string-repr like "['As', 'Fe']" or "[As, Fe]" or "As, Fe"
    s = str(x).strip()
    # if it looks like a python list, literal_eval it
    if s.startswith("[") and s.endswith("]"):
        try:
            return list(ast.literal_eval(s))
        except Exception:
            pass
    # split on commas
    parts = [p.strip() for p in re.split(r",\s*", s) if p.strip() != ""]
    return parts


def normalize_metal_name(m: str) -> str:
    """Optional normalization of metal/electrolyte label (strip extra whitespace)."""
    if m is None:
        return ""
    return str(m).strip()


def explode_metals(df: pd.DataFrame,
                   metals_col: str = "consolidated_metals",
                   metal_out_col: str = "metal",
                   keep_original_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Explode the consolidated_metals column into one metal per row.
    Handles lists or string representations of lists.
    Returns a new DataFrame with column metal_out_col.
    """
    if metals_col not in df.columns:
        raise KeyError(f"{metals_col} not found in dataframe columns")

    df = df.copy()
    # convert to lists
    df[metals_col] = df[metals_col].apply(_ensure_list)
    # explode
    df_expl = df.explode(metals_col, ignore_index=True)
    # normalize names
    df_expl[metal_out_col] = df_expl[metals_col].apply(normalize_metal_name)
    # optional: drop rows where metal is empty string
    df_expl = df_expl[df_expl[metal_out_col] != ""].reset_index(drop=True)
    return df_expl


def compute_metal_counts(df_exploded: pd.DataFrame,
                         groupby_cols: List[str] = ["protein_name"],
                         metal_col: str = "metal",
                         distinct: bool = True) -> pd.DataFrame:
    """
    Compute number of distinct metals per group (default per protein_name).
    Returns df with columns groupby cols + metal_count + selectivity_label.
    """
    if metal_col not in df_exploded.columns:
        raise KeyError(f"{metal_col} not found in exploded dataframe")
    agg = df_exploded.groupby(groupby_cols)[metal_col].nunique().reset_index()
    agg = agg.rename(columns={metal_col: "metal_count"})
    # label: selective if metal_count == 1, promiscuous if >1
    agg["selectivity"] = np.where(agg["metal_count"] == 1, "SELECTIVE", "PROMISCUOUS")
    return agg


def presence_pivot(df_exploded: pd.DataFrame,
                   index_cols: List[str] = ["Genus", "protein_name"],
                   metal_col: str = "metal") -> pd.DataFrame:
    """
    Create a presence (0/1) pivot table: index=index_cols joined with ' - ' and columns=metals.
    Each cell indicates presence (1) of that metal for that protein/genus pair.
    """
    df = df_exploded.copy()
    # create a unique index label (protein-genus pair)
    df["_index_label"] = df[index_cols].astype(str).agg(" - ".join, axis=1)
    # presence
    df["_presence"] = 1
    pivot = df.drop_duplicates(subset=["_index_label", metal_col]).pivot_table(
        index="_index_label",
        columns=metal_col,
        values="_presence",
        aggfunc="max",
        fill_value=0
    )
    # sort metals by total prevalence descending
    pivot = pivot.loc[:, pivot.sum(axis=0).sort_values(ascending=False).index]
    return pivot


def abundance_pivot(df_exploded: pd.DataFrame,
                    value_col: str = "mean_cat1",
                    index_cols: List[str] = ["Genus", "protein_name"],
                    metal_col: str = "metal",
                    aggfunc = "median") -> pd.DataFrame:
    """
    Create pivot table of abundances (value_col) aggregated by index and metal.
    aggfunc can be 'mean', 'median', 'sum' or a callable.
    """
    if value_col not in df_exploded.columns:
        raise KeyError(f"{value_col} not in dataframe")
    df = df_exploded.copy()
    df["_index_label"] = df[index_cols].astype(str).agg(" - ".join, axis=1)
    pivot = df.pivot_table(index="_index_label", columns=metal_col, values=value_col, aggfunc=aggfunc, fill_value=0)
    # Optionally order columns
    pivot = pivot.loc[:, pivot.max(axis=0).sort_values(ascending=False).index]
    return pivot


# PLOTTING HELPERS

def plot_metal_counts(df_exploded: pd.DataFrame, metal_col: str = "metal", top_n: int = 30, figsize=(8,5)):
    """Bar chart: number of distinct proteins (or pairs) associated with each metal."""
    counts = df_exploded.groupby(metal_col)["_index_label" if "_index_label" in df_exploded.columns else metal_col].nunique()
    counts = counts.sort_values(ascending=False).head(top_n)
    plt.figure(figsize=figsize)
    sns.barplot(x=counts.values, y=counts.index, palette="viridis")
    plt.xlabel("Distinct protein / protein-Genus pairs")
    plt.ylabel("Metal")
    plt.title("Top metals by number of associated protein (or pair) entries")
    plt.tight_layout()


def plot_presence_heatmap(presence_df: pd.DataFrame, top_proteins: int = 50, figsize=(12,8), cmap="YlGnBu"):
    """
    presence_df: pivot table returned by presence_pivot (index = protein-genus label, columns = metal)
    Shows a heatmap of presence (0/1). Limit to top_proteins rows for readability.
    """
    if presence_df.shape[0] == 0:
        print("Empty presence dataframe")
        return
    sub = presence_df.copy()
    if sub.shape[0] > top_proteins:
        # choose the rows with highest metal diversity
        row_sums = sub.sum(axis=1).sort_values(ascending=False)
        sub = sub.loc[row_sums.index[:top_proteins]]
    plt.figure(figsize=figsize)
    sns.heatmap(sub, cmap=cmap, cbar_kws={'label': 'presence (1)/absence (0)'}, vmin=0, vmax=1)
    plt.xlabel("Metal")
    plt.ylabel("protein Genus - protein_name")
    plt.title(f"Presence heatmap ({sub.shape[0]} proteins shown)")
    plt.tight_layout()


def plot_abundance_dotplot(df_exploded: pd.DataFrame,
                           value_col: str = "mean_cat1",
                           groupby_col: str = "metal",
                           color_by: Optional[str] = None,
                           top_n_metals: int = 20,
                           figsize=(10,6)):
    """
    Dot plot with distribution of `value_col` per metal (uses seaborn stripplot + box/violin).
    """
    if value_col not in df_exploded.columns:
        raise KeyError(f"{value_col} not in dataframe")
    # take top metals by count
    metal_counts = df_exploded[groupby_col].value_counts().head(top_n_metals).index
    sub = df_exploded[df_exploded[groupby_col].isin(metal_counts)].copy()
    plt.figure(figsize=figsize)
    sns.boxplot(x=groupby_col, y=value_col, data=sub, order=metal_counts, showfliers=False, color="lightgray")
    sns.stripplot(x=groupby_col, y=value_col, data=sub, order=metal_counts, jitter=True, size=3,
                  hue=color_by, dodge=False, palette="tab10" if color_by else None)
    plt.xticks(rotation=45, ha="right")
    plt.title(f"Distribution of {value_col} by metal")
    plt.tight_layout()
    if color_by:
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')


def plot_selectivity_vs_abundance(selectivity_df: pd.DataFrame,
                                  exploded_df: pd.DataFrame,
                                  value_col: str = "mean_cat1",
                                  groupby_cols: List[str] = ["Genus", "protein_name"],
                                  figsize=(8,5)):
    """
    Compare abundance distributions for SELECTIVE vs PROMISCUOUS proteins.
    selectivity_df is output of compute_metal_counts.
    exploded_df should be used to map values.
    """
    # merge selectivity label back to exploded rows per protein (or per pair)
    if len(groupby_cols) == 1 and groupby_cols[0] == "protein_name":
        merged = exploded_df.merge(selectivity_df[["protein_name","selectivity"]], on="protein_name", how="left")
    else:
        merged = exploded_df.copy()
        merged["_index_label"] = merged[groupby_cols].astype(str).agg(" - ".join, axis=1)
        sel = selectivity_df.copy()
        sel["_index_label"] = sel[groupby_cols].astype(str).agg(" - ".join, axis=1)
        merged = merged.merge(sel[["_index_label","selectivity"]], on="_index_label", how="left")
    if value_col not in merged.columns:
        raise KeyError(f"{value_col} not in dataframe")
    plt.figure(figsize=figsize)
    sns.violinplot(x="selectivity", y=value_col, data=merged, inner="quartile", palette="Set2")
    plt.title(f"{value_col} distribution: SELECTIVE vs PROMISCUOUS")
    plt.tight_layout()


def interactive_sunburst(df_exploded: pd.DataFrame,
                         path_cols: List[str] = ["Genus", "protein_name", "metal"],
                         values_col: Optional[str] = None):
    """
    Interactive sunburst using plotly (if installed). Shows hierarchy Genus -> protein -> metal.
    values_col can be e.g., 'norm_abund_contri' to size slices.
    """
    if not _HAS_PLOTLY:
        raise RuntimeError("plotly is not installed. pip install plotly to use interactive_sunburst")
    df = df_exploded.copy()
    if values_col and values_col in df.columns:
        agg = df.groupby(path_cols)[values_col].sum().reset_index()
        fig = px.sunburst(agg, path=path_cols, values=values_col)
    else:
        agg = df.groupby(path_cols).size().reset_index(name="count")
        fig = px.sunburst(agg, path=path_cols, values="count")
    fig.update_layout(margin=dict(t=10, l=10, r=10, b=10))
    return fig


# Utility: list numeric columns in dataframe
def numeric_columns(df: pd.DataFrame) -> List[str]:
    """Return list of columns in df that are numeric and suitable for plotting as numeric data."""
    return df.select_dtypes(include=[np.number]).columns.tolist()


if __name__ == "__main__":
    print("metal_selectivity_analysis module loaded.")