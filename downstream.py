"""Downstream analysis utilities for the bulk RNA-seq workflow."""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
from plotly.graph_objs import Figure
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import linregress
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from . import config

try:  # pragma: no cover - optional dependency
    import gseapy as gp
except ImportError:  # pragma: no cover - handled gracefully
    gp = None


@dataclass
class Deseq2Result:
    result_file: Path
    normalized_counts: Path
    rlog_counts: Path


class DownstreamAnalysisError(RuntimeError):
    """Raised when downstream analysis fails."""


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_counts(dataset: str) -> pd.DataFrame:
    counts_path = config.PROCESSED_DATA_DIR / dataset / "counts.csv"
    if not counts_path.exists():
        raise DownstreamAnalysisError(
            f"Counts matrix was not found at '{counts_path}'. Run raw processing first."
        )
    return pd.read_csv(counts_path, index_col=0)


def load_metadata(dataset: str) -> pd.DataFrame:
    metadata_path = config.RAW_DATA_DIR / dataset / "metadata.csv"
    if not metadata_path.exists():
        raise DownstreamAnalysisError(
            f"Metadata file '{metadata_path}' was not found."
        )
    return pd.read_csv(metadata_path, index_col=0)


# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------

def run_deseq2(
    dataset: str,
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design_formula: str,
    condition_column: str,
    control_group: str,
    treatment_group: str,
    output_dir: Optional[Path] = None,
) -> Deseq2Result:
    if output_dir is None:
        output_dir = config.RESULTS_DIR / dataset / "deseq2"
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile("w", suffix="_counts.csv", delete=False) as counts_handle:
        counts.to_csv(counts_handle.name)
        counts_path = Path(counts_handle.name)

    with tempfile.NamedTemporaryFile("w", suffix="_metadata.csv", delete=False) as metadata_handle:
        metadata.to_csv(metadata_handle.name)
        metadata_path = Path(metadata_handle.name)

    result_file = output_dir / "deseq2_results.csv"
    normalized_counts = output_dir / "normalized_counts.csv"
    rlog_counts = output_dir / "rlog_counts.csv"

    cmd = [
        "Rscript",
        str(config.DESEQ2_SCRIPT),
        "--counts",
        str(counts_path),
        "--metadata",
        str(metadata_path),
        "--design",
        design_formula,
        "--condition",
        condition_column,
        "--control",
        control_group,
        "--treatment",
        treatment_group,
        "--output",
        str(result_file),
        "--normalized",
        str(normalized_counts),
        "--rlog",
        str(rlog_counts),
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:  # pragma: no cover - error path
        raise DownstreamAnalysisError(
            f"DESeq2 failed with exit code {exc.returncode}. Check R installation and inputs."
        ) from exc
    finally:
        counts_path.unlink(missing_ok=True)
        metadata_path.unlink(missing_ok=True)

    if not result_file.exists():  # pragma: no cover - defensive
        raise DownstreamAnalysisError("DESeq2 did not produce the expected result file.")

    return Deseq2Result(result_file=result_file, normalized_counts=normalized_counts, rlog_counts=rlog_counts)


# ---------------------------------------------------------------------------
# Visualization helpers
# ---------------------------------------------------------------------------

def perform_pca(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    top_n: int = 500,
    scale: bool = True,
    color_by: Optional[str] = None,
    hover_data: Optional[List[str]] = None,
) -> Figure:
    counts = counts.loc[counts.var(axis=1).sort_values(ascending=False).index[:top_n]]
    matrix = counts.to_numpy().T
    if scale:
        matrix = StandardScaler().fit_transform(matrix)

    pca = PCA(n_components=2)
    components = pca.fit_transform(matrix)
    explained = pca.explained_variance_ratio_ * 100

    pca_df = metadata.copy()
    pca_df["PC1"] = components[:, 0]
    pca_df["PC2"] = components[:, 1]
    hover_cols = hover_data if hover_data else metadata.columns.tolist()

    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color=color_by if color_by in metadata.columns else None,
        hover_data=hover_cols,
        labels={"PC1": f"PC1 ({explained[0]:.1f}% explained)", "PC2": f"PC2 ({explained[1]:.1f}% explained)"},
    )
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color="DarkSlateGrey")))
    return fig


def generate_volcano_plot(
    deseq_results: pd.DataFrame,
    log2fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
    max_log2fc: Optional[float] = None,
    color_up: str = "#d73027",
    color_down: str = "#4575b4",
    color_ns: str = "#808080",
) -> Figure:
    df = deseq_results.copy()
    df["-log10(padj)"] = -np.log10(df["padj"].replace(0, np.nan))
    if max_log2fc is not None:
        df["log2FoldChange"] = df["log2FoldChange"].clip(-max_log2fc, max_log2fc)

    conditions = [
        (df["padj"] < padj_threshold) & (df["log2FoldChange"] > log2fc_threshold),
        (df["padj"] < padj_threshold) & (df["log2FoldChange"] < -log2fc_threshold),
    ]
    choices = ["up", "down"]
    df["regulation"] = np.select(conditions, choices, default="ns")

    color_discrete_map = {"up": color_up, "down": color_down, "ns": color_ns}
    fig = px.scatter(
        df.reset_index().rename(columns={"index": "gene"}),
        x="log2FoldChange",
        y="-log10(padj)",
        color="regulation",
        hover_name="gene",
        hover_data={"padj": True, "log2FoldChange": True},
        color_discrete_map=color_discrete_map,
    )
    fig.add_vline(x=log2fc_threshold, line_dash="dash", line_color="black")
    fig.add_vline(x=-log2fc_threshold, line_dash="dash", line_color="black")
    fig.add_hline(y=-np.log10(padj_threshold), line_dash="dash", line_color="black")
    fig.update_layout(xaxis_title="log2 Fold Change", yaxis_title="-log10 adjusted p-value")
    return fig


def generate_heatmap(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    top_n: int = 50,
    color_scale: str = "RdBu",
) -> Figure:
    top_genes = counts.var(axis=1).sort_values(ascending=False).index[:top_n]
    sub_counts = counts.loc[top_genes]
    z = StandardScaler().fit_transform(sub_counts.T).T

    row_linkage = linkage(pdist(z), method="average")
    row_order = leaves_list(row_linkage)
    ordered_counts = pd.DataFrame(z[row_order], index=sub_counts.index[row_order], columns=sub_counts.columns)
    shared_samples = [col for col in metadata.index if col in ordered_counts.columns]
    if shared_samples:
        ordered_counts = ordered_counts[shared_samples]

    figure = px.imshow(
        ordered_counts,
        color_continuous_scale=color_scale,
        aspect="auto",
        labels=dict(x="Sample", y="Gene", color="Z-score"),
    )
    return figure


# ---------------------------------------------------------------------------
# Enrichment analysis
# ---------------------------------------------------------------------------

def run_gsea(
    dataset: str,
    ranking: pd.Series,
    gene_set: str,
    outdir: Optional[Path] = None,
    min_size: int = 15,
    max_size: int = 500,
) -> Path:
    if gp is None:  # pragma: no cover - optional dependency
        raise DownstreamAnalysisError("gseapy is not installed. Install it to perform GSEA analysis.")

    if outdir is None:
        outdir = config.RESULTS_DIR / dataset / "gsea"
    outdir.mkdir(parents=True, exist_ok=True)

    gene_set_path = config.GENE_SET_DIR / gene_set
    if not gene_set_path.exists():
        raise DownstreamAnalysisError(
            f"Gene set '{gene_set}' was not found in '{config.GENE_SET_DIR}'."
        )

    gp.prerank(
        rnk=ranking,
        gene_sets=str(gene_set_path),
        outdir=str(outdir),
        min_size=min_size,
        max_size=max_size,
        seed=42,
        verbose=True,
    )
    return outdir


# ---------------------------------------------------------------------------
# Trend analysis
# ---------------------------------------------------------------------------

def compute_trends(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    numeric_column: str,
    top_n: int = 50,
) -> pd.DataFrame:
    if numeric_column not in metadata.columns:
        raise DownstreamAnalysisError(
            f"Metadata column '{numeric_column}' is not available for trend analysis."
        )
    values = metadata[numeric_column].astype(float)
    trends: List[Dict[str, float]] = []
    for gene, expression in counts.iterrows():
        slope, intercept, rvalue, pvalue, stderr = linregress(values, expression)
        trends.append(
            {
                "gene": gene,
                "slope": slope,
                "rvalue": rvalue,
                "pvalue": pvalue,
                "stderr": stderr,
            }
        )

    trend_df = pd.DataFrame(trends).sort_values("pvalue").head(top_n)
    return trend_df


__all__ = [
    "Deseq2Result",
    "DownstreamAnalysisError",
    "load_counts",
    "load_metadata",
    "run_deseq2",
    "perform_pca",
    "generate_volcano_plot",
    "generate_heatmap",
    "run_gsea",
    "compute_trends",
]
