"""
visualization.py
----------------
Reusable plotting functions for the AML TCGA analysis.
Covers mutation landscape, PCA, clustering, and general utilities.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# ── Color palette ─────────────────────────────────────────────────────────────
PALETTE = {
    "primary": "#2E75B6",
    "danger": "#C0392B",
    "success": "#27AE60",
    "warning": "#E67E22",
    "neutral": "#7F8C8D",
    "light_bg": "#F4F6F9",
}

MUTATION_TYPE_COLORS = {
    "Missense_Mutation": "#3498DB",
    "Nonsense_Mutation": "#E74C3C",
    "Frame_Shift_Del": "#E67E22",
    "Frame_Shift_Ins": "#F39C12",
    "Splice_Site": "#9B59B6",
    "In_Frame_Del": "#1ABC9C",
    "In_Frame_Ins": "#2ECC71",
    "Translation_Start_Site": "#34495E",
}


def plot_top_mutated_genes(mutations_df: pd.DataFrame,
                            gene_col: str = "Hugo_Symbol",
                            patient_col: str = "Tumor_Sample_Barcode",
                            n_genes: int = 15,
                            n_patients: int = None,
                            ax=None):
    """
    Bar chart of the top N most frequently mutated genes.
    Y-axis shows mutation frequency (%) across the cohort.

    Parameters
    ----------
    n_patients : total number of patients in cohort (for frequency calculation)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))

    mutations_df = mutations_df.copy()
    mutations_df["PATIENT_ID"] = mutations_df[patient_col].str[:12]

    gene_counts = (
        mutations_df.groupby(gene_col)["PATIENT_ID"]
        .nunique()
        .sort_values(ascending=False)
        .head(n_genes)
    )

    if n_patients:
        freq = (gene_counts / n_patients * 100).round(1)
        ylabel = "Mutation frequency (%)"
    else:
        freq = gene_counts
        ylabel = "Number of patients"

    colors = [PALETTE["primary"]] * len(freq)
    bars = ax.bar(freq.index, freq.values, color=colors, edgecolor="white", linewidth=0.5)

    # Annotate bars
    for bar, val in zip(bars, freq.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                f"{val:.0f}{'%' if n_patients else ''}", ha="center", va="bottom", fontsize=9)

    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"Top {n_genes} Mutated Genes — AML TCGA Cohort", fontsize=14, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor(PALETTE["light_bg"])

    return ax


def plot_vaf_distribution(mutations_df: pd.DataFrame,
                           vaf_col: str = "VAF",
                           gene_col: str = "Hugo_Symbol",
                           top_genes: list = None,
                           ax=None):
    """
    Violin/box plot of VAF distribution across top mutated genes.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))

    df = mutations_df.dropna(subset=[vaf_col]).copy()

    if top_genes is None:
        top_genes = df[gene_col].value_counts().head(10).index.tolist()

    df = df[df[gene_col].isin(top_genes)]

    sns.boxplot(data=df, x=gene_col, y=vaf_col, order=top_genes,
                palette="Blues", ax=ax, linewidth=0.8)
    sns.stripplot(data=df, x=gene_col, y=vaf_col, order=top_genes,
                  color="black", alpha=0.3, size=2.5, ax=ax)

    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Variant Allele Frequency (VAF)", fontsize=12)
    ax.set_title("VAF Distribution by Gene", fontsize=14, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    ax.axhline(0.4, linestyle="--", color="red", alpha=0.5, label="VAF = 0.4 (clonal threshold)")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax


def plot_pca(expression_matrix: pd.DataFrame,
             color_labels: pd.Series = None,
             n_components: int = 2,
             ax=None,
             title: str = "PCA — Gene Expression Profiles"):
    """
    PCA scatter plot on gene expression data.

    Parameters
    ----------
    expression_matrix : genes × patients DataFrame (already log2-transformed, HVG-filtered)
    color_labels : optional pd.Series (patient-indexed) for point coloring

    Returns
    -------
    pca object, ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # Transpose: samples as rows, genes as columns
    X = expression_matrix.T
    X.columns = X.columns.astype(str)
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=max(n_components, 10))
    coords = pca.fit_transform(X_scaled)

    explained = pca.explained_variance_ratio_ * 100

    if color_labels is not None:
        # Align labels to sample order
        aligned = color_labels.reindex(X.index)
        unique_labels = aligned.dropna().unique()
        cmap = plt.cm.get_cmap("tab10", len(unique_labels))
        label_to_color = {lbl: cmap(i) for i, lbl in enumerate(unique_labels)}

        for lbl in unique_labels:
            mask = aligned == lbl
            ax.scatter(coords[mask, 0], coords[mask, 1],
                       label=str(lbl), alpha=0.7, s=40, color=label_to_color[lbl])
        ax.legend(title=color_labels.name, fontsize=9, title_fontsize=10)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], alpha=0.6, s=40, color=PALETTE["primary"])

    ax.set_xlabel(f"PC1 ({explained[0]:.1f}% variance)", fontsize=12)
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}% variance)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return pca, ax


def plot_mutation_count_distribution(mutation_counts: pd.Series, ax=None):
    """
    Histogram of per-patient mutation counts with outlier annotation.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))

    mean = mutation_counts.mean()
    std = mutation_counts.std()
    threshold = mean + 3 * std

    ax.hist(mutation_counts, bins=30, color=PALETTE["primary"], edgecolor="white", linewidth=0.5)
    ax.axvline(mean, color="grey", linestyle="--", label=f"Mean = {mean:.1f}")
    ax.axvline(threshold, color=PALETTE["danger"], linestyle="--",
               label=f"Outlier threshold (+3σ = {threshold:.0f})")

    n_outliers = (mutation_counts > threshold).sum()
    ax.text(0.97, 0.92, f"Hypermutated outliers: {n_outliers}",
            transform=ax.transAxes, ha="right", fontsize=10,
            color=PALETTE["danger"], fontweight="bold")

    ax.set_xlabel("Number of somatic mutations per patient", fontsize=12)
    ax.set_ylabel("Number of patients", fontsize=12)
    ax.set_title("Somatic Mutational Burden — AML TCGA Cohort", fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax
