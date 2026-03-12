"""
data_preprocessing.py
---------------------
Data loading, cleaning, and feature engineering for the AML TCGA Pan-Can Atlas 2018 dataset.

Real file structure (verified):
- data_clinical_patient.txt : PATIENT_ID, SUBTYPE, OS_STATUS, OS_MONTHS only (AGE/FAB not available)
- data_clinical_sample.txt  : PATIENT_ID, SAMPLE_ID, TMB_NONSYNONYMOUS, ANEUPLOIDY_SCORE, etc.
- data_mutations.txt        : MAF format — Hugo_Symbol, Tumor_Sample_Barcode, t_alt_count, t_ref_count
- data_mrna_seq_v2_rsem.txt : genes x patients RSEM expression matrix
"""

import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"


def load_clinical(filepath: str = None) -> pd.DataFrame:
    """
    Load and clean the clinical patient metadata file.

    Available columns in the Pan-Can Atlas version:
        PATIENT_ID, SUBTYPE, OS_STATUS, OS_MONTHS
    Note: AGE, SEX, FAB_TYPE are NOT available in this dataset version.

    Returns a DataFrame with binary OS_STATUS column.
    """
    path = filepath or DATA_DIR / "data_clinical_patient.txt"
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    df.columns = [c.strip().upper() for c in df.columns]

    # Parse OS_STATUS to binary (0 = alive, 1 = deceased)
    df["OS_STATUS"] = df["OS_STATUS"].apply(
        lambda x: 1 if str(x).startswith("1") or "DECEASED" in str(x).upper() else 0
    )

    # Parse OS_MONTHS — TCGA may encode missing as [Not Available]
    df["OS_MONTHS"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")

    # Drop rows with missing survival data
    df = df.dropna(subset=["OS_MONTHS", "OS_STATUS"])
    df = df[df["OS_MONTHS"] > 0].reset_index(drop=True)

    return df


def load_sample(filepath: str = None) -> pd.DataFrame:
    """
    Load the clinical sample metadata file.

    Available columns: PATIENT_ID, SAMPLE_ID, TMB_NONSYNONYMOUS,
                       ANEUPLOIDY_SCORE, ONCOTREE_CODE, etc.
    """
    path = filepath or DATA_DIR / "data_clinical_sample.txt"
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    df.columns = [c.strip().upper() for c in df.columns]

    for col in ["TMB_NONSYNONYMOUS", "ANEUPLOIDY_SCORE"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def load_clinical_full(patient_path: str = None, sample_path: str = None) -> pd.DataFrame:
    """
    Load and merge clinical patient + sample files.
    Recommended entry point — gives OS, SUBTYPE, TMB, and aneuploidy score.
    """
    patient = load_clinical(patient_path)
    sample  = load_sample(sample_path)

    keep = [c for c in ["PATIENT_ID", "TMB_NONSYNONYMOUS", "ANEUPLOIDY_SCORE",
                         "ONCOTREE_CODE", "SAMPLE_ID"] if c in sample.columns]
    return patient.merge(sample[keep], on="PATIENT_ID", how="left")


def load_mutations(filepath: str = None) -> pd.DataFrame:
    """
    Load somatic mutation data (MAF format).
    - Computes VAF from allele counts
    - Removes Silent mutations
    - Adds normalized PATIENT_ID (12-char barcode)
    """
    path = filepath or DATA_DIR / "data_mutations.txt"
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
    df["t_ref_count"] = pd.to_numeric(df["t_ref_count"], errors="coerce")
    total = df["t_alt_count"] + df["t_ref_count"]
    df["VAF"] = df["t_alt_count"] / total.replace(0, np.nan)

    df["PATIENT_ID"] = df["Tumor_Sample_Barcode"].str[:12]
    df = df[df["Variant_Classification"] != "Silent"].copy()

    return df


def add_mutation_flags(clinical_df: pd.DataFrame,
                       mutations_df: pd.DataFrame,
                       genes: list = None) -> pd.DataFrame:
    """
    Add binary {GENE}_mut columns to clinical DataFrame.
    Defaults to top 10 AML driver genes.
    """
    if genes is None:
        genes = ["FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2",
                 "RUNX1", "TP53", "CEBPA", "TET2", "NRAS"]

    df = clinical_df.copy()
    for gene in genes:
        pts = mutations_df[mutations_df["Hugo_Symbol"] == gene]["PATIENT_ID"].unique()
        df[f"{gene}_mut"] = df["PATIENT_ID"].isin(pts).astype(int)
    return df


def get_mutation_counts(mutations_df: pd.DataFrame) -> pd.Series:
    """Count non-silent mutations per patient. Returns pd.Series indexed by PATIENT_ID."""
    if "PATIENT_ID" not in mutations_df.columns:
        mutations_df = mutations_df.copy()
        mutations_df["PATIENT_ID"] = mutations_df["Tumor_Sample_Barcode"].str[:12]
    return mutations_df.groupby("PATIENT_ID").size().rename("mutation_count")


def get_gene_mutation_matrix(mutations_df: pd.DataFrame,
                              gene_col: str = "Hugo_Symbol") -> pd.DataFrame:
    """Binary patient x gene mutation matrix. Rows = patients, columns = genes."""
    if "PATIENT_ID" not in mutations_df.columns:
        mutations_df = mutations_df.copy()
        mutations_df["PATIENT_ID"] = mutations_df["Tumor_Sample_Barcode"].str[:12]
    mutations_df = mutations_df.copy()
    mutations_df["mutated"] = 1
    return (mutations_df
            .groupby(["PATIENT_ID", gene_col])["mutated"]
            .max()
            .unstack(fill_value=0))


def load_expression(filepath: str = None, log_transform: bool = True) -> pd.DataFrame:
    """
    Load RNA-seq RSEM gene expression. Genes as rows, patients as columns.
    Applies log2(x+1) by default.
    """
    path = filepath or DATA_DIR / "data_mrna_seq_v2_rsem.txt"
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    if "Entrez_Gene_Id" in df.columns:
        df = df.drop(columns=["Entrez_Gene_Id"])
    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if log_transform:
        df = np.log2(df + 1)
    return df


def select_highly_variable_genes(expression_df: pd.DataFrame,
                                  n_genes: int = 2000) -> pd.DataFrame:
    """Filter to top n most variable genes. Standard step before PCA/clustering."""
    return expression_df.loc[expression_df.var(axis=1).nlargest(n_genes).index]
