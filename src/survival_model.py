"""
survival_model.py
-----------------
Survival analysis utilities for the AML TCGA dataset.
Implements Kaplan-Meier estimation, log-rank testing, and Cox PH regression.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test


def kaplan_meier_overall(clinical_df: pd.DataFrame,
                          duration_col: str = "OS_MONTHS",
                          event_col: str = "OS_STATUS",
                          ax=None,
                          title: str = "Overall Survival — AML TCGA Cohort"):
    """
    Fit and plot a Kaplan-Meier curve for the full cohort.

    Parameters
    ----------
    clinical_df : DataFrame with duration and event columns
    ax : matplotlib Axes (created if None)

    Returns
    -------
    KaplanMeierFitter object, matplotlib Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5))

    kmf = KaplanMeierFitter()
    kmf.fit(
        durations=clinical_df[duration_col],
        event_observed=clinical_df[event_col],
        label="All patients"
    )
    kmf.plot_survival_function(ax=ax, ci_show=True, color="#2E75B6")

    median_os = kmf.median_survival_time_
    ax.axvline(median_os, linestyle="--", color="grey", alpha=0.6,
               label=f"Median OS = {median_os:.1f} months")
    ax.set_xlabel("Time (months)", fontsize=12)
    ax.set_ylabel("Survival probability", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return kmf, ax


def kaplan_meier_stratified(clinical_df: pd.DataFrame,
                             group_col: str,
                             duration_col: str = "OS_MONTHS",
                             event_col: str = "OS_STATUS",
                             group_labels: dict = None,
                             colors: list = None,
                             ax=None,
                             title: str = None):
    """
    Fit and plot stratified Kaplan-Meier curves with log-rank test.

    Parameters
    ----------
    group_col : binary column (0/1) used to split the cohort
    group_labels : dict mapping group values to display names
    colors : list of colors for each group

    Returns
    -------
    dict of KaplanMeierFitter objects, log-rank test result, matplotlib Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    groups = clinical_df[group_col].dropna().unique()
    colors = colors or ["#2E75B6", "#C0392B", "#27AE60", "#8E44AD"]
    kmf_dict = {}

    for i, group in enumerate(sorted(groups)):
        subset = clinical_df[clinical_df[group_col] == group]
        label = group_labels.get(group, str(group)) if group_labels else str(group)

        kmf = KaplanMeierFitter()
        kmf.fit(
            durations=subset[duration_col],
            event_observed=subset[event_col],
            label=f"{label} (n={len(subset)})"
        )
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[i % len(colors)])
        kmf_dict[group] = kmf

    # Log-rank test (binary groups only)
    if len(groups) == 2:
        g0, g1 = sorted(groups)
        s0 = clinical_df[clinical_df[group_col] == g0]
        s1 = clinical_df[clinical_df[group_col] == g1]
        lr = logrank_test(
            s0[duration_col], s1[duration_col],
            event_observed_A=s0[event_col],
            event_observed_B=s1[event_col]
        )
        p = lr.p_value
        p_text = f"p = {p:.4f}" if p >= 0.0001 else "p < 0.0001"
        ax.text(0.65, 0.85, f"Log-rank {p_text}", transform=ax.transAxes,
                fontsize=11, bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    else:
        lr = None

    ax.set_xlabel("Time (months)", fontsize=12)
    ax.set_ylabel("Survival probability", fontsize=12)
    ax.set_title(title or f"Survival by {group_col}", fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return kmf_dict, lr, ax


def run_cox_regression(clinical_df: pd.DataFrame,
                        covariates: list,
                        duration_col: str = "OS_MONTHS",
                        event_col: str = "OS_STATUS") -> CoxPHFitter:
    """
    Fit a Cox Proportional Hazards model.

    Parameters
    ----------
    covariates : list of column names to include as predictors
    duration_col, event_col : survival outcome columns

    Returns
    -------
    Fitted CoxPHFitter object (call .print_summary() for results)
    """
    cols = [duration_col, event_col] + covariates
    df = clinical_df[cols].dropna()

    cph = CoxPHFitter()
    cph.fit(df, duration_col=duration_col, event_col=event_col)
    return cph


def plot_cox_forest(cph: CoxPHFitter, ax=None, title: str = "Cox PH Model — Hazard Ratios"):
    """
    Plot a forest plot of hazard ratios from a fitted Cox model.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, max(4, len(cph.params_) * 0.6)))

    summary = cph.summary
    coefs = summary["exp(coef)"]
    lower = summary["exp(coef) lower 95%"]
    upper = summary["exp(coef) upper 95%"]
    p_vals = summary["p"]

    y_pos = range(len(coefs))
    ax.errorbar(
        coefs.values, list(y_pos),
        xerr=[coefs.values - lower.values, upper.values - coefs.values],
        fmt="o", color="#2E75B6", ecolor="#95A5A6", capsize=4, markersize=7
    )
    ax.axvline(1.0, linestyle="--", color="grey", alpha=0.7, label="HR = 1 (no effect)")

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels([f"{v}  (p={p:.3f})" for v, p in zip(coefs.index, p_vals)], fontsize=10)
    ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax
