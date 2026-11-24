"""EMT scores and transfer of TCGA-derived EMT states to new cohorts."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C


def zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Per-gene z-score across samples (within one cohort)."""
    return (df - df.mean()) / (df.std(ddof=1) + 1e-8)


def emt_score(expr: pd.DataFrame, e_genes=C.E_ANCHOR, m_genes=C.M_ANCHOR, normalise=True) -> pd.Series:
    """mean z(M anchors) − mean z(E anchors); optionally min–max scaled to [0, 1] (0 = E, 1 = M)."""
    e = [g for g in e_genes if g in expr.columns]
    m = [g for g in m_genes if g in expr.columns]
    z = zscore(expr[e + m])
    s = z[m].mean(axis=1) - z[e].mean(axis=1)
    if normalise:
        s = (s - s.min()) / (s.max() - s.min() + 1e-12)
    return s.rename("EMT_score")


class CentroidStateClassifier:
    """Nearest-centroid transfer of the TCGA NMF-derived EMT states.

    Centroids are the mean per-gene z-scores of each state on the EMT panel in TCGA. A new
    cohort is z-scored within itself and every sample is assigned to the centroid with the
    highest Pearson correlation, using only the genes measured in that cohort. This keeps the
    transferred labels tied to the multi-dimensional NMF states rather than to quantiles of a
    one-dimensional score.
    """

    def fit(self, expr_panel: pd.DataFrame, states: pd.Series):
        z = zscore(expr_panel)
        self.centroids_ = z.groupby(states.loc[z.index].values).mean()
        return self

    def predict(self, expr_panel: pd.DataFrame) -> pd.Series:
        genes = [g for g in self.centroids_.columns if g in expr_panel.columns]
        z = zscore(expr_panel[genes]).values
        c = self.centroids_[genes].values
        zc = (z - z.mean(1, keepdims=True)) / (z.std(1, keepdims=True) + 1e-12)
        cc = (c - c.mean(1, keepdims=True)) / (c.std(1, keepdims=True) + 1e-12)
        corr = zc @ cc.T / len(genes)
        return pd.Series(self.centroids_.index[corr.argmax(1)], index=expr_panel.index, name="EMT_state")
