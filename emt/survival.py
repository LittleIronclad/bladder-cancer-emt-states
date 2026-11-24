"""Survival utilities: LASSO Cox with honest (nested) validation, C-index, log-rank tests, KM plots."""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from sklearn.model_selection import KFold, RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from statsmodels.stats.multitest import multipletests

from . import config as C


# ── C-index ───────────────────────────────────────────────────────────────────
def cindex(time, event, risk) -> float:
    """Harrell's C (higher risk ↔ shorter survival)."""
    return concordance_index_censored(np.asarray(event, bool), np.asarray(time, float),
                                      np.asarray(risk, float))[0]


def bootstrap_cindex(time, event, risk, n=C.N_BOOTSTRAP, seed=C.SEED):
    rng = np.random.default_rng(seed)
    time, event, risk = map(np.asarray, (time, event, risk))
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(time), len(time))
        if event[i].sum() > 1:
            vals.append(cindex(time[i], event[i], risk[i]))
    return np.percentile(vals, [2.5, 97.5])


# ── LASSO Cox ─────────────────────────────────────────────────────────────────
def _y(time, event):
    return np.array(list(zip(np.asarray(event, bool), np.asarray(time, float))),
                    dtype=[("event", bool), ("time", float)])


def _fit(X, y, alpha):
    return CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[alpha], max_iter=1000).fit(X, y)


def select_alpha(X, y, alphas=C.LASSO_ALPHAS, folds=C.LASSO_CV_FOLDS, seed=C.SEED):
    """Inner K-fold CV over the alpha grid; returns (best alpha, mean CV C-index per alpha)."""
    scores = np.full((folds, len(alphas)), np.nan)
    for f, (tr, va) in enumerate(KFold(folds, shuffle=True, random_state=seed).split(X)):
        sc = StandardScaler().fit(X[tr])
        Xtr, Xva = sc.transform(X[tr]), sc.transform(X[va])
        for j, a in enumerate(alphas):
            try:
                m = _fit(Xtr, y[tr], a)
                scores[f, j] = cindex(y[va]["time"], y[va]["event"], m.predict(Xva))
            except Exception:
                pass
    mean = np.nanmean(scores, axis=0)
    return alphas[int(np.nanargmax(mean))], mean


class LassoCox:
    """Standardise → L1-penalised Cox; alpha chosen by inner CV (max mean C-index)."""

    def fit(self, X: pd.DataFrame, time, event, seed=C.SEED):
        self.genes_ = list(X.columns)
        y = _y(time, event)
        self.alpha_, self.cv_curve_ = select_alpha(X.values, y, seed=seed)
        self.scaler_ = StandardScaler().fit(X.values)
        self.model_ = _fit(self.scaler_.transform(X.values), y, self.alpha_)
        self.coef_ = pd.Series(self.model_.coef_[:, 0], index=self.genes_)
        return self

    def risk(self, X: pd.DataFrame) -> np.ndarray:
        return self.model_.predict(self.scaler_.transform(X[self.genes_].values))


def nested_cv(X: pd.DataFrame, time, event, folds=5, repeats=C.NESTED_CV_REPEATS, seed=0):
    """Outer repeated stratified K-fold; alpha re-selected inside every outer training fold.

    Returns (per-fold C-indices, out-of-fold risk of the first repeat). The out-of-fold risk is
    a cross-fitted score: each patient's risk comes from a model that never saw that patient.
    """
    time, event = np.asarray(time, float), np.asarray(event, int)
    y = _y(time, event)
    outer, oof = [], np.full(len(X), np.nan)
    rskf = RepeatedStratifiedKFold(n_splits=folds, n_repeats=repeats, random_state=seed)
    for k, (tr, te) in enumerate(rskf.split(X, event)):
        a, _ = select_alpha(X.values[tr], y[tr], seed=k)
        sc = StandardScaler().fit(X.values[tr])
        m = _fit(sc.transform(X.values[tr]), y[tr], a)
        r = m.predict(sc.transform(X.values[te]))
        outer.append(cindex(time[te], event[te], r))
        if k < folds:                                   # first repeat → one OOF prediction per patient
            oof[te] = r
    return np.array(outer), oof


def transfer_risk(model: LassoCox, expr: pd.DataFrame) -> tuple[np.ndarray, list]:
    """Apply TCGA coefficients to a new cohort without refitting.

    Genes are z-scored within the new cohort (platform-independent scale); genes missing on
    the platform are dropped from the linear predictor.
    """
    nz = model.coef_[model.coef_ != 0]
    genes = [g for g in nz.index if g in expr.columns]
    z = (expr[genes] - expr[genes].mean()) / (expr[genes].std(ddof=0) + 1e-8)
    return z.values @ nz.loc[genes].values, genes


# ── Group comparisons ─────────────────────────────────────────────────────────
def global_logrank(df, group_col, time="time", event="event"):
    r = multivariate_logrank_test(df[time], df[group_col], df[event])
    return float(r.p_value)


def pairwise_logrank(df, group_col, order, time="time", event="event"):
    rows = []
    for a, b in combinations(order, 2):
        A, B = df[df[group_col] == a], df[df[group_col] == b]
        if len(A) < 5 or len(B) < 5:
            continue
        r = logrank_test(A[time], B[time], event_observed_A=A[event], event_observed_B=B[event])
        rows.append({"group_1": a, "group_2": b, "p": r.p_value})
    out = pd.DataFrame(rows)
    if len(out):
        out["p_BH"] = multipletests(out["p"], method="fdr_bh")[1]
    return out


def median_survival_table(df, group_col, order, time="time", event="event"):
    from lifelines.utils import median_survival_times
    rows = []
    for g in order:
        s = df[df[group_col] == g]
        if len(s) < 5:
            continue
        km = KaplanMeierFitter().fit(s[time], s[event])
        ci = median_survival_times(km.confidence_interval_).values[0]
        rows.append({"group": g, "n": len(s), "events": int(s[event].sum()),
                     "median_years": km.median_survival_time_, "ci_low": ci[0], "ci_high": ci[1]})
    return pd.DataFrame(rows)


# ── Plotting ──────────────────────────────────────────────────────────────────
SHORT = {"Epithelial + Partial-E": "E + PE", "Partial-M + Mesenchymal": "PM + M", "Epithelial": "E",
         "Partial-E": "PE", "Partial-M": "PM", "Mesenchymal": "M"}


def km_plot(ax, df, group_col, order, colors, time="time", event="event", at_risk=True,
            linestyles=None, xmax=C.MAX_FOLLOWUP_YEARS):
    fits = []
    for g in order:
        s = df[df[group_col] == g]
        if len(s) < 5:
            continue
        km = KaplanMeierFitter(label=SHORT.get(g, g)).fit(s[time], s[event])
        sf, ci = km.survival_function_.iloc[:, 0], km.confidence_interval_
        ls = (linestyles or {}).get(g, "-")
        ax.step(sf.index, sf.values, where="post", color=colors[g], lw=2.2, ls=ls, label=f"{g} (n={len(s)})")
        ax.fill_between(ci.index, ci.iloc[:, 0], ci.iloc[:, 1], step="post", color=colors[g], alpha=0.10)
        fits.append(km)
    ax.set(xlim=(0, xmax), ylim=(0, 1.03), xlabel="Years", ylabel="Survival probability")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    if at_risk and fits:
        from lifelines.plotting import add_at_risk_counts
        add_at_risk_counts(*fits, ax=ax, rows_to_show=["At risk"], fontsize=8)
    return fits


def p_fmt(p):
    return "p < 0.001" if p < 1e-3 else f"p = {p:.3f}"
