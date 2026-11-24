"""Step 3 — Prognostic value of EMT states and a LASSO Cox model on the 42-gene panel (TCGA-BLCA).

* Kaplan–Meier / log-rank for the five states and three functional groups
  (Epithelial + Partial-E | Hybrid | Partial-M + Mesenchymal), with BH-adjusted pairwise tests.
* LASSO Cox: alpha by 5-fold CV; performance reported both as the apparent (training) C-index
  and as an honest nested-CV C-index (5 × 5 folds, alpha re-tuned inside every training fold).
* Risk groups and the multivariable Cox model use the *cross-fitted* (out-of-fold) risk score,
  so the hazard ratio is not inflated by evaluating the model on its own training data.
* Landmark analysis at 2 years.

Output: work/tcga_survival.csv, results/lasso_cox_model.json, results/lasso_cox_coefficients.csv,
        results/survival_summary.json, results/tcga_*.csv, figures/03_*.png
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines import CoxPHFitter  # noqa: E402
from lifelines.statistics import logrank_test  # noqa: E402

from emt import banner, config as C, savefig  # noqa: E402
from emt import survival as S  # noqa: E402

banner("Step 3 · Survival analysis and LASSO Cox model (TCGA-BLCA)")
panel = pd.read_csv(C.WORK_DIR / "tcga_panel.csv", index_col=0)
clin = pd.read_csv(C.WORK_DIR / "tcga_clinical.csv", index_col=0)
states = pd.read_csv(C.WORK_DIR / "tcga_emt_states.csv", index_col=0)
df = clin.join(states[["EMT_state", "EMT_score"]])
df["group"] = df["EMT_state"].map(C.broad_group)
order5 = [s for s in C.STATE_NAMES[5] if s in set(df.EMT_state)]
print(f"{len(df)} patients, {int(df.event.sum())} deaths within 10 years")

# ── EMT states and survival ───────────────────────────────────────────────────
p5 = S.global_logrank(df, "EMT_state")
p3 = S.global_logrank(df, "group")
pw5 = S.pairwise_logrank(df, "EMT_state", order5)
pw3 = S.pairwise_logrank(df, "group", C.BROAD_ORDER)
med5 = S.median_survival_table(df, "EMT_state", order5)
med3 = S.median_survival_table(df, "group", C.BROAD_ORDER)
print(f"\nFive states: global log-rank {S.p_fmt(p5)}")
print(med5.round(2).to_string(index=False))
print(f"\nThree groups: global log-rank {S.p_fmt(p3)}")
print(med3.round(2).to_string(index=False))
print(pw3.round(4).to_string(index=False))
pd.concat([med5.assign(level="state"), med3.assign(level="group")]).to_csv(
    C.RESULTS_DIR / "tcga_median_survival.csv", index=False)
pd.concat([pw5.assign(level="state"), pw3.assign(level="group")]).to_csv(
    C.RESULTS_DIR / "tcga_pairwise_logrank.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"wspace": 0.3})
S.km_plot(axes[0], df, "EMT_state", order5, C.STATE_COLORS, linestyles={"Hybrid": "--"})
axes[0].set_title(f"Five EMT states — log-rank {S.p_fmt(p5)}")
S.km_plot(axes[1], df, "group", C.BROAD_ORDER, C.BROAD_COLORS, linestyles={"Hybrid": "--"})
axes[1].set_title(f"Three functional groups — log-rank {S.p_fmt(p3)}")
fig.suptitle("TCGA-BLCA overall survival by EMT state", y=1.02)
savefig(fig, "03_km_emt_states.png")

# ── LASSO Cox ─────────────────────────────────────────────────────────────────
X = panel.loc[df.index]
model = S.LassoCox().fit(X, df.time, df.event)
coef = model.coef_[model.coef_ != 0].sort_values(ascending=False)
risk_train = model.risk(X)
c_app = S.cindex(df.time, df.event, risk_train)
print(f"\nLASSO Cox: alpha = {model.alpha_:.4g}, {len(coef)}/{X.shape[1]} genes with non-zero coefficients")
print(f"  apparent C-index (train = test): {c_app:.3f}")

print("  nested cross-validation (5 × 5 folds, alpha re-tuned inside each fold) …")
outer, oof = S.nested_cv(X, df.time, df.event)
c_oof = S.cindex(df.time, df.event, oof)
print(f"  nested-CV C-index: {outer.mean():.3f} ± {outer.std():.3f}  (out-of-fold, single repeat: {c_oof:.3f})")

coef_df = pd.DataFrame({"coef": coef, "HR_per_SD": np.exp(coef)})
coef_df.index.name = "gene"
coef_df.to_csv(C.RESULTS_DIR / "lasso_cox_coefficients.csv")

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), gridspec_kw={"width_ratios": [1, 1.2]})
alphas = np.array(C.LASSO_ALPHAS)
axes[0].semilogx(alphas, model.cv_curve_, "o-", ms=3)
axes[0].axvline(model.alpha_, color="crimson", ls="--", lw=1)
axes[0].set(xlabel="alpha (L1 penalty)", ylabel="mean CV C-index", title="LASSO Cox: penalty selection")
top = coef.reindex(coef.abs().sort_values(ascending=False).index[:20]).sort_values()
axes[1].barh(top.index, top.values, color=["#B2182B" if v > 0 else "#2166AC" for v in top.values])
axes[1].axvline(0, color="black", lw=0.8)
axes[1].set(xlabel="coefficient (log HR per SD)",
            title=f"Top-20 genes · apparent C = {c_app:.2f}, nested-CV C = {outer.mean():.2f}")
savefig(fig, "03_lasso_cox.png")

# ── Cross-fitted risk groups and multivariable Cox ────────────────────────────
df["risk_oof"] = oof
df["risk_group"] = np.where(oof >= np.median(oof), "High risk", "Low risk")
hi = df.risk_group == "High risk"
lr = logrank_test(df.time[hi], df.time[~hi], df.event[hi], df.event[~hi])
fig, ax = plt.subplots(figsize=(6.8, 5.2))
S.km_plot(ax, df, "risk_group", ["Low risk", "High risk"], {"Low risk": "#2166AC", "High risk": "#B2182B"})
ax.set_title(f"Cross-fitted LASSO risk score (median split)\nlog-rank {S.p_fmt(lr.p_value)}, C = {c_oof:.3f}")
savefig(fig, "03_km_risk_groups.png")

stage_map = {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4, "Stage IVA": 4, "Stage IVB": 4}
mv = df.assign(stage_num=df.stage.map(stage_map))[["time", "event", "risk_oof", "age", "stage_num"]].dropna()
mv_z = mv.copy()
mv_z[["risk_oof", "age", "stage_num"]] = (mv[["risk_oof", "age", "stage_num"]] - mv[["risk_oof", "age", "stage_num"]].mean()) \
    / mv[["risk_oof", "age", "stage_num"]].std()
cph = CoxPHFitter().fit(mv_z, "time", "event")
summ = cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
summ.columns = ["HR_per_SD", "ci_low", "ci_high", "p"]
summ.index = ["LASSO risk (cross-fitted)", "Age", "Stage (I–IV)"]
print(f"\nMultivariable Cox (n = {len(mv)}, {int(mv.event.sum())} events), HR per SD:")
print(summ.round(4).to_string())
summ.to_csv(C.RESULTS_DIR / "tcga_multivariable_cox.csv")

fig, ax = plt.subplots(figsize=(6.5, 2.6))
y = np.arange(len(summ))[::-1]
ax.errorbar(summ.HR_per_SD, y, xerr=[summ.HR_per_SD - summ.ci_low, summ.ci_high - summ.HR_per_SD],
            fmt="s", color="black", capsize=4)
ax.axvline(1, color="grey", ls="--")
ax.set(yticks=y, yticklabels=summ.index, xlabel="Hazard ratio per SD (95% CI)",
       title=f"Multivariable Cox, TCGA-BLCA (n = {len(mv)})")
for yi, (_, r) in zip(y, summ.iterrows()):
    ax.text(r.ci_high + 0.05, yi, f"{r.HR_per_SD:.2f} ({S.p_fmt(r.p)})", va="center", fontsize=8)
savefig(fig, "03_multivariable_cox.png")

# ── Landmark analysis (2 years) ───────────────────────────────────────────────
lm = df[df.time > C.LANDMARK_YEARS].assign(time=lambda d: d.time - C.LANDMARK_YEARS)
p_lm = S.global_logrank(lm, "group")
pw_lm = S.pairwise_logrank(lm, "group", C.BROAD_ORDER)
hyb = lm[lm.group.isin(["Hybrid", "Partial-M + Mesenchymal"])].assign(hybrid=lambda d: (d.group == "Hybrid").astype(int))
hr_lm = CoxPHFitter().fit(hyb[["time", "event", "hybrid"]], "time", "event").summary.loc["hybrid"]
print(f"\nLandmark 2 y: n = {len(lm)}, events = {int(lm.event.sum())}, global log-rank {S.p_fmt(p_lm)}")
print(pw_lm.round(4).to_string(index=False))
print(f"  HR Hybrid vs Partial-M + Mesenchymal = {hr_lm['exp(coef)']:.2f} "
      f"({hr_lm['exp(coef) lower 95%']:.2f}–{hr_lm['exp(coef) upper 95%']:.2f}), {S.p_fmt(hr_lm['p'])}")
pw_lm.to_csv(C.RESULTS_DIR / "tcga_landmark_pairwise_logrank.csv", index=False)
fig, ax = plt.subplots(figsize=(6.8, 5.2))
S.km_plot(ax, lm, "group", C.BROAD_ORDER, C.BROAD_COLORS, linestyles={"Hybrid": "--"},
          xmax=C.MAX_FOLLOWUP_YEARS - C.LANDMARK_YEARS)
ax.set(xlabel=f"Years after the {C.LANDMARK_YEARS:.0f}-year landmark",
       title=f"Landmark analysis: patients alive at 2 years (n = {len(lm)})\nlog-rank {S.p_fmt(p_lm)}")
savefig(fig, "03_landmark_2y.png")

# ── Save ──────────────────────────────────────────────────────────────────────
df.to_csv(C.WORK_DIR / "tcga_survival.csv")
model_out = {
    "genes": list(coef.index), "coef": [float(v) for v in coef.values],
    "alpha": float(model.alpha_), "tcga_median_risk": float(np.median(risk_train)),
    "note": "Risk = sum(coef * z), z = per-gene z-score within the cohort (ddof = 0).",
}
(C.RESULTS_DIR / "lasso_cox_model.json").write_text(json.dumps(model_out, indent=2))
ci_app = S.bootstrap_cindex(df.time, df.event, risk_train)
summary = {
    "n": int(len(df)), "events": int(df.event.sum()),
    "logrank_p_5_states": p5, "logrank_p_3_groups": p3, "logrank_p_landmark_2y": p_lm,
    "landmark_n": int(len(lm)), "landmark_events": int(lm.event.sum()),
    "landmark_hr_hybrid_vs_pm_m": [round(float(hr_lm["exp(coef)"]), 2), round(float(hr_lm["exp(coef) lower 95%"]), 2),
                                   round(float(hr_lm["exp(coef) upper 95%"]), 2), float(hr_lm["p"])],
    "median_os_years": med3.set_index("group")["median_years"].round(2).to_dict(),
    "lasso_alpha": float(model.alpha_), "lasso_n_genes": int(len(coef)),
    "c_index_apparent": round(c_app, 3), "c_index_apparent_bootstrap_ci": [round(v, 3) for v in ci_app],
    "c_index_nested_cv_mean": round(float(outer.mean()), 3), "c_index_nested_cv_sd": round(float(outer.std()), 3),
    "c_index_out_of_fold": round(c_oof, 3), "risk_groups_logrank_p_oof": float(lr.p_value),
    "multivariable_hr_risk_per_sd": round(float(summ.iloc[0].HR_per_SD), 2),
    "multivariable_hr_ci": [round(float(summ.iloc[0].ci_low), 2), round(float(summ.iloc[0].ci_high), 2)],
    "multivariable_p": float(summ.iloc[0].p),
}
(C.RESULTS_DIR / "survival_summary.json").write_text(json.dumps(summary, indent=2))
print("\nSaved: results/survival_summary.json, lasso_cox_model.json, lasso_cox_coefficients.csv")
