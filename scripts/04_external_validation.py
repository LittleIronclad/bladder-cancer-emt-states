"""Step 4 — External validation in four independent GEO cohorts (no refitting).

Cohorts: GSE32894 (Sjödahl 2012), GSE48075 (Choi 2014, cT2+), GSE31684 (Riester 2012), GSE13507 (Kim 2010).
Meta-GEO = GSE13507 + GSE31684 + GSE48075 (z-scored within each cohort before pooling).

For every cohort:
  A. LASSO risk score with the frozen TCGA coefficients → C-index (bootstrap 95% CI) and
     log-rank for high vs low risk (TCGA median threshold).
  B. EMT states transferred with a nearest-centroid classifier trained on the TCGA NMF states
     → log-rank across the three functional groups.

Output: results/validation_summary.csv, figures/04_*.png
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines.statistics import logrank_test  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

from emt import banner, config as C, savefig  # noqa: E402
from emt import survival as S  # noqa: E402
from emt.data import geo_cohort  # noqa: E402
from emt.scoring import CentroidStateClassifier  # noqa: E402

banner("Step 4 · External validation (GEO)")
panel = pd.read_csv(C.WORK_DIR / "tcga_panel.csv", index_col=0)
tcga = pd.read_csv(C.WORK_DIR / "tcga_survival.csv", index_col=0)
model = json.loads((C.RESULTS_DIR / "lasso_cox_model.json").read_text())
coef = pd.Series(model["coef"], index=model["genes"])
threshold = model["tcga_median_risk"]

# Nearest-centroid transfer of EMT states: check agreement with NMF labels inside TCGA (5-fold CV)
agree = []
for tr, te in StratifiedKFold(5, shuffle=True, random_state=C.SEED).split(panel, tcga.EMT_state):
    clf = CentroidStateClassifier().fit(panel.iloc[tr], tcga.EMT_state.iloc[tr])
    agree.append((clf.predict(panel.iloc[te]) == tcga.EMT_state.iloc[te]).mean())
clf = CentroidStateClassifier().fit(panel, tcga.EMT_state)
print(f"Centroid classifier reproduces TCGA NMF states in held-out tumours: {np.mean(agree):.1%} "
      f"(chance ≈ {1 / tcga.EMT_state.nunique():.0%})")


def risk_score(expr):
    genes = [g for g in coef.index if g in expr.columns]
    z = (expr[genes] - expr[genes].mean()) / (expr[genes].std(ddof=0) + 1e-8)
    return z.values @ coef.loc[genes].values, len(genes)


def evaluate(name, expr, surv):
    d = surv.copy()
    d["risk"], n_genes = risk_score(expr)
    d["risk_group"] = np.where(d.risk >= threshold, "High risk", "Low risk")
    d["EMT_state"] = clf.predict(expr)
    d["group"] = d.EMT_state.map(C.broad_group)
    c = S.cindex(d.time, d.event, d.risk)
    lo, hi = S.bootstrap_cindex(d.time, d.event, d.risk)
    h = d.risk_group == "High risk"
    p_risk = logrank_test(d.time[h], d.time[~h], d.event[h], d.event[~h]).p_value
    p_grp = S.global_logrank(d, "group")
    pw = S.pairwise_logrank(d, "group", C.BROAD_ORDER).set_index(["group_1", "group_2"])["p"]
    med = S.median_survival_table(d, "group", C.BROAD_ORDER).set_index("group")["median_years"]
    row = {"cohort": name, "endpoint": "/".join(sorted(set(d.endpoint))), "n": len(d), "events": int(d.event.sum()),
           "lasso_genes": f"{n_genes}/{len(coef)}", "c_index": c, "ci_low": lo, "ci_high": hi,
           "logrank_p_risk": p_risk, "logrank_p_emt_groups": p_grp,
           "p_hybrid_vs_epithelial": pw.get(("Epithelial + Partial-E", "Hybrid"), np.nan),
           "p_hybrid_vs_mesenchymal": pw.get(("Hybrid", "Partial-M + Mesenchymal"), np.nan),
           **{f"median_os_{g}": med.get(g, np.nan) for g in C.BROAD_ORDER},
           **{f"n_{g}": int((d.group == g).sum()) for g in C.BROAD_ORDER}}
    print(f"\n{name}: n={len(d)}, events={row['events']} ({row['endpoint']}), genes {row['lasso_genes']}")
    print(f"  C-index {c:.3f} ({lo:.3f}–{hi:.3f}); risk groups {S.p_fmt(p_risk)}; EMT groups {S.p_fmt(p_grp)}")
    print("  states:", d.EMT_state.value_counts().to_dict())
    return d, row


cohorts, rows, frames = {}, [], {}
for gse in ["GSE32894", "GSE48075", "GSE31684", "GSE13507"]:
    expr, surv = geo_cohort(gse)
    cohorts[gse] = (expr, surv)
    frames[gse], row = evaluate(gse, expr, surv)
    rows.append(row)

# Meta-GEO: genes z-scored within each cohort, then pooled
meta_ids = ["GSE13507", "GSE31684", "GSE48075"]
genes = sorted(set.intersection(*[set(cohorts[g][0].columns) for g in meta_ids]))
z = pd.concat([(cohorts[g][0][genes] - cohorts[g][0][genes].mean()) / cohorts[g][0][genes].std() for g in meta_ids])
s = pd.concat([cohorts[g][1].assign(cohort=g) for g in meta_ids])
frames["Meta-GEO"], row = evaluate("Meta-GEO (GSE13507 + GSE31684 + GSE48075)", z, s)
rows.append(row)

# TCGA reference rows (apparent and nested-CV)
tsum = json.loads((C.RESULTS_DIR / "survival_summary.json").read_text())
ref = pd.DataFrame([
    {"cohort": "TCGA-BLCA (training, apparent)", "endpoint": "OS", "n": tsum["n"], "events": tsum["events"],
     "c_index": tsum["c_index_apparent"], "ci_low": tsum["c_index_apparent_bootstrap_ci"][0],
     "ci_high": tsum["c_index_apparent_bootstrap_ci"][1], "logrank_p_emt_groups": tsum["logrank_p_3_groups"]},
    {"cohort": "TCGA-BLCA (nested CV)", "endpoint": "OS", "n": tsum["n"], "events": tsum["events"],
     "c_index": tsum["c_index_nested_cv_mean"], "ci_low": tsum["c_index_nested_cv_mean"] - 1.96 * tsum["c_index_nested_cv_sd"],
     "ci_high": tsum["c_index_nested_cv_mean"] + 1.96 * tsum["c_index_nested_cv_sd"]},
])
table = pd.concat([ref, pd.DataFrame(rows)], ignore_index=True)
table.to_csv(C.RESULTS_DIR / "validation_summary.csv", index=False)
print("\n" + table[["cohort", "endpoint", "n", "events", "c_index", "ci_low", "ci_high",
                    "logrank_p_risk", "logrank_p_emt_groups"]].round(4).to_string(index=False))

# ── Figures ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7.5, 3.6))
y = np.arange(len(table))[::-1]
colors = ["#888888", "#444444"] + ["#2166AC"] * 4 + ["#B2182B"]
for yi, (_, r), col in zip(y, table.iterrows(), colors):
    ax.errorbar(r.c_index, yi, xerr=[[r.c_index - r.ci_low], [r.ci_high - r.c_index]], fmt="o", color=col, capsize=4)
    ax.text(r.ci_high + 0.01, yi, f"{r.c_index:.3f}  (n={int(r.n)}, {int(r.events)} ev.)", va="center", fontsize=8)
ax.axvline(0.5, color="grey", ls=":")
ax.set(yticks=y, yticklabels=[c.replace(" (GSE13507 + GSE31684 + GSE48075)", "") for c in table.cohort],
       xlim=(0.35, 1.0), xlabel="Harrell's C-index (95% CI)",
       title="LASSO Cox risk score: internal and external performance")
savefig(fig, "04_cindex_all_cohorts.png")

meta = frames["Meta-GEO"]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"wspace": 0.3})
S.km_plot(axes[0], meta, "risk_group", ["Low risk", "High risk"], {"Low risk": "#2166AC", "High risk": "#B2182B"})
r = table.iloc[-1]
axes[0].set_title(f"Meta-GEO: frozen TCGA risk score\nC = {r.c_index:.3f}, log-rank {S.p_fmt(r.logrank_p_risk)}")
S.km_plot(axes[1], meta, "group", C.BROAD_ORDER, C.BROAD_COLORS, linestyles={"Hybrid": "--"})
axes[1].set_title(f"Meta-GEO: transferred EMT states\nlog-rank {S.p_fmt(r.logrank_p_emt_groups)}")
savefig(fig, "04_km_meta_geo.png")
print("\nSaved: results/validation_summary.csv")
