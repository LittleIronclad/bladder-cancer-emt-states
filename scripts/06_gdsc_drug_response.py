"""Step 6 — EMT score of bladder cancer cell lines vs drug sensitivity (GDSC2).

Cell-line RNA-seq (Cell Model Passports, RSEM TPM) → log2(TPM+1) → anchor-marker EMT score
(z-scored across bladder lines). Spearman correlation with LN_IC50 per drug, Benjamini–Hochberg
FDR across drugs, and Kruskal–Wallis across EMT tertiles.

Input (data/gdsc/): GDSC2_fitted_dose_response_*.xlsx, model_list_*.csv, rnaseq_merged_*.zip
Output: results/gdsc_*.csv, figures/06_*.png
"""
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

from emt import banner, config as C, savefig  # noqa: E402
from emt.scoring import emt_score  # noqa: E402

banner("Step 6 · GDSC2: EMT status vs drug sensitivity in bladder cancer cell lines")


def newest(pattern):
    files = sorted(C.GDSC_DIR.glob(pattern))
    if not files:
        raise FileNotFoundError(f"data/gdsc/{pattern} not found — see data/README.md")
    return files[-1]


models = pd.read_csv(newest("model_list_*.csv"))
blca = models[models["tissue"].str.contains("bladder", case=False, na=False)]
id_to_name = blca.set_index("model_id")["model_name"].to_dict()
print(f"Bladder cancer models in Cell Model Passports: {len(blca)}")

# RNA-seq: rows 0–3 are metadata (model_id, model_name, data_source, column names); genes from row 4.
zpath = newest("rnaseq_merged_*.zip")
with zipfile.ZipFile(zpath) as z:
    member = next(n for n in z.namelist() if "rsem_tpm" in n and n.endswith(".csv"))
    with z.open(member) as f:
        head = pd.read_csv(f, nrows=4, header=None, low_memory=False)
    model_ids = head.iloc[0, 3:].astype(str).values
    keep = [i + 3 for i, m in enumerate(model_ids) if m in id_to_name]
    with z.open(member) as f:
        raw = pd.read_csv(f, header=None, skiprows=4, usecols=[0] + keep, low_memory=False)
raw.columns = ["gene"] + [model_ids[i - 3] for i in keep]
raw = raw.dropna(subset=["gene"]).set_index("gene").apply(pd.to_numeric, errors="coerce")
raw = raw.loc[raw.var(axis=1).sort_values(ascending=False).index]
raw = raw[~raw.index.duplicated(keep="first")]                       # duplicated symbol → most variable
expr = np.log2(raw.T.clip(lower=0) + 1)
expr.index = [id_to_name[m] for m in expr.index]
print(f"RNA-seq: {expr.shape[0]} bladder cell lines × {expr.shape[1]} genes")

score = emt_score(expr)
tert = pd.qcut(score, 3, labels=["Epithelial", "Intermediate", "Mesenchymal"])
lines = pd.DataFrame({"EMT_score": score, "EMT_tertile": tert})
lines.index.name = "cell_line"
lines.sort_values("EMT_score").to_csv(C.RESULTS_DIR / "gdsc_cell_line_emt_scores.csv")

ic50 = pd.read_excel(newest("GDSC2_fitted_dose_response_*.xlsx"))
ic50 = ic50[ic50["CELL_LINE_NAME"].isin(lines.index) & ic50["DRUG_NAME"].isin(C.GDSC_DRUGS)]
ic50 = ic50.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False)["LN_IC50"].mean()   # one value per line × drug
d = ic50.merge(lines, left_on="CELL_LINE_NAME", right_index=True)

rows = []
for drug, g in d.groupby("DRUG_NAME"):
    if len(g) < C.GDSC_MIN_LINES:
        continue
    rho, p = stats.spearmanr(g.EMT_score, g.LN_IC50)
    kw = stats.kruskal(*[x.LN_IC50 for _, x in g.groupby("EMT_tertile", observed=True) if len(x) > 1]).pvalue
    rows.append({"drug": drug, "n_lines": len(g), "spearman_rho": rho, "p": p, "kruskal_p_tertiles": kw})
res = pd.DataFrame(rows)
res["FDR_BH"] = multipletests(res.p, method="fdr_bh")[1]
res = res.sort_values("p")
res.to_csv(C.RESULTS_DIR / "gdsc_emt_ic50_correlations.csv", index=False)
print(res.round(4).to_string(index=False))
print(f"\nPositive ρ = more mesenchymal lines are more resistant (higher IC50). "
      f"{(res.FDR_BH < 0.05).sum()} drugs at FDR < 0.05.")

fig, ax = plt.subplots(figsize=(7, 3.8))
r = res.sort_values("spearman_rho")
bars = ax.barh(r.drug, r.spearman_rho, color=np.where(r.spearman_rho > 0, "#B2182B", "#2166AC"))
for b, (_, row) in zip(bars, r.iterrows()):
    ax.text(row.spearman_rho + (0.02 if row.spearman_rho >= 0 else -0.02), b.get_y() + b.get_height() / 2,
            f"p={row.p:.3f}, FDR={row.FDR_BH:.2f}", va="center", ha="left" if row.spearman_rho >= 0 else "right",
            fontsize=7)
ax.axvline(0, color="black", lw=0.8)
ax.set(xlim=(-1, 1), xlabel="Spearman ρ (EMT score vs LN IC50)",
       title=f"GDSC2 bladder cancer lines (n = {res.n_lines.min()}–{res.n_lines.max()} per drug)")
savefig(fig, "06_gdsc_correlations.png")

top = res.iloc[0]
g = d[d.DRUG_NAME == top.drug]
fig, ax = plt.subplots(figsize=(4.8, 4))
ax.scatter(g.EMT_score, g.LN_IC50, color="#444444")
for _, row in g.iterrows():
    ax.annotate(row.CELL_LINE_NAME, (row.EMT_score, row.LN_IC50), fontsize=6, alpha=0.7)
ax.set(xlabel="EMT score (0 = E, 1 = M)", ylabel="LN IC50",
       title=f"{top.drug}: ρ = {top.spearman_rho:.2f}, p = {top.p:.3f}, FDR = {top.FDR_BH:.2f}")
savefig(fig, "06_gdsc_top_drug.png")
print("Saved: results/gdsc_*.csv")
