"""Step 5 — Pathway enrichment (GSEA prerank) of each EMT state versus all other tumours.

Ranking metric per gene: sign(t) · −log10(p) from Welch's t-test (state vs rest) on the full
TCGA transcriptome. Gene set libraries are the complete Hallmark (MSigDB 2020), KEGG 2021 and
GO Biological Process 2023 collections from Enrichr; they are downloaded once by gseapy and
cached as GMT files in data/genesets/ (place GMT files there to run offline).

Output: results/gsea_all.csv, results/gsea_significant.csv, figures/05_*.png
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import gseapy as gp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy import stats  # noqa: E402

from emt import banner, config as C, savefig  # noqa: E402

banner("Step 5 · GSEA of EMT states (one state vs the rest)")
states = pd.read_csv(C.WORK_DIR / "tcga_emt_states.csv", index_col=0)["EMT_state"]
expr = pd.read_pickle(C.WORK_DIR / "tcga_expr_symbols.pkl").loc[states.index]
expr = expr.loc[:, expr.mean() >= 0.5]           # drop essentially unexpressed genes
expr.columns = expr.columns.str.upper()
order = [s for s in C.STATE_NAMES[5] if s in set(states)]
print(f"{expr.shape[0]} tumours × {expr.shape[1]} expressed genes")


def load_library(name):
    C.GENESET_DIR.mkdir(parents=True, exist_ok=True)
    gmt = C.GENESET_DIR / f"{name}.gmt"
    if not gmt.exists():
        lib = gp.get_library(name, organism="Human")
        gmt.write_text("\n".join("\t".join([k, ""] + v) for k, v in lib.items()) + "\n")
    lib = {}
    for line in gmt.read_text().splitlines():
        parts = line.split("\t")
        lib[parts[0]] = [g.upper() for g in parts[2:] if g]
    return lib


libraries = {db: load_library(name) for db, name in C.GSEA_LIBRARIES.items()}
print("Gene set libraries:", {db: len(v) for db, v in libraries.items()})


def ranking(state):
    m = (states == state).values
    t, p = stats.ttest_ind(expr.values[m], expr.values[~m], equal_var=False, axis=0)
    score = np.sign(t) * -np.log10(np.clip(p, 1e-300, 1))
    r = pd.Series(score, index=expr.columns).dropna()
    return r[~r.index.duplicated()].sort_values(ascending=False)


rows = []
for state in order:
    rnk = ranking(state)
    for db, lib in libraries.items():
        res = gp.prerank(rnk=rnk, gene_sets=lib, min_size=C.GSEA_MIN_SIZE, max_size=C.GSEA_MAX_SIZE,
                         permutation_num=C.GSEA_PERMUTATIONS, seed=C.SEED, threads=4, outdir=None,
                         verbose=False).res2d
        res = res.assign(state=state, library=db)
        rows.append(res[["state", "library", "Term", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]])
    print(f"  {state}: done")

gsea = pd.concat(rows, ignore_index=True)
for col in ["NES", "NOM p-val", "FDR q-val"]:
    gsea[col] = pd.to_numeric(gsea[col], errors="coerce")
gsea = gsea.rename(columns={"NOM p-val": "p", "FDR q-val": "FDR"})
gsea.to_csv(C.RESULTS_DIR / "gsea_all.csv", index=False)
sig = gsea[gsea.FDR < C.GSEA_FDR].sort_values(["state", "NES"], ascending=[True, False])
sig.drop(columns="Lead_genes").to_csv(C.RESULTS_DIR / "gsea_significant.csv", index=False)
n_tests = len(gsea)
print(f"\n{len(sig)} of {n_tests} state × gene-set tests at FDR < {C.GSEA_FDR} "
      f"({(gsea.FDR < 0.05).sum()} at FDR < 0.05)")

# ── Hallmark NES heatmap ──────────────────────────────────────────────────────
hall = gsea[gsea.library == "Hallmark"]
keep = hall[hall.FDR < 0.05].Term.value_counts()
keep = keep[keep >= 1].index
nes = hall.pivot(index="Term", columns="state", values="NES").loc[keep, order]
fdr = hall.pivot(index="Term", columns="state", values="FDR").loc[keep, order]
nes = nes.loc[nes.abs().max(axis=1).sort_values(ascending=False).index[:30]]
nes = nes.loc[nes.sort_values("Hybrid" if "Hybrid" in order else order[0]).index]
fig, ax = plt.subplots(figsize=(7.5, max(5, 0.28 * len(nes))))
sns.heatmap(nes, cmap="RdBu_r", center=0, vmin=-3, vmax=3, ax=ax, linewidths=0.3,
            annot=np.where(fdr.loc[nes.index] < 0.05, "•", ""), fmt="", cbar_kws={"label": "NES"})
ax.set(xlabel="", ylabel="", title="Hallmark pathways by EMT state (• FDR < 0.05)")
savefig(fig, "05_hallmark_nes_heatmap.png")

if "Hybrid" in order:
    h = gsea[(gsea.state == "Hybrid") & (gsea.FDR < 0.05)]
    top = pd.concat([h.nlargest(10, "NES"), h.nsmallest(10, "NES")]).drop_duplicates("Term").sort_values("NES")
    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.3 * len(top))))
    ax.barh(top.Term.str.slice(0, 60), top.NES, color=np.where(top.NES > 0, "#B2182B", "#2166AC"))
    ax.axvline(0, color="black", lw=0.8)
    ax.set(xlabel="NES", title="Hybrid state vs rest: top enriched / depleted gene sets (FDR < 0.05)")
    savefig(fig, "05_hybrid_top_pathways.png")
print("Saved: results/gsea_all.csv, results/gsea_significant.csv")
