"""Step 2 — EMT score and data-driven EMT states in TCGA-BLCA.

1. One-dimensional EMT score from anchor markers (reference axis E → M).
2. NMF on the 3,000 most variable genes (unsupervised transcriptional programmes).
3. Consensus K-means on the NMF sample weights (W) with k = 5 (one cluster per NMF programme);
   silhouette and EMT-score spread are reported for k = 4–8 as a sensitivity check.
4. Clusters are ordered by their anchor-marker profile and named Epithelial → Mesenchymal.

Output: work/tcga_emt_states.csv, results/clustering_summary.json, figures/02_*.png
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
import umap  # noqa: E402
from scipy.cluster.hierarchy import fcluster, linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from scipy.stats import kruskal, spearmanr  # noqa: E402
from sklearn.cluster import KMeans  # noqa: E402
from sklearn.decomposition import NMF  # noqa: E402
from sklearn.metrics import silhouette_score  # noqa: E402

from emt import banner, config as C, savefig  # noqa: E402
from emt.scoring import emt_score, zscore  # noqa: E402

banner("Step 2 · EMT states (NMF + consensus clustering)")
panel = pd.read_csv(C.WORK_DIR / "tcga_panel.csv", index_col=0)
full = pd.read_pickle(C.WORK_DIR / "tcga_expr_symbols.pkl").loc[panel.index]

score = emt_score(panel)
print(f"EMT score: median {score.median():.3f}, IQR {score.quantile(.25):.3f}–{score.quantile(.75):.3f}")

# ── NMF ───────────────────────────────────────────────────────────────────────
hvg = full.var().nlargest(C.N_HVG).index
X = full[hvg].values
print(f"NMF input: {X.shape[0]} tumours × {X.shape[1]} most variable genes")

def nmf(k, seed, max_iter=500):
    return NMF(n_components=k, init="nndsvda", random_state=seed, max_iter=max_iter,
               alpha_W=0.01, alpha_H=0.01, l1_ratio=0.1)

err = {k: [nmf(k, s).fit(X).reconstruction_err_ for s in range(5)] for k in C.NMF_K_RANGE}
fig, ax = plt.subplots(figsize=(6, 3.8))
ks = list(err)
ax.errorbar(ks, [np.mean(err[k]) for k in ks], yerr=[np.std(err[k]) for k in ks], marker="o", capsize=4)
ax.axvline(C.NMF_COMPONENTS, color="crimson", ls="--", lw=1)
ax.set(xlabel="NMF rank", ylabel="Reconstruction error", title="NMF rank selection (elbow)")
savefig(fig, "02_nmf_elbow.png")

best = min((nmf(C.NMF_COMPONENTS, s, 2000).fit(X) for s in range(C.NMF_RESTARTS)),
           key=lambda m: m.reconstruction_err_)
W = best.transform(X)
print(f"NMF rank {C.NMF_COMPONENTS}: best of {C.NMF_RESTARTS} restarts, error {best.reconstruction_err_:.1f}")

# ── Consensus K-means on W ────────────────────────────────────────────────────
def consensus(W, k, runs=C.CONSENSUS_RUNS, frac=C.CONSENSUS_SUBSAMPLE, seed=C.SEED):
    n = len(W)
    rng = np.random.RandomState(seed)
    together, sampled = np.zeros((n, n)), np.zeros((n, n))
    for r in range(runs):
        idx = rng.choice(n, int(n * frac), replace=False)
        lab = KMeans(k, n_init=3, random_state=r, max_iter=200).fit_predict(W[idx])
        onehot = np.eye(k)[lab]
        together[np.ix_(idx, idx)] += onehot @ onehot.T
        sampled[np.ix_(idx, idx)] += 1
    cons = np.divide(together, sampled, out=np.zeros_like(together), where=sampled > 0)
    dist = 1 - cons
    np.fill_diagonal(dist, 0)
    labels = fcluster(linkage(squareform((dist + dist.T) / 2, checks=False), "average"), k, "maxclust") - 1
    return cons, labels

metrics = {}
for k in C.CLUSTER_K_RANGE:
    cons, lab = consensus(W, k)
    means = [score.values[lab == c].mean() for c in range(k)]
    metrics[k] = {"silhouette": silhouette_score(W, lab, metric="cosine"),
                  "emt_range": max(means) - min(means), "labels": lab, "consensus": cons}
    print(f"  k={k}: silhouette={metrics[k]['silhouette']:.3f}, EMT-score range={metrics[k]['emt_range']:.3f}")

k_best = C.CLUSTER_K
labels = metrics[k_best]["labels"]
k_sil = max(metrics, key=lambda k: metrics[k]["silhouette"])
print(f"Using pre-specified k = {k_best} (highest silhouette: k = {k_sil})")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
ks = list(metrics)
for ax, key, title in [(axes[0], "silhouette", "Silhouette (cosine)"), (axes[1], "emt_range", "EMT-score range")]:
    ax.plot(ks, [metrics[k][key] for k in ks], "o-")
    ax.axvline(k_best, color="crimson", ls="--", lw=1)
    ax.set(xlabel="k", title=title, xticks=ks)
savefig(fig, "02_cluster_selection.png")

# ── Order clusters along E → M by anchor markers ──────────────────────────────
z = zscore(panel)
rank_score = {c: z.loc[labels == c, C.M_ANCHOR].values.mean() - z.loc[labels == c, C.E_ANCHOR].values.mean()
              for c in range(k_best)}
order = sorted(rank_score, key=rank_score.get)
names = C.STATE_NAMES.get(k_best, [f"State_{i + 1}" for i in range(k_best)])
cluster_to_state = {c: names[i] for i, c in enumerate(order)}
states = pd.Series([cluster_to_state[c] for c in labels], index=panel.index, name="EMT_state")

kw_p = kruskal(*[score[states == s] for s in names]).pvalue
rho = spearmanr(states.map({s: i for i, s in enumerate(names)}), score).correlation
print("\n" + states.value_counts().reindex(names).to_string())
print(f"EMT score differs between states: Kruskal–Wallis p = {kw_p:.1e}; Spearman ρ(state rank, score) = {rho:.2f}")

# ── Figures ───────────────────────────────────────────────────────────────────
emb = umap.UMAP(n_neighbors=30, min_dist=0.3, metric="cosine", random_state=C.SEED).fit_transform(W)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
sc = axes[0].scatter(emb[:, 0], emb[:, 1], c=score, cmap="RdBu_r", s=14)
fig.colorbar(sc, ax=axes[0], label="EMT score (0 = E, 1 = M)")
for s in names:
    m = (states == s).values
    axes[1].scatter(emb[m, 0], emb[m, 1], s=14, color=C.STATE_COLORS[s], label=f"{s} (n={m.sum()})")
axes[1].legend(frameon=False, fontsize=8, loc="best")
for ax in axes:
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", xticks=[], yticks=[])
axes[0].set_title("EMT score"); axes[1].set_title("NMF-derived EMT states")
savefig(fig, "02_umap_emt_states.png")

fig, ax = plt.subplots(figsize=(7, 3.6))
sns.violinplot(x=states.values, y=score.values, order=names, palette=C.STATE_COLORS, hue=states.values,
               legend=False, ax=ax, cut=0, inner="quartile")
ax.set(ylabel="EMT score", title=f"EMT score by state (Kruskal–Wallis p = {kw_p:.0e})")
savefig(fig, "02_emt_score_by_state.png")

genes = ([g for g in C.GENE_PANEL["epithelial"]] + [g for g in C.GENE_PANEL["intermediate"]]
         + [g for g in C.GENE_PANEL["mesenchymal"]] + C.GENE_PANEL["stemness"])
heat = z[genes].groupby(states).mean().reindex(names)
fig, ax = plt.subplots(figsize=(14, 3.4))
sns.heatmap(heat, cmap="RdBu_r", center=0, vmin=-1.2, vmax=1.2, ax=ax, linewidths=0.3,
            cbar_kws={"label": "mean z-score"})
for x in np.cumsum([len(C.GENE_PANEL[g]) for g in ["epithelial", "intermediate", "mesenchymal"]]):
    ax.axvline(x, color="black", lw=1.5)
ax.set(xlabel="", ylabel="", title="EMT panel expression by state  (epithelial | intermediate | mesenchymal | stemness)")
savefig(fig, "02_panel_heatmap.png")

cons = metrics[k_best]["consensus"]
idx = np.argsort([order.index(c) for c in labels], kind="stable")
fig, ax = plt.subplots(figsize=(5, 4.4))
im = ax.imshow(cons[np.ix_(idx, idx)], cmap="Blues", vmin=0, vmax=1)
fig.colorbar(im, ax=ax, label="co-clustering frequency")
ax.set(xticks=[], yticks=[], title=f"Consensus matrix (k = {k_best})")
savefig(fig, "02_consensus_matrix.png")

# ── Save ──────────────────────────────────────────────────────────────────────
out = pd.DataFrame({"EMT_state": states, "EMT_score": score, "UMAP_1": emb[:, 0], "UMAP_2": emb[:, 1]})
for i in range(W.shape[1]):
    out[f"NMF_{i + 1}"] = W[:, i]
out.to_csv(C.WORK_DIR / "tcga_emt_states.csv")
summary = {
    "n_tumours": int(len(states)), "n_hvg": C.N_HVG, "nmf_rank": C.NMF_COMPONENTS, "k": int(k_best),
    "silhouette": {int(k): round(float(v["silhouette"]), 3) for k, v in metrics.items()},
    "emt_score_range": {int(k): round(float(v["emt_range"]), 3) for k, v in metrics.items()},
    "mean_within_state_consensus": {cluster_to_state[c]: round(float(
        metrics[k_best]["consensus"][np.ix_(labels == c, labels == c)].mean()), 3) for c in order},
    "state_sizes": states.value_counts().reindex(names).astype(int).to_dict(),
    "state_anchor_profile": {cluster_to_state[c]: round(float(rank_score[c]), 3) for c in order},
    "kruskal_p_emt_score": float(kw_p), "spearman_rho_state_vs_score": round(float(rho), 3),
}
(C.RESULTS_DIR / "clustering_summary.json").write_text(json.dumps(summary, indent=2))
print("\nSaved: work/tcga_emt_states.csv, results/clustering_summary.json")
