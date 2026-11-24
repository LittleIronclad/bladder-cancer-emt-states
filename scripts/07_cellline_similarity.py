"""Step 7 — Map chemoresistant cell sublines onto the TCGA EMT states (wet-lab ↔ dry-lab bridge).

Input: data/cellline_qpcr_log2fc.csv — RT-qPCR log2 fold changes (2^-ΔΔCt, GAPDH/ACTB geometric
mean) of each resistant subline relative to its parental line, columns:
    subline, CDH1, CDH2, ZEB1, ZEB2, SNAI1, SNAI2, VIM

Each subline's log2FC vector is compared with directional reference profiles of the five states
(+1 up, −1 down, 0 variable) by cosine similarity, and summarised by
    Δ_EMT = mean(log2FC of mesenchymal genes) − log2FC(CDH1).
This is a coarse, 7-gene, direction-only comparison — a hypothesis generator, not a classifier.

Output: results/cellline_cosine_similarity.csv, figures/07_cellline_similarity.png
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from emt import banner, config as C, savefig  # noqa: E402

banner("Step 7 · Resistant sublines vs TCGA EMT states (cosine similarity)")
path = C.DATA_DIR / "cellline_qpcr_log2fc.csv"
if not path.exists():
    print(f"{path} not found — skipping (the qPCR measurements are not distributed with this repository).")
    sys.exit(0)

GENES = ["CDH1", "CDH2", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "VIM"]
REFERENCE = pd.DataFrame(
    [[+1, -1, -1, -1, -1, -1, -1],
     [+1, -1, -1, -1, 0, 0, -1],
     [0, +1, 0, 0, +1, +1, +1],
     [-1, +1, +1, +1, +1, +1, +1],
     [-1, +1, +1, +1, +1, +1, +1]],
    index=C.STATE_NAMES[5], columns=GENES)

fc = pd.read_csv(path, index_col="subline")[GENES]
A, B = fc.values, REFERENCE.values
cos = pd.DataFrame(A @ B.T / (np.linalg.norm(A, axis=1, keepdims=True) * np.linalg.norm(B, axis=1)),
                   index=fc.index, columns=REFERENCE.index)
delta = fc[GENES[1:]].mean(axis=1) - fc["CDH1"]
out = cos.round(3).assign(best_match=cos.idxmax(axis=1), cosine=cos.max(axis=1).round(3), delta_EMT=delta.round(2))
out.to_csv(C.RESULTS_DIR / "cellline_cosine_similarity.csv")
print(out.to_string())

fig, axes = plt.subplots(1, 2, figsize=(12, 0.6 * len(fc) + 1.8), gridspec_kw={"width_ratios": [1.2, 1]})
sns.heatmap(fc, cmap="RdBu_r", center=0, annot=True, fmt=".1f", ax=axes[0], cbar_kws={"label": "log2 FC vs parental"})
axes[0].set(title="RT-qPCR: EMT markers in resistant sublines", ylabel="")
sns.heatmap(cos, cmap="PuOr_r", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", ax=axes[1],
            cbar_kws={"label": "cosine similarity"})
axes[1].set(title="Similarity to TCGA EMT-state profiles", ylabel="")
savefig(fig, "07_cellline_similarity.png")
