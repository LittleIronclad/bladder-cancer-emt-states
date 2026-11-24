"""Step 1 — TCGA-BLCA preprocessing.

Input  (data/tcga/): UCSC Xena STAR TPM matrix, GDC clinical.tsv, GENCODE v36 probemap
Output (work/):     tcga_expr_symbols.pkl  primary tumours × gene symbols, log2(TPM+1)
                    tcga_panel.csv         primary tumours × 42 EMT panel genes
                    tcga_clinical.csv      OS (years, censored at 10 y), event, age, stage
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402

from emt import banner, config as C, data  # noqa: E402

banner("Step 1 · TCGA-BLCA preprocessing")
paths = data.tcga_paths()

print("Loading expression:", paths["expr"].name)
expr = data.load_tcga_expression(paths["expr"])
print(f"  {expr.shape[0]} aliquots × {expr.shape[1]} genes")
expr = data.primary_tumour_one_per_patient(expr)
print(f"  primary tumours, one per patient: {len(expr)}")

clin = data.load_gdc_clinical(paths["clinical"])
common = expr.index.intersection(clin.index)
expr, clin = expr.loc[common], clin.loc[common]
print(f"  with overall-survival data: {len(common)} patients, {int(clin.event.sum())} deaths")

surv = data.survival_frame(clin["os_days"] / C.DAYS_PER_YEAR, clin["event"])
clin_out = pd.concat([surv, clin[["age", "stage", "gender"]]], axis=1)
clin_out.index.name = "patient"

missing = [g for g, e in C.ENSEMBL_MAP.items() if e not in expr.columns]
assert not missing, f"panel genes missing from matrix: {missing}"
panel = expr[[C.ENSEMBL_MAP[g] for g in C.PANEL_GENES]]
panel.columns = C.PANEL_GENES
panel.index.name = "patient"

sym = data.to_symbols(expr, data.ensembl_to_symbol(paths["probemap"]))
print(f"  full matrix with gene symbols: {sym.shape[1]} unique genes")

sym.to_pickle(C.WORK_DIR / "tcga_expr_symbols.pkl")
panel.to_csv(C.WORK_DIR / "tcga_panel.csv")
clin_out.to_csv(C.WORK_DIR / "tcga_clinical.csv")
print(f"\nSaved to {C.WORK_DIR}/: tcga_expr_symbols.pkl, tcga_panel.csv, tcga_clinical.csv")
print(f"Median observed time: {surv.time.median():.2f} years")
