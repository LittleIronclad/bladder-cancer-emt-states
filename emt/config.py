"""Central configuration: paths, gene panel, model parameters and plotting colours."""
import os
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("EMT_DATA_DIR", ROOT / "data"))
WORK_DIR = Path(os.environ.get("EMT_WORK_DIR", ROOT / "work"))        # large intermediates (git-ignored)
RESULTS_DIR = ROOT / "results"                                          # small tables (committed)
FIG_DIR = ROOT / "figures"                                              # figures (committed)

TCGA_DIR = DATA_DIR / "tcga"
GEO_DIR = DATA_DIR / "geo"
GDSC_DIR = DATA_DIR / "gdsc"
GENESET_DIR = DATA_DIR / "genesets"

for _d in (WORK_DIR, RESULTS_DIR, FIG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ── General ───────────────────────────────────────────────────────────────────
SEED = 42
DAYS_PER_YEAR = 365.25
MAX_FOLLOWUP_YEARS = 10.0       # administrative censoring horizon
LANDMARK_YEARS = 2.0

# ── EMT gene panel (42 genes) ─────────────────────────────────────────────────
GENE_PANEL = {
    "epithelial": ["CDH1", "OCLN", "EPCAM", "TJP1", "CLDN4", "KRT8", "KRT18", "KRT19", "GRHL2", "ESRP1"],
    "mesenchymal": ["CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2", "FN1",
                    "COL1A1", "COL3A1", "MMP2", "MMP9", "FOXC2"],
    "intermediate": ["SFN", "NRG1", "CD44", "CTNNB1", "PDPN", "ANXA1", "AXL", "MET", "EGFR",
                     "NOTCH1", "TGFB1", "WNT5A", "HMGA2", "FOSL1"],
    "stemness": ["ALDH1A1", "PROM1", "SOX2", "NANOG"],
}
PANEL_GENES = [g for genes in GENE_PANEL.values() for g in genes]

# Anchor markers used for the one-dimensional EMT score and for ordering clusters along E→M.
E_ANCHOR = ["CDH1", "OCLN", "EPCAM", "KRT18", "GRHL2", "ESRP1", "KRT8", "KRT19"]
M_ANCHOR = ["VIM", "CDH2", "ZEB1", "ZEB2", "FN1", "COL1A1", "SNAI2", "TWIST1"]

# Ensembl IDs (GENCODE v36) — used to pull the panel from the TCGA matrix without ambiguity.
ENSEMBL_MAP = {
    "CDH1": "ENSG00000039068", "CDH2": "ENSG00000170558", "OCLN": "ENSG00000197822",
    "EPCAM": "ENSG00000119888", "TJP1": "ENSG00000104067", "CLDN4": "ENSG00000189143",
    "KRT8": "ENSG00000170421", "KRT18": "ENSG00000111057", "KRT19": "ENSG00000171345",
    "GRHL2": "ENSG00000083307", "ESRP1": "ENSG00000104413", "VIM": "ENSG00000026025",
    "SNAI1": "ENSG00000124216", "SNAI2": "ENSG00000019549", "ZEB1": "ENSG00000148516",
    "ZEB2": "ENSG00000169554", "TWIST1": "ENSG00000122691", "TWIST2": "ENSG00000233608",
    "FN1": "ENSG00000115414", "COL1A1": "ENSG00000108821", "COL3A1": "ENSG00000168542",
    "MMP2": "ENSG00000087245", "MMP9": "ENSG00000100985", "FOXC2": "ENSG00000176692",
    "SFN": "ENSG00000175793", "NRG1": "ENSG00000157168", "CD44": "ENSG00000026508",
    "CTNNB1": "ENSG00000168036", "PDPN": "ENSG00000162493", "ANXA1": "ENSG00000135046",
    "AXL": "ENSG00000167601", "MET": "ENSG00000105976", "EGFR": "ENSG00000146648",
    "NOTCH1": "ENSG00000148400", "TGFB1": "ENSG00000105329", "WNT5A": "ENSG00000114251",
    "HMGA2": "ENSG00000149948", "FOSL1": "ENSG00000175592", "ALDH1A1": "ENSG00000165092",
    "PROM1": "ENSG00000007062", "SOX2": "ENSG00000181449", "NANOG": "ENSG00000111704",
}

# ── EMT state discovery (NMF + consensus clustering) ──────────────────────────
N_HVG = 3000                    # most variable genes used for NMF
NMF_K_RANGE = range(2, 11)      # elbow curve
NMF_COMPONENTS = 5              # chosen from the elbow curve
NMF_RESTARTS = 30
CLUSTER_K_RANGE = range(4, 9)  # evaluated and reported
CLUSTER_K = 5                   # pre-specified: one cluster per NMF programme (as in the 2026 version).
                                # Silhouette slightly prefers k = 4, which splits the Hybrid state
                                # between Partial-E and Partial-M — see README, "Sensitivity".
CONSENSUS_RUNS = 100
CONSENSUS_SUBSAMPLE = 0.8

STATE_NAMES = {
    4: ["Epithelial", "Partial-E", "Partial-M", "Mesenchymal"],
    5: ["Epithelial", "Partial-E", "Hybrid", "Partial-M", "Mesenchymal"],
    6: ["Epithelial", "Partial-E", "Hybrid-E", "Hybrid-M", "Partial-M", "Mesenchymal"],
}
STATE_COLORS = {
    "Epithelial": "#2166AC", "Partial-E": "#74ADD1", "Hybrid-E": "#ABD9E9", "Hybrid": "#7F7F7F",
    "Hybrid-M": "#FDAE61", "Partial-M": "#F46D43", "Mesenchymal": "#B2182B",
}
# Three functional groups used for survival analysis
BROAD_ORDER = ["Epithelial + Partial-E", "Hybrid", "Partial-M + Mesenchymal"]
BROAD_COLORS = dict(zip(BROAD_ORDER, ["#2166AC", "#7F7F7F", "#B2182B"]))


def broad_group(state: str) -> str:
    if state in ("Epithelial", "Partial-E"):
        return "Epithelial + Partial-E"
    if state == "Hybrid":
        return "Hybrid"
    return "Partial-M + Mesenchymal"


# ── LASSO Cox ─────────────────────────────────────────────────────────────────
LASSO_ALPHAS = tuple(10 ** x for x in [i * 5 / 59 - 4 for i in range(60)])   # logspace(-4, 1, 60)
LASSO_CV_FOLDS = 5
NESTED_CV_REPEATS = 5
N_BOOTSTRAP = 1000

# ── GSEA ──────────────────────────────────────────────────────────────────────
GSEA_LIBRARIES = {                       # Enrichr library names (downloaded by gseapy)
    "Hallmark": "MSigDB_Hallmark_2020",
    "KEGG": "KEGG_2021_Human",
    "GO_BP": "GO_Biological_Process_2023",
}
GSEA_PERMUTATIONS = 1000
GSEA_MIN_SIZE, GSEA_MAX_SIZE = 15, 500
GSEA_FDR = 0.25

# ── GDSC ──────────────────────────────────────────────────────────────────────
GDSC_DRUGS = ["Cisplatin", "Gemcitabine", "Doxorubicin", "Paclitaxel", "Carboplatin", "Vinblastine",
              "Vorinostat", "Trametinib", "Bortezomib", "Erlotinib", "Gefitinib", "Lapatinib"]
GDSC_MIN_LINES = 10
