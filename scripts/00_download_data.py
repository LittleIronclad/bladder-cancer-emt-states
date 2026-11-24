"""Step 0 — Fetch the public input data into data/ (skips files that are already present).

Automatic: TCGA-BLCA expression + GENCODE probemap (UCSC Xena GDC hub), four GEO series (GEOparse).
Manual (license / portal export): TCGA clinical table, GDSC2 drug response, Cell Model Passports.
See data/README.md for details.
"""
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from emt import banner, config as C  # noqa: E402

XENA = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
FILES = {C.TCGA_DIR / "TCGA-BLCA.star_tpm.tsv.gz": XENA + "TCGA-BLCA.star_tpm.tsv.gz",
         C.TCGA_DIR / "gencode.v36.annotation.gtf.gene.probemap": XENA + "gencode.v36.annotation.gtf.gene.probemap"}
GEO = ["GSE32894", "GSE48075", "GSE31684", "GSE13507"]

banner("Step 0 · Download public data")
for dest, url in FILES.items():
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or (dest.name.startswith("TCGA-BLCA.star_tpm") and (C.TCGA_DIR / "expression_matrix.tsv").exists()):
        print(f"  ✓ {dest.relative_to(C.ROOT)} (present)")
        continue
    print(f"  ↓ {url}")
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:  # noqa: BLE001
        print(f"    failed ({e}); download it manually from https://xenabrowser.net "
              "(GDC TCGA Bladder Cancer → 'gene expression RNAseq - STAR - TPM').")

try:
    import GEOparse
    for gse in GEO:
        d = C.GEO_DIR / gse
        d.mkdir(parents=True, exist_ok=True)
        if (d / f"{gse}_family.soft.gz").exists():
            print(f"  ✓ data/geo/{gse} (present)")
        else:
            print(f"  ↓ {gse} from GEO")
            GEOparse.get_GEO(gse, destdir=str(d), silent=True)
except ImportError:
    print("  GEOparse not installed — pip install -r requirements.txt")

missing = []
if not (C.TCGA_DIR / "clinical.tsv").exists():
    missing.append("data/tcga/clinical.tsv — GDC portal → project TCGA-BLCA → Clinical → TSV (unzip clinical.tsv)")
if not list(C.GDSC_DIR.glob("GDSC2_fitted_dose_response_*.xlsx")):
    missing.append("data/gdsc/GDSC2_fitted_dose_response_*.xlsx — https://www.cancerrxgene.org/downloads/bulk_download")
if not list(C.GDSC_DIR.glob("model_list_*.csv")) or not list(C.GDSC_DIR.glob("rnaseq_merged_*.zip")):
    missing.append("data/gdsc/model_list_*.csv and rnaseq_merged_*.zip — https://cellmodelpassports.sanger.ac.uk/downloads")
if missing:
    print("\nManual downloads still needed:\n  - " + "\n  - ".join(missing))
else:
    print("\nAll inputs present.")
