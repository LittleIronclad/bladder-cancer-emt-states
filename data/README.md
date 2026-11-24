# Input data

All inputs are public. They are not tracked in git. `python scripts/00_download_data.py` fetches what can be downloaded automatically and lists the rest.

```
data/
├── tcga/
│   ├── TCGA-BLCA.star_tpm.tsv.gz                 UCSC Xena, GDC hub — log2(TPM+1), STAR, GENCODE v36
│   ├── gencode.v36.annotation.gtf.gene.probemap  UCSC Xena, GDC hub — Ensembl → gene symbol
│   └── clinical.tsv                              GDC portal → TCGA-BLCA → Clinical → TSV
├── geo/
│   └── GSE32894/, GSE48075/, GSE31684/, GSE13507/   *_family.soft.gz (downloaded by GEOparse)
├── gdsc/
│   ├── GDSC2_fitted_dose_response_27Oct23.xlsx   cancerrxgene.org → Downloads → Bulk data
│   ├── model_list_<date>.csv                     cellmodelpassports.sanger.ac.uk → Downloads
│   └── rnaseq_merged_<date>.zip                  cellmodelpassports.sanger.ac.uk → Downloads (RNA-seq, all models)
└── genesets/                                     Enrichr libraries cached by step 5 (or your own GMT files)
```

| Dataset | Reference | Used for |
|---|---|---|
| TCGA-BLCA | Robertson et al., *Cell* 2017 | discovery cohort (n = 400 primary tumours with OS) |
| GSE32894 | Sjödahl et al., *Clin Cancer Res* 2012 | validation (disease-specific death) |
| GSE48075 | Choi et al., *Cancer Cell* 2014 | validation (OS; cT2+ only) |
| GSE31684 | Riester et al., *Clin Cancer Res* 2012 | validation (disease-specific death) |
| GSE13507 | Kim et al., *Mol Cancer* 2010 | validation (OS; primary tumours) |
| GDSC2 / Cell Model Passports | Iorio et al., *Cell* 2016; van der Meer et al., *NAR* 2019 | cell-line drug response |

`expression_matrix.tsv` (the file name used by the 2026 version) is accepted in place of `TCGA-BLCA.star_tpm.tsv.gz`.
