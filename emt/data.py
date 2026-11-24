"""Loaders for TCGA-BLCA (UCSC Xena / GDC) and GEO bladder cancer cohorts."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C


def _first_existing(*paths: Path) -> Path:
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("None of these files exist: " + ", ".join(map(str, paths))
                            + "\nSee data/README.md for download instructions.")


# ── TCGA-BLCA ─────────────────────────────────────────────────────────────────
def tcga_paths():
    d = C.TCGA_DIR
    return {
        "expr": _first_existing(d / "TCGA-BLCA.star_tpm.tsv.gz", d / "TCGA-BLCA.star_tpm.tsv",
                                d / "expression_matrix.tsv"),
        "clinical": _first_existing(d / "clinical.tsv"),
        "probemap": _first_existing(d / "gencode.v36.annotation.gtf.gene.probemap"),
    }


def load_tcga_expression(path: Path) -> pd.DataFrame:
    """UCSC Xena GDC-hub STAR TPM matrix → samples × Ensembl IDs (version suffix removed).

    The Xena file is *already* log2(TPM+1). Raw TPM is detected by its range and transformed;
    log-scale input is left untouched (the original 2026 code applied log2 a second time here).
    """
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.str.replace(r"\.\d+$", "", regex=True)
    expr = expr.groupby(level=0).mean()
    vmax = np.nanmax(expr.values)
    if vmax > 50:                       # raw TPM (log2(TPM+1) never exceeds ~20)
        print(f"  input looks like raw TPM (max={vmax:.0f}) → log2(TPM+1)")
        expr = np.log2(expr.clip(lower=0) + 1)
    else:
        print(f"  input is already log2(TPM+1) (max={vmax:.1f}) → no transform")
    return expr.T


def primary_tumour_one_per_patient(expr: pd.DataFrame) -> pd.DataFrame:
    """Keep primary tumours (sample type 01), one aliquot per patient (01A before 01B)."""
    types = expr.index.str.split("-").str[3].str[:2]
    prim = expr[types == "01"].sort_index()
    patient = prim.index.str[:12]
    prim = prim[~patient.duplicated(keep="first")]
    prim.index = prim.index.str[:12]
    return prim


def load_gdc_clinical(path: Path) -> pd.DataFrame:
    """GDC portal clinical.tsv export → one row per patient with OS (days), event, age, stage."""
    clin = pd.read_csv(path, sep="\t", low_memory=False, na_values=["'--"])
    n_valid = clin.notna().sum(axis=1)
    clin = clin.loc[n_valid.sort_values(ascending=False, kind="stable").index]
    clin = clin.drop_duplicates("cases.submitter_id")
    out = pd.DataFrame(index=clin["cases.submitter_id"].str[:12].values)
    vital = clin["demographic.vital_status"].str.strip().str.lower().map({"dead": 1, "alive": 0}).values
    d_death = pd.to_numeric(clin["demographic.days_to_death"], errors="coerce").values
    d_fu = pd.to_numeric(clin["diagnoses.days_to_last_follow_up"], errors="coerce").values
    out["event"] = vital
    out["os_days"] = np.where(vital == 1, d_death, d_fu)
    out["age"] = pd.to_numeric(clin["demographic.age_at_index"], errors="coerce").values
    out["stage"] = clin["diagnoses.ajcc_pathologic_stage"].values
    out["gender"] = clin["demographic.gender"].values
    out = out[~out.index.duplicated()]
    return out.dropna(subset=["event", "os_days"]).query("os_days > 0")


def ensembl_to_symbol(probemap: Path) -> dict:
    pm = pd.read_csv(probemap, sep="\t")
    pm["ens"] = pm["id"].str.replace(r"\.\d+$", "", regex=True)
    return pm.set_index("ens")["gene"].to_dict()


def to_symbols(expr: pd.DataFrame, e2g: dict) -> pd.DataFrame:
    """Ensembl columns → gene symbols; duplicated symbols keep the most highly expressed gene."""
    cols = [c for c in expr.columns if c in e2g]
    sym = expr[cols].copy()
    sym.columns = [e2g[c] for c in cols]
    sym = sym.iloc[:, np.argsort(-sym.mean().values, kind="stable")]     # highest mean first
    return sym.loc[:, ~sym.columns.duplicated(keep="first")]


# ── Survival helpers ──────────────────────────────────────────────────────────
def survival_frame(time_years: pd.Series, event: pd.Series, horizon=C.MAX_FOLLOWUP_YEARS) -> pd.DataFrame:
    """Administrative censoring at `horizon` years: events after the horizon become censored."""
    t = pd.to_numeric(time_years, errors="coerce")
    e = pd.to_numeric(event, errors="coerce")
    df = pd.DataFrame({"time": t.clip(upper=horizon), "event": e})
    df.loc[t > horizon, "event"] = 0
    return df


# ── GEO ───────────────────────────────────────────────────────────────────────
SYMBOL_COLUMNS = ["Symbol", "Gene Symbol", "GENE_SYMBOL", "gene_symbol", "ILMN_Gene", "Gene symbol"]


def load_geo(gse_id: str):
    """Return (expression samples × gene symbols, sample metadata) for a GEO series.

    Uses the cached SOFT file in data/geo/<GSE>/ (downloaded on first use by GEOparse).
    Probes are collapsed to genes by mean; data are log2-transformed if not already on log scale.
    """
    import GEOparse

    dest = C.GEO_DIR / gse_id
    dest.mkdir(parents=True, exist_ok=True)
    soft = dest / f"{gse_id}_family.soft.gz"
    gse = (GEOparse.get_GEO(filepath=str(soft), silent=True) if soft.exists()
           else GEOparse.get_GEO(gse_id, destdir=str(dest), silent=True))

    meta, expr = {}, {}
    for name, gsm in gse.gsms.items():
        chars = {}
        for ch in gsm.metadata.get("characteristics_ch1", []):
            if ":" in ch:
                k, v = ch.split(":", 1)
                chars[k.strip()] = v.strip()
        meta[name] = {"title": gsm.metadata.get("title", [""])[0],
                      "source": gsm.metadata.get("source_name_ch1", [""])[0], **chars}
        if gsm.table is not None and not gsm.table.empty:
            expr[name] = gsm.table.set_index("ID_REF")["VALUE"]
    meta = pd.DataFrame(meta).T
    raw = pd.DataFrame(expr).T.apply(pd.to_numeric, errors="coerce")
    raw.columns = raw.columns.astype(str)

    gpl = next(iter(gse.gpls.values())).table
    sym_col = next((c for c in SYMBOL_COLUMNS if c in gpl.columns), None)
    if sym_col is None:
        raise ValueError(f"{gse_id}: no gene symbol column in platform table ({list(gpl.columns)[:10]})")
    sym = gpl.set_index(gpl["ID"].astype(str))[sym_col].astype(str).str.strip()
    sym = sym.str.split(r" /// | // |/|\||;").str[0].str.strip()
    sym = sym[~sym.isin(["", "nan", "---", "NA"])]
    cols = [c for c in raw.columns if c in sym.index]
    genes = raw[cols].T.groupby(sym.loc[cols].values).mean().T

    if np.nanmedian(genes.values) > 20:
        genes = np.log2(genes.clip(lower=1))
    return genes, meta


def _years(series, unit="months"):
    v = pd.to_numeric(series, errors="coerce")
    return v / 12.0 if unit == "months" else v


def geo_cohort(gse_id: str):
    """Expression + survival for the four validation cohorts, with cohort-specific parsing.

    Returns (expr, surv) where surv has columns time (years, capped), event, endpoint.
    """
    expr, meta = load_geo(gse_id)
    m = meta.copy()

    if gse_id == "GSE32894":               # Sjödahl 2012, Illumina HT-12; disease-specific death
        excl = m["dod_excluded_(yes/no)"].str.lower().str.strip().eq("yes")
        m = m[~excl]
        t = _years(m["time_to_dod_(months)"])
        e = m["dod_event_(yes/no)"].str.lower().str.strip().map({"yes": 1, "no": 0})
        endpoint = "DSS"
    elif gse_id == "GSE48075":             # Choi 2014, Illumina HT-12 v3; overall survival; cT2+ only
        mibc = m["cstage"].astype(str).str.contains(r"T[234]", case=False, na=False)
        m = m[mibc]
        t = _years(m["survival (mo)"])
        e = m["os censor"].str.lower().str.strip().map({"uncensored": 1, "censored": 0})
        endpoint = "OS"
    elif gse_id == "GSE31684":             # Riester 2012, Affymetrix U133 Plus 2; cystectomy (DOD = event)
        t = _years(m["survival.months"])
        status = m["last known status"].str.upper().str.strip()
        e = status.map(lambda s: 1 if s == "DOD" else (0 if s in ("NED", "DOC", "AWD") else np.nan))
        endpoint = "DSS"
    elif gse_id == "GSE13507":             # Kim 2010, Illumina Human-6 v2; primary tumours (NMIBC + MIBC)
        prim = m["source"].str.contains("primary", case=False, na=False)
        m = m[prim]
        t = _years(m["survival month"])
        e = m["overall survival"].str.lower().str.strip().map({"death": 1, "survival": 0})
        endpoint = "OS"
    else:
        raise ValueError(gse_id)

    surv = survival_frame(t, e).dropna()
    surv = surv[surv["time"] > 0]
    surv["endpoint"] = endpoint
    common = expr.index.intersection(surv.index)
    return expr.loc[common], surv.loc[common]
