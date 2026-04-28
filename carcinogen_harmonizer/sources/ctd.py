from __future__ import annotations

from pathlib import Path
from typing import Iterable
import re

import pandas as pd

try:
    from ..normalize import clean_text
except ImportError:
    from normalize import clean_text


CANCER_TERMS = re.compile(r"cancer|carcinoma|neoplasm|tumou?r|leukemia|lymphoma|melanoma|sarcoma", re.I)


def normalize_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", clean_text(value).lower()).strip()


def read_ctd_tsv(path: str | Path) -> pd.DataFrame:
    """Read CTD TSV/CSV exports, including gzip files and comment headers."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"CTD file not found: {p}")
    sep = "\t" if p.suffix.lower() in {".tsv", ".gz"} else ","
    return pd.read_csv(p, sep=sep, comment="#", dtype=str, keep_default_na=False, compression="infer")


def _first_present(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def summarize_ctd_diseases(seed: pd.DataFrame, diseases_path: str | Path, max_items: int = 20) -> pd.DataFrame:
    ctd = read_ctd_tsv(diseases_path)
    chem_col = _first_present(ctd, ["ChemicalName", "Chemical", "chemical_name"])
    cas_col = _first_present(ctd, ["CasRN", "CASRN", "CAS", "casrn"])
    disease_col = _first_present(ctd, ["DiseaseName", "Disease", "disease_name"])
    pubmed_col = _first_present(ctd, ["PubMedIDs", "PubMedID", "pubmed_ids"])
    if not chem_col or not disease_col:
        raise ValueError("CTD disease file must contain ChemicalName and DiseaseName-like columns")

    by_name: dict[str, pd.DataFrame] = {k: g for k, g in ctd.groupby(ctd[chem_col].map(normalize_key))}
    by_cas: dict[str, pd.DataFrame] = {}
    if cas_col:
        by_cas = {k: g for k, g in ctd.groupby(ctd[cas_col].astype(str).str.strip()) if k}

    additions = []
    for _, row in seed.iterrows():
        frames = []
        cas = str(row.get("casrn", "")).strip()
        if cas and cas in by_cas:
            frames.append(by_cas[cas])
        for key in [normalize_key(row.get("preferred_query_name", "")), normalize_key(row.get("display_label", ""))]:
            if key in by_name:
                frames.append(by_name[key])
        if frames:
            subset = pd.concat(frames, ignore_index=True).drop_duplicates()
            diseases = list(dict.fromkeys(subset[disease_col].astype(str).tolist()))[:max_items]
            cancers = [d for d in diseases if CANCER_TERMS.search(d)][:max_items]
            pubmed_ids = []
            if pubmed_col:
                for val in subset[pubmed_col].astype(str).tolist():
                    for pmid in re.split(r"[|;,]", val):
                        pmid = pmid.strip()
                        if pmid and pmid not in pubmed_ids:
                            pubmed_ids.append(pmid)
                        if len(pubmed_ids) >= max_items:
                            break
                    if len(pubmed_ids) >= max_items:
                        break
            additions.append({
                "ctd_disease_match_status": "matched",
                "associated_diseases": "; ".join(diseases),
                "associated_cancers": "; ".join(cancers),
                "ctd_disease_pubmed_ids": "; ".join(pubmed_ids),
            })
        else:
            additions.append({
                "ctd_disease_match_status": "no_match",
                "associated_diseases": "",
                "associated_cancers": "",
                "ctd_disease_pubmed_ids": "",
            })
    return pd.concat([seed.reset_index(drop=True), pd.DataFrame(additions)], axis=1)


def summarize_ctd_genes(seed: pd.DataFrame, genes_path: str | Path, max_items: int = 30) -> pd.DataFrame:
    ctd = read_ctd_tsv(genes_path)
    chem_col = _first_present(ctd, ["ChemicalName", "Chemical", "chemical_name"])
    cas_col = _first_present(ctd, ["CasRN", "CASRN", "CAS", "casrn"])
    gene_col = _first_present(ctd, ["GeneSymbol", "Gene", "gene_symbol"])
    if not chem_col or not gene_col:
        raise ValueError("CTD gene file must contain ChemicalName and GeneSymbol-like columns")

    by_name: dict[str, pd.DataFrame] = {k: g for k, g in ctd.groupby(ctd[chem_col].map(normalize_key))}
    by_cas: dict[str, pd.DataFrame] = {}
    if cas_col:
        by_cas = {k: g for k, g in ctd.groupby(ctd[cas_col].astype(str).str.strip()) if k}

    additions = []
    for _, row in seed.iterrows():
        frames = []
        cas = str(row.get("casrn", "")).strip()
        if cas and cas in by_cas:
            frames.append(by_cas[cas])
        for key in [normalize_key(row.get("preferred_query_name", "")), normalize_key(row.get("display_label", ""))]:
            if key in by_name:
                frames.append(by_name[key])
        if frames:
            subset = pd.concat(frames, ignore_index=True).drop_duplicates()
            genes = list(dict.fromkeys(subset[gene_col].astype(str).tolist()))[:max_items]
            additions.append({
                "ctd_gene_match_status": "matched",
                "associated_genes": "; ".join(genes),
            })
        else:
            additions.append({
                "ctd_gene_match_status": "no_match",
                "associated_genes": "",
            })
    return pd.concat([seed.reset_index(drop=True), pd.DataFrame(additions)], axis=1)
