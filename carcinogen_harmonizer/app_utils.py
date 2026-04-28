from __future__ import annotations

from pathlib import Path

import pandas as pd


IARC_GROUP_ORDER = ["Group 1", "Group 2A", "Group 2B", "Group 3", "Not classified"]
SEARCH_COLUMNS = [
    "carcinogen_id",
    "standard_node_id",
    "display_label",
    "preferred_query_name",
    "pubchem_lookup_name",
    "pubchem_synonyms",
    "iarc_agent_name",
    "node_label",
    "iupac_name",
    "casrn",
    "iarc_casrn",
    "pubchem_cid",
    "inchikey",
]
DEFAULT_RESULT_COLUMNS = [
    "carcinogen_id",
    "display_label",
    "preferred_query_name",
    "iarc_group_normalized",
    "carcinogen_class_normalized",
    "agent_entity_type",
    "pubchem_cid",
    "casrn",
    "pubchem_match_status",
    "confidence_score",
    "needs_manual_review",
]


def read_static_dataset(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Static dataset not found: {p}")
    return pd.read_csv(p, dtype=str, keep_default_na=False)


def available_columns(df: pd.DataFrame, candidates: list[str]) -> list[str]:
    return [col for col in candidates if col in df.columns]


def iarc_group_summary(df: pd.DataFrame) -> pd.DataFrame:
    if "iarc_group_normalized" not in df.columns:
        return pd.DataFrame(columns=["iarc_group_normalized", "count"])
    counts = df["iarc_group_normalized"].fillna("").astype(str).value_counts(dropna=False)
    rows = [{"iarc_group_normalized": group, "count": int(count)} for group, count in counts.items()]
    order = {group: idx for idx, group in enumerate(IARC_GROUP_ORDER)}
    return (
        pd.DataFrame(rows)
        .assign(_order=lambda d: d["iarc_group_normalized"].map(order).fillna(len(order)))
        .sort_values(["_order", "iarc_group_normalized"])
        .drop(columns=["_order"])
        .reset_index(drop=True)
    )


def sorted_unique_values(df: pd.DataFrame, column: str) -> list[str]:
    if column not in df.columns:
        return []
    values = df[column].fillna("").astype(str).str.strip()
    return sorted(v for v in values.unique() if v)


def filter_dataset(
    df: pd.DataFrame,
    query: str = "",
    *,
    groups: list[str] | None = None,
    classes: list[str] | None = None,
    entity_types: list[str] | None = None,
    manual_review: str = "All",
) -> pd.DataFrame:
    out = df.copy()
    query = str(query or "").strip().lower()
    if query:
        cols = available_columns(out, SEARCH_COLUMNS)
        if cols:
            haystack = out[cols].fillna("").astype(str).agg(" ".join, axis=1).str.lower()
            terms = [term for term in query.split() if term]
            for term in terms:
                out = out[haystack.loc[out.index].str.contains(term, regex=False)]

    if groups and "iarc_group_normalized" in out.columns:
        out = out[out["iarc_group_normalized"].isin(groups)]
    if classes and "carcinogen_class_normalized" in out.columns:
        out = out[out["carcinogen_class_normalized"].isin(classes)]
    if entity_types and "agent_entity_type" in out.columns:
        out = out[out["agent_entity_type"].isin(entity_types)]
    if manual_review != "All" and "needs_manual_review" in out.columns:
        expected = manual_review == "Needs review"
        review = out["needs_manual_review"].fillna("").astype(str).str.lower().isin(["true", "1", "yes"])
        out = out[review == expected]
    return out.reset_index(drop=True)
