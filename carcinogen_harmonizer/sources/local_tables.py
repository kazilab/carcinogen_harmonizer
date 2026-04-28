from __future__ import annotations

from pathlib import Path
from typing import Iterable
import re

import pandas as pd

try:
    from ..normalize import clean_text
except ImportError:
    from normalize import clean_text


def normalize_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", clean_text(value).lower()).strip()


def load_optional_table(path: str | Path | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Reference file not found: {p}")
    if p.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    return pd.read_csv(p, dtype=str, keep_default_na=False)


def build_name_lookup(df: pd.DataFrame, name_columns: Iterable[str]) -> dict[str, pd.Series]:
    lookup: dict[str, pd.Series] = {}
    if df.empty:
        return lookup
    for _, row in df.iterrows():
        for col in name_columns:
            if col in df.columns and str(row.get(col, "")).strip():
                lookup.setdefault(normalize_key(row[col]), row)
    return lookup


def left_enrich_by_names(
    seed: pd.DataFrame,
    ref: pd.DataFrame,
    seed_name_columns: list[str],
    ref_name_columns: list[str],
    output_prefix: str,
    keep_columns: list[str] | None = None,
) -> pd.DataFrame:
    if ref.empty:
        return seed
    lookup = build_name_lookup(ref, ref_name_columns)
    keep_columns = keep_columns or list(ref.columns)
    rows = []
    for _, row in seed.iterrows():
        match = None
        match_key = ""
        for col in seed_name_columns:
            key = normalize_key(row.get(col, ""))
            if key in lookup:
                match = lookup[key]
                match_key = key
                break
        extras = {}
        if match is not None:
            for col in keep_columns:
                if col in ref.columns:
                    extras[f"{output_prefix}_{col}"] = match.get(col, "")
            extras[f"{output_prefix}_match_key"] = match_key
            extras[f"{output_prefix}_match_status"] = "matched"
        else:
            for col in keep_columns:
                extras[f"{output_prefix}_{col}"] = ""
            extras[f"{output_prefix}_match_key"] = ""
            extras[f"{output_prefix}_match_status"] = "no_match"
        rows.append(extras)
    return pd.concat([seed.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
