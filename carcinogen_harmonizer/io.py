from __future__ import annotations

from pathlib import Path
from typing import Iterable
import json

import pandas as pd


REQUIRED_SEED_COLUMNS = [
    "node_id",
    "node_label",
    "node_type",
    "carcinogen_class",
    "iarc",
    "Status",
]


def read_seed_csv(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_SEED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Seed CSV missing required columns: {missing}")
    for col in REQUIRED_SEED_COLUMNS:
        df[col] = df[col].astype(str).str.strip()
    return df


def write_jsonl(records: Iterable[dict], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False, sort_keys=True) + "\n")


def write_json(obj: dict, path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def write_csv(df: pd.DataFrame, path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)


def write_dataframe_json(df: pd.DataFrame, path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(
        json.dumps(df.to_dict(orient="records"), ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
