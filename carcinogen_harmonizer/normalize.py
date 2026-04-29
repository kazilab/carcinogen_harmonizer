from __future__ import annotations

from html import unescape
from importlib import resources
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd


_ALIAS_MAP_RESOURCE = "alias_map.json"


def _read_alias_map_json() -> str:
    package = __package__ or "carcinogen_harmonizer"
    try:
        return (
            resources.files(package)
            .joinpath("data", _ALIAS_MAP_RESOURCE)
            .read_text(encoding="utf-8")
        )
    except (AttributeError, FileNotFoundError, ModuleNotFoundError, TypeError):
        return (Path(__file__).with_name("data") / _ALIAS_MAP_RESOURCE).read_text(encoding="utf-8")


def load_alias_map() -> dict[str, str]:
    """Load externally sourced alias overrides for PubChem-friendly lookups."""
    data = json.loads(_read_alias_map_json())
    aliases = data.get("aliases", {})
    out: dict[str, str] = {}
    for alias, entry in aliases.items():
        if isinstance(entry, str):
            preferred_query_name = entry
        else:
            preferred_query_name = str(entry.get("preferred_query_name", "")).strip()
        alias = str(alias).strip()
        if alias and preferred_query_name:
            out[alias] = preferred_query_name
    return out


ALIAS_MAP: dict[str, str] = load_alias_map()

CLASS_RULES: list[tuple[str, str]] = [
    (r"benzo\[a\]pyrene|dimethylbenz\[a\]anthracene|polycyclic aromatic|\bPAH\b", "PAH"),
    (r"nitroso|nitrosamine|NNK|NDMA|NDEA", "Nitrosamine"),
    (r"aldehyde|formaldehyde|acetaldehyde|acrolein|crotonaldehyde|furfural|malondialdehyde|hydroxynonenal", "Aldehyde"),
    (r"PCB|chlorobiphenyl", "PCB"),
    (r"dioxin|dibenzofuran|TCDD|PeCDF", "Dioxin"),
    (r"arsenic|cadmium|beryllium|chromium|nickel|uranium|cobalt|lead", "Heavy_Metal"),
    (r"aflatoxin", "Mycotoxin"),
    (r"aminobiphenyl|benzidine|aromatic amine", "Aromatic_Amine"),
    (r"estradiol|estrogen", "Estrogen"),
    (r"testosterone|dihydrotestosterone|androgen", "Androgen"),
    (r"chlordane|DDT|DDE|heptachlor|hexachlorobenzene|lindane|pentachlorophenol|toxaphene", "Organochlorine"),
    (r"ethylene oxide|sulfur mustard|busulfan|chlorambucil|cyclophosphamide|temozolomide|acrylamide|glycidamide|nitrosourea", "Alkylating"),
    (r"trichloroethylene|tetrachloroethylene", "Chlorinated_Solvent"),
    (r"benzene|vinyl chloride", "Solvent"),
    (r"ethanol|ethyl carbamate|urethane", "Alcohol"),
    (r"MeIQx|PhIP|imidazo", "HCA"),
]

_HTML_TAG_RE = re.compile(r"<[^>]+>")
_TRAILING_PAREN_GROUPS_RE = re.compile(r"(?:\s*\([^()]*\)\s*)+$")


def normalize_iarc_group(value: str) -> str:
    v = str(value or "").strip()
    if not v:
        return ""
    if re.search(r"not classified", v, flags=re.I):
        return "Not classified"
    match = re.search(r"Group\s*(1|2A|2B|3)", v, flags=re.I)
    return f"Group {match.group(1).upper()}" if match else v


def clean_text(value: Any) -> str:
    text = unescape(str(value or "")).replace("\xa0", " ")
    text = _HTML_TAG_RE.sub("", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def strip_trailing_parenthetical_groups(text: str) -> str:
    cleaned = clean_text(text)
    stripped = _TRAILING_PAREN_GROUPS_RE.sub("", cleaned).strip()
    return stripped or cleaned


def is_mixture_or_exposure(label: str, iarc: str, agent_entity_type: str = "") -> bool:
    if str(agent_entity_type or "") in {"Mixture", "Exposure", "Occupation", "Process"}:
        return True
    text = f"{label} {iarc}".lower()
    keywords = [
        "alcoholic beverages",
        "mixture",
        "technical",
        "compounds",
        "inorganic",
        "smoke",
        "exhaust",
        "pollution",
        "occupational",
        "emissions",
        "dust",
        "fumes",
    ]
    return any(k in text for k in keywords)


def is_endogenous(label: str, iarc: str) -> bool:
    text = f"{label} {iarc}".lower()
    return any(k in text for k in ["endogenous", "malondialdehyde", "hydroxynonenal", "4-hne"])


def _coerce_chemical_like(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    return text not in {"false", "0", "no"}


def search_name(row: pd.Series) -> str:
    node_id = clean_text(row.get("node_id", "")).strip()
    label = clean_text(row.get("display_label", row.get("node_label", ""))).strip()
    if node_id in ALIAS_MAP:
        return ALIAS_MAP[node_id]
    if label in ALIAS_MAP:
        return ALIAS_MAP[label]
    cleaned = strip_trailing_parenthetical_groups(label)
    return cleaned or label or node_id


def infer_carcinogen_class(label: str, seed_class: str = "") -> str:
    klass, _source = infer_carcinogen_class_with_source(label, seed_class)
    return klass


def infer_carcinogen_class_with_source(label: str, seed_class: str = "") -> tuple[str, str]:
    """Return ``(class, source)`` where source is one of:

    - ``seed_csv``: the seed/IARC row already had a non-empty class string
    - ``rules``: matched one of the regex patterns in ``CLASS_RULES``
    - ``unknown``: nothing matched
    """
    if str(seed_class or "").strip():
        return str(seed_class).strip(), "seed_csv"
    haystack = label or ""
    for pattern, klass in CLASS_RULES:
        if re.search(pattern, haystack, flags=re.I):
            return klass, "rules"
    return "Unknown", "unknown"


def normalize_seed(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy().fillna("")
    out["carcinogen_id"] = out["node_id"].astype(str).str.strip()
    out["standard_node_id"] = out["carcinogen_id"]
    out["srandard_node_id"] = out["carcinogen_id"]
    out["display_label"] = out["node_label"].map(clean_text)
    out["preferred_query_name"] = out.apply(search_name, axis=1)
    out["iarc_group_normalized"] = out["iarc"].map(normalize_iarc_group)
    out["curation_status"] = out["Status"].astype(str).map({"1": "included", "0": "seed_only"}).fillna("unknown")
    if "expansion_source" in out.columns:
        out.loc[out["expansion_source"].astype(str).eq("iarc"), "curation_status"] = "iarc_expanded"
    if "agent_entity_type" not in out.columns:
        out["agent_entity_type"] = "Chemical"
    out["agent_entity_type"] = out["agent_entity_type"].replace("", "Chemical")
    if "is_chemical_like" not in out.columns:
        out["is_chemical_like"] = True
    out["is_chemical_like"] = out["is_chemical_like"].map(_coerce_chemical_like)
    out["mixture_or_exposure_flag"] = out.apply(
        lambda r: is_mixture_or_exposure(r["display_label"], r["iarc"], r.get("agent_entity_type", "")), axis=1
    )
    out["endogenous_flag"] = out.apply(lambda r: is_endogenous(r["display_label"], r["iarc"]), axis=1)
    def _classify(row: pd.Series) -> tuple[str, str]:
        klass, source = infer_carcinogen_class_with_source(
            row["preferred_query_name"], row.get("carcinogen_class", "")
        )
        # Distinguish curated-seed values from IARC-expanded entity-type
        # values that flowed in via the IARC->seed conversion. Both end
        # up in the ``carcinogen_class`` column, but their provenance is
        # different.
        if source == "seed_csv" and str(row.get("expansion_source", "")).strip() == "iarc":
            return klass, "iarc_entity_type"
        return klass, source

    class_results = out.apply(_classify, axis=1)
    out["carcinogen_class_normalized"] = [c for c, _ in class_results]
    out["carcinogen_class_source"] = [s for _, s in class_results]
    if "source_seed" not in out.columns:
        out["source_seed"] = "node_type_carcinogen.csv"
    out["source_seed"] = out["source_seed"].replace("", "node_type_carcinogen.csv")
    return out


def provenance_for_seed(row: pd.Series) -> dict[str, Any]:
    return {
        "carcinogen_id": row.get("carcinogen_id"),
        "source": row.get("source_seed", "seed_csv"),
        "fields": {
            "node_id": row.get("node_id"),
            "node_label": row.get("node_label"),
            "node_type": row.get("node_type"),
            "carcinogen_class": row.get("carcinogen_class"),
            "iarc": row.get("iarc"),
            "Status": row.get("Status"),
            "agent_entity_type": row.get("agent_entity_type"),
            "iarc_volume": row.get("iarc_volume"),
            "iarc_year": row.get("iarc_year"),
            "iarc_source_url": row.get("iarc_source_url"),
        },
    }
