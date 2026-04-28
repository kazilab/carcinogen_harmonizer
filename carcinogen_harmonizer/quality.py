from __future__ import annotations

from typing import Any

import pandas as pd

try:
    from .config import PipelineConfig
except ImportError:
    from config import PipelineConfig


PUBCHEM_FAILURE_STATUSES = {"api_error", "no_pubchem_match", "empty_query"}


def add_quality_flags(df: pd.DataFrame, cfg: PipelineConfig) -> pd.DataFrame:
    out = df.copy()
    id_cols = [c for c in cfg.quality.required_identifier_any if c in out.columns]
    if id_cols:
        out["has_stable_identifier"] = out[id_cols].fillna("").astype(str).apply(lambda r: any(x.strip() for x in r), axis=1)
    else:
        out["has_stable_identifier"] = False
    out["iarc_group_valid"] = out["iarc_group_normalized"].isin(cfg.quality.controlled_iarc_groups)
    out["needs_manual_review"] = False
    if "pubchem_match_status" in out.columns:
        is_chemical_like = out.get("is_chemical_like", pd.Series(True, index=out.index)).astype(str).str.lower().isin(["true", "1", "yes"])
        pubchem_failed = out["pubchem_match_status"].isin(PUBCHEM_FAILURE_STATUSES)
        out["needs_manual_review"] = out["needs_manual_review"] | (is_chemical_like & pubchem_failed)
    out["needs_manual_review"] = out["needs_manual_review"] | (~out["iarc_group_valid"])
    out["needs_manual_review"] = out["needs_manual_review"] | (out["mixture_or_exposure_flag"] == True)
    out["needs_manual_review"] = out["needs_manual_review"] | (out["endogenous_flag"] == True)
    out["confidence_score"] = out.apply(_confidence_score, axis=1)
    return out


def _confidence_score(row: pd.Series) -> str:
    score = 0
    if bool(row.get("has_stable_identifier", False)):
        score += 1
    if bool(row.get("iarc_group_valid", False)):
        score += 1
    pubchem_status = row.get("pubchem_match_status", "matched")
    if pubchem_status == "matched":
        score += 1
    elif pubchem_status == "skipped_non_chemical_agent":
        score += 0
    if bool(row.get("mixture_or_exposure_flag", False)) or bool(row.get("endogenous_flag", False)):
        score -= 1
    if score >= 3:
        return "high"
    if score >= 1:
        return "medium"
    return "low"


def make_qa_report(df: pd.DataFrame) -> dict[str, Any]:
    identifier_cols = [
        "standard_node_id",
        "srandard_node_id",
        "carcinogen_id",
        "pubchem_cid",
        "inchikey",
        "casrn",
        "iarc_casrn",
        "dtxsid",
        "chebi_id",
    ]
    report: dict[str, Any] = {
        "row_count": int(len(df)),
        "columns": list(df.columns),
        "identifier_coverage": _coverage_summary(df, identifier_cols),
        "iarc_group_counts": df.get("iarc_group_normalized", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "agent_entity_type_counts": df.get("agent_entity_type", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "carcinogen_class_counts": df.get("carcinogen_class_normalized", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "carcinogen_class_source_counts": df.get("carcinogen_class_source", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "chemont_class_counts": df.get("chemont_class", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "chemont_superclass_counts": df.get("chemont_superclass", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "classyfire_match_counts": df.get("classyfire_match_status", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "classyfire_validated_count": int(df.get("classyfire_validated", pd.Series(dtype=bool)).fillna(False).astype(bool).sum()),
        "curation_status_counts": df.get("curation_status", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "pubchem_match_counts": df.get("pubchem_match_status", pd.Series(dtype=str)).value_counts(dropna=False).to_dict(),
        "manual_review_count": int(df.get("needs_manual_review", pd.Series(dtype=bool)).sum()),
        "mixture_or_exposure_count": int(df.get("mixture_or_exposure_flag", pd.Series(dtype=bool)).sum()),
        "endogenous_count": int(df.get("endogenous_flag", pd.Series(dtype=bool)).sum()),
    }
    if "standard_node_id" in df.columns:
        report["standard_node_id_unique_count"] = int(df["standard_node_id"].astype(str).nunique(dropna=False))
        report["standard_node_id_duplicate_count"] = int(df["standard_node_id"].astype(str).duplicated(keep=False).sum())
    if "carcinogen_id" in df.columns:
        report["carcinogen_id_unique_count"] = int(df["carcinogen_id"].astype(str).nunique(dropna=False))
        report["carcinogen_id_duplicate_count"] = int(df["carcinogen_id"].astype(str).duplicated(keep=False).sum())
    if "needs_manual_review" in df.columns:
        report["manual_review_ids"] = df.loc[df["needs_manual_review"], "carcinogen_id"].astype(str).tolist()
    return report


def _coverage_summary(df: pd.DataFrame, columns: list[str]) -> dict[str, dict[str, float | int]]:
    total = int(len(df))
    summary: dict[str, dict[str, float | int]] = {}
    for col in columns:
        if col not in df.columns:
            summary[col] = {"present": 0, "missing": total, "coverage": 0.0}
            continue
        series = df[col].fillna("").astype(str).str.strip()
        present = int((series != "").sum())
        summary[col] = {
            "present": present,
            "missing": total - present,
            "coverage": round(present / total, 4) if total else 0.0,
        }
    return summary


def graph_tables(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    node_id = df["carcinogen_id"]
    standard_node_id = df.get("standard_node_id", node_id)
    srandard_node_id = df.get("srandard_node_id", standard_node_id)
    node_label = df.get("display_label", df.get("node_label", ""))
    node_type = df.get("node_type", "Carcinogen")
    carcinogen_class = df.get("carcinogen_class_normalized", "")
    carcinogen_class_source = df.get("carcinogen_class_source", "")
    nodes = pd.DataFrame({
        "standard_node_id": standard_node_id,
        "srandard_node_id": srandard_node_id,
        "node_id": node_id,
        "node_label": node_label,
        "node_type": node_type,
        "id": node_id,
        "label": node_label,
        "type": node_type,
        "agent_entity_type": df.get("agent_entity_type", ""),
        "norm_carcinogen_class": df.get("norm_carcinogen_class", carcinogen_class),
        "carcinogen_class_normalized": carcinogen_class,
        "carcinogen_class_source": carcinogen_class_source,
        "class": carcinogen_class,
        "class_source": carcinogen_class_source,
        "chemont_class": df.get("chemont_class", ""),
        "chemont_superclass": df.get("chemont_superclass", ""),
        "chemont_kingdom": df.get("chemont_kingdom", ""),
        "iarc_group": df.get("iarc_group_normalized", ""),
        "pubchem_cid": df.get("pubchem_cid", ""),
        "inchikey": df.get("inchikey", ""),
        "casrn": df.get("casrn", df.get("iarc_casrn", "")),
    })
    edges = []
    for _, row in df.iterrows():
        cid = row["carcinogen_id"]
        if row.get("iarc_group_normalized"):
            edges.append({
                "source_standard_node_id": cid,
                "source_srandard_node_id": cid,
                "source_node_id": cid,
                "source": cid,
                "predicate": "HAS_IARC_CLASS",
                "target_standard_node_id": row["iarc_group_normalized"],
                "target_srandard_node_id": row["iarc_group_normalized"],
                "target_node_id": row["iarc_group_normalized"],
                "target": row["iarc_group_normalized"],
                "provenance": row.get("iarc_source_url", "seed/IARC") or "seed/IARC",
            })
        if row.get("carcinogen_class_normalized"):
            class_source = str(row.get("carcinogen_class_source", "") or "manual_rules")
            edges.append({
                "source_standard_node_id": cid,
                "source_srandard_node_id": cid,
                "source_node_id": cid,
                "source": cid,
                "predicate": "HAS_AGENT_CLASS",
                "target_standard_node_id": row["carcinogen_class_normalized"],
                "target_srandard_node_id": row["carcinogen_class_normalized"],
                "target_node_id": row["carcinogen_class_normalized"],
                "target": row["carcinogen_class_normalized"],
                "provenance": f"carcinogen_class_source:{class_source}",
            })
        if row.get("agent_entity_type"):
            edges.append({
                "source_standard_node_id": cid,
                "source_srandard_node_id": cid,
                "source_node_id": cid,
                "source": cid,
                "predicate": "HAS_ENTITY_TYPE",
                "target_standard_node_id": row["agent_entity_type"],
                "target_srandard_node_id": row["agent_entity_type"],
                "target_node_id": row["agent_entity_type"],
                "target": row["agent_entity_type"],
                "provenance": "IARC_entity_type_rules",
            })
        if row.get("pubchem_cid"):
            edges.append({
                "source_standard_node_id": cid,
                "source_srandard_node_id": cid,
                "source_node_id": cid,
                "source": cid,
                "predicate": "SAME_AS",
                "target_standard_node_id": f"PubChem:{row['pubchem_cid']}",
                "target_srandard_node_id": f"PubChem:{row['pubchem_cid']}",
                "target_node_id": f"PubChem:{row['pubchem_cid']}",
                "target": f"PubChem:{row['pubchem_cid']}",
                "provenance": "PubChem PUG-REST",
            })
        if row.get("chemont_class"):
            edges.append({
                "source_standard_node_id": cid,
                "source_srandard_node_id": cid,
                "source_node_id": cid,
                "source": cid,
                "predicate": "HAS_CHEMONT_CLASS",
                "target_standard_node_id": row["chemont_class"],
                "target_srandard_node_id": row["chemont_class"],
                "target_node_id": row["chemont_class"],
                "target": row["chemont_class"],
                "provenance": row.get("classyfire_reference_url") or row.get("classyfire_source_url") or "ClassyFire",
            })
        if row.get("chemont_superclass"):
            edges.append({
                "source_standard_node_id": cid,
                "source_srandard_node_id": cid,
                "source_node_id": cid,
                "source": cid,
                "predicate": "HAS_CHEMONT_SUPERCLASS",
                "target_standard_node_id": row["chemont_superclass"],
                "target_srandard_node_id": row["chemont_superclass"],
                "target_node_id": row["chemont_superclass"],
                "target": row["chemont_superclass"],
                "provenance": row.get("classyfire_source_url") or "ClassyFire",
            })
        if row.get("chemont_direct_parent"):
            edges.append({
                "source_standard_node_id": cid,
                "source_srandard_node_id": cid,
                "source_node_id": cid,
                "source": cid,
                "predicate": "HAS_CHEMONT_DIRECT_PARENT",
                "target_standard_node_id": row["chemont_direct_parent"],
                "target_srandard_node_id": row["chemont_direct_parent"],
                "target_node_id": row["chemont_direct_parent"],
                "target": row["chemont_direct_parent"],
                "provenance": row.get("classyfire_source_url") or "ClassyFire",
            })
    return {"nodes_carcinogens": nodes, "edges_carcinogen_facts": pd.DataFrame(edges)}
