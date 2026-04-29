from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

try:
    from .config import PipelineConfig
    from .normalize import clean_text, normalize_seed, provenance_for_seed
    from .sources.pubchem import PubChemClient
    from .sources.classyfire import (
        ClassyFireClient,
        chemont_reference_url,
        chemont_to_class_bucket,
    )
    from .sources.local_tables import load_optional_table, left_enrich_by_names
    from .sources.ctd import summarize_ctd_diseases, summarize_ctd_genes
except ImportError:
    from config import PipelineConfig
    from normalize import clean_text, normalize_seed, provenance_for_seed
    from sources.pubchem import PubChemClient
    from sources.classyfire import (
        ClassyFireClient,
        chemont_reference_url,
        chemont_to_class_bucket,
    )
    from sources.local_tables import load_optional_table, left_enrich_by_names
    from sources.ctd import summarize_ctd_diseases, summarize_ctd_genes


PUBCHEM_COLUMN_MAP = {
    "query": "pubchem_lookup_name",
    "pubchem_cid": "pubchem_cid",
    "molecular_formula": "molecular_formula",
    "molecular_weight": "molecular_weight",
    "canonical_smiles": "canonical_smiles",
    "isomeric_smiles": "isomeric_smiles",
    "inchi": "inchi",
    "inchikey": "inchikey",
    "iupac_name": "iupac_name",
    "casrn": "casrn",
    "synonyms": "pubchem_synonyms",
    "xlogp": "xlogp",
    "tpsa": "tpsa",
    "match_status": "pubchem_match_status",
    "error": "pubchem_error",
    "source_url": "pubchem_source_url",
}


def _format_pubchem_value(value: Any) -> str:
    if isinstance(value, list):
        return "; ".join(str(v) for v in value if str(v).strip())
    return str(value or "")


def _pubchem_query_candidates(row: pd.Series) -> list[str]:
    raw_candidates = [
        row.get("preferred_query_name", ""),
        row.get("display_label", ""),
        row.get("node_label", ""),
        row.get("iarc_agent_name", ""),
        row.get("casrn", ""),
        row.get("iarc_casrn", ""),
    ]
    candidates: list[str] = []
    seen: set[str] = set()
    for raw in raw_candidates:
        candidate = clean_text(raw)
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        candidates.append(candidate)
    return candidates


def _empty_pubchem_record(status: str, error: str = "") -> dict[str, str]:
    rec = {out: "" for out in PUBCHEM_COLUMN_MAP.values()}
    rec["pubchem_match_status"] = status
    rec["pubchem_error"] = error
    return rec


def enrich_pubchem(df: pd.DataFrame, cfg: PipelineConfig, cache_dir: str | Path) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if not cfg.pubchem.enabled:
        return df, []
    client = PubChemClient(
        cache_dir=cache_dir,
        timeout=cfg.pubchem.timeout_seconds,
        rate_limit_seconds=cfg.pubchem.rate_limit_seconds,
    )
    records: list[dict[str, Any]] = []
    provenance: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        candidates = _pubchem_query_candidates(row)
        query = candidates[0] if candidates else ""
        is_chemical_like = str(row.get("is_chemical_like", "True")).lower() not in {"false", "0", "no"}
        if not is_chemical_like:
            records.append(_empty_pubchem_record("skipped_non_chemical_agent"))
            provenance.append({
                "carcinogen_id": row.get("carcinogen_id"),
                "source": "PubChem PUG-REST",
                "query": query,
                "query_candidates": candidates,
                "match_status": "skipped_non_chemical_agent",
                "pubchem_cid": "",
                "source_url": "",
                "error": "",
            })
            continue
        result = None
        used_query = ""
        for candidate in candidates or [""]:
            candidate_result = client.enrich_name(
                candidate,
                properties=cfg.pubchem.properties,
                max_synonyms=cfg.pubchem.max_synonyms,
            )
            result = candidate_result
            used_query = candidate
            if candidate_result.match_status == "matched":
                break
            if candidate_result.match_status in {"no_pubchem_match", "empty_query"}:
                continue
            break
        rec = result.to_record()
        out_rec = {
            out: _format_pubchem_value(rec.get(src, ""))
            for src, out in PUBCHEM_COLUMN_MAP.items()
        }
        out_rec["pubchem_lookup_name"] = used_query
        if not out_rec.get("casrn") and row.get("iarc_casrn"):
            out_rec["casrn"] = row.get("iarc_casrn")
        records.append(out_rec)
        provenance.append({
            "carcinogen_id": row.get("carcinogen_id"),
            "source": "PubChem PUG-REST",
            "query": used_query,
            "query_candidates": candidates,
            "match_status": result.match_status,
            "pubchem_cid": result.pubchem_cid,
            "source_url": result.source_url,
            "error": result.error,
        })
    enriched = pd.concat([df.reset_index(drop=True), pd.DataFrame(records)], axis=1)
    return enriched, provenance


CLASSYFIRE_COLUMNS = [
    "chemont_kingdom",
    "chemont_kingdom_id",
    "chemont_superclass",
    "chemont_superclass_id",
    "chemont_class",
    "chemont_class_id",
    "chemont_subclass",
    "chemont_subclass_id",
    "chemont_direct_parent",
    "chemont_direct_parent_id",
    "chemont_alternative_parents",
    "classyfire_match_status",
    "classyfire_source_url",
    "classyfire_reference_url",
    "classyfire_error",
]


def _empty_classyfire_record(status: str, error: str = "") -> dict[str, str]:
    rec: dict[str, str] = {col: "" for col in CLASSYFIRE_COLUMNS}
    rec["classyfire_match_status"] = status
    rec["classyfire_error"] = error
    return rec


def enrich_classyfire(
    df: pd.DataFrame, cfg: PipelineConfig, cache_dir: str | Path
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Look up ClassyFire (ChemOnt) chemical taxonomy by InChIKey.

    Adds the ``chemont_*`` and ``classyfire_*`` columns and returns
    provenance entries describing the API call status. Rows without an
    InChIKey or that are not chemical-like are skipped with a status
    code so they are still traceable in the QA report.
    """
    if not cfg.classyfire.enabled:
        empty = pd.DataFrame([_empty_classyfire_record("disabled") for _ in range(len(df))])
        merged = pd.concat([df.reset_index(drop=True), empty], axis=1)
        return merged, []

    client = ClassyFireClient(
        cache_dir=cache_dir,
        timeout=cfg.classyfire.timeout_seconds,
        rate_limit_seconds=cfg.classyfire.rate_limit_seconds,
        max_retries=cfg.classyfire.max_retries,
        backoff_base_seconds=cfg.classyfire.backoff_base_seconds,
    )
    records: list[dict[str, str]] = []
    provenance: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        cid = row.get("carcinogen_id", "")
        is_chemical_like = (
            str(row.get("is_chemical_like", "True")).lower()
            not in {"false", "0", "no"}
        )
        inchikey = str(row.get("inchikey", "") or "").strip()
        if not is_chemical_like:
            records.append(_empty_classyfire_record("skipped_non_chemical_agent"))
            provenance.append(
                {
                    "carcinogen_id": cid,
                    "source": "ClassyFire",
                    "query_inchikey": "",
                    "match_status": "skipped_non_chemical_agent",
                    "chemont_class": "",
                    "source_url": "",
                    "error": "",
                }
            )
            continue
        if not inchikey:
            records.append(_empty_classyfire_record("skipped_no_inchikey"))
            provenance.append(
                {
                    "carcinogen_id": cid,
                    "source": "ClassyFire",
                    "query_inchikey": "",
                    "match_status": "skipped_no_inchikey",
                    "chemont_class": "",
                    "source_url": "",
                    "error": "",
                }
            )
            continue
        result = client.lookup_inchikey(inchikey)
        rec = {
            "chemont_kingdom": result.kingdom,
            "chemont_kingdom_id": result.kingdom_id,
            "chemont_superclass": result.superclass,
            "chemont_superclass_id": result.superclass_id,
            "chemont_class": result.chemical_class,
            "chemont_class_id": result.chemical_class_id,
            "chemont_subclass": result.subclass,
            "chemont_subclass_id": result.subclass_id,
            "chemont_direct_parent": result.direct_parent,
            "chemont_direct_parent_id": result.direct_parent_id,
            "chemont_alternative_parents": "; ".join(result.alternative_parents),
            "classyfire_match_status": result.match_status,
            "classyfire_source_url": result.source_url,
            "classyfire_reference_url": chemont_reference_url(
                result.chemical_class_id or result.direct_parent_id or result.superclass_id
            ),
            "classyfire_error": result.error,
        }
        records.append(rec)
        provenance.append(
            {
                "carcinogen_id": cid,
                "source": "ClassyFire (Wishart Lab)",
                "query_inchikey": inchikey,
                "match_status": result.match_status,
                "chemont_class": result.chemical_class,
                "chemont_direct_parent": result.direct_parent,
                "source_url": result.source_url,
                "error": result.error,
            }
        )
    enriched = pd.concat([df.reset_index(drop=True), pd.DataFrame(records)], axis=1)
    return enriched, provenance


def finalize_carcinogen_class(df: pd.DataFrame, cfg: PipelineConfig) -> pd.DataFrame:
    """Promote ClassyFire results into ``carcinogen_class_normalized``.

    Behaviour:

    - Rows whose source is ``seed_csv`` or ``iarc_entity_type`` are
      left untouched (curated values win).
    - Rows whose source is ``rules`` keep their regex-derived class but
      gain a ``classyfire_validated`` boolean indicating whether the
      ChemOnt taxonomy independently agrees with the assigned bucket.
    - Rows whose source is ``unknown`` are upgraded:

      1. to a project bucket if a ``CHEMONT_CLASS_RULES`` entry matched
         (source → ``classyfire``);
      2. otherwise to the raw ChemOnt class / direct_parent string if
         ClassyFire matched the structure (source → ``chemont_class``);
      3. otherwise they stay ``Unknown``.
    """
    out = df.copy()
    if "carcinogen_class_normalized" not in out.columns:
        return out
    if "carcinogen_class_source" not in out.columns:
        out["carcinogen_class_source"] = "unknown"

    chemont_buckets: list[str] = []
    classyfire_validated: list[bool] = []
    new_classes: list[str] = []
    new_sources: list[str] = []
    for _, row in out.iterrows():
        alt = row.get("chemont_alternative_parents", "") or ""
        if isinstance(alt, str):
            alt_list = [p.strip() for p in alt.split(";") if p.strip()]
        else:
            alt_list = list(alt)
        chemont_record = {
            "kingdom": row.get("chemont_kingdom", ""),
            "superclass": row.get("chemont_superclass", ""),
            "chemical_class": row.get("chemont_class", ""),
            "subclass": row.get("chemont_subclass", ""),
            "direct_parent": row.get("chemont_direct_parent", ""),
            "alternative_parents": alt_list,
        }
        bucket = chemont_to_class_bucket(chemont_record)
        chemont_buckets.append(bucket)

        current_class = str(row.get("carcinogen_class_normalized", "") or "").strip()
        current_source = str(row.get("carcinogen_class_source", "") or "").strip()
        cf_status = str(row.get("classyfire_match_status", "") or "").strip()
        chemont_class = str(row.get("chemont_class", "") or "").strip()
        chemont_direct_parent = str(row.get("chemont_direct_parent", "") or "").strip()

        if current_source in {"seed_csv", "iarc_entity_type"} or (
            cfg.classyfire.upgrade_unknown_only and current_source != "unknown"
        ):
            new_classes.append(current_class)
            new_sources.append(current_source or "unknown")
            classyfire_validated.append(bool(bucket) and bucket == current_class)
            continue

        if bucket:
            new_classes.append(bucket)
            new_sources.append("classyfire")
            classyfire_validated.append(True)
        elif cf_status == "matched" and (chemont_direct_parent or chemont_class):
            # No project bucket matched, but we have a structure-based
            # ChemOnt classification. Surface that as the class so the
            # row is no longer ``Unknown``.
            new_classes.append(chemont_direct_parent or chemont_class)
            new_sources.append("chemont_class")
            classyfire_validated.append(True)
        else:
            new_classes.append(current_class or "Unknown")
            new_sources.append(current_source or "unknown")
            classyfire_validated.append(False)

    out["chemont_class_bucket"] = chemont_buckets
    out["classyfire_validated"] = classyfire_validated
    out["carcinogen_class_normalized"] = new_classes
    out["carcinogen_class_source"] = new_sources
    out["norm_carcinogen_class"] = out["carcinogen_class_normalized"]
    return out


def enrich_local_references(
    df: pd.DataFrame,
    iarc_csv: str | Path | None = None,
    ntp_csv: str | Path | None = None,
    oehha_csv: str | Path | None = None,
    ctd_diseases: str | Path | None = None,
    ctd_genes: str | Path | None = None,
) -> pd.DataFrame:
    if iarc_csv:
        iarc = load_optional_table(iarc_csv)
        df = left_enrich_by_names(
            df,
            iarc,
            seed_name_columns=["display_label", "preferred_query_name", "iarc_agent_name"],
            ref_name_columns=["Agent", "agent", "Substance", "substance", "Name", "name"],
            output_prefix="iarc_ref",
            keep_columns=[c for c in ["Agent", "agent", "Group", "group", "Volume", "volume", "Year", "year", "CAS", "cas"] if c in iarc.columns],
        )
    if ntp_csv:
        ntp = load_optional_table(ntp_csv)
        df = left_enrich_by_names(
            df,
            ntp,
            seed_name_columns=["display_label", "preferred_query_name", "casrn", "iarc_casrn"],
            ref_name_columns=["Name", "name", "Substance", "substance", "CASRN", "CAS", "casrn", "cas"],
            output_prefix="ntp",
            keep_columns=[c for c in ["Name", "Substance", "Listing", "status", "CASRN", "CAS"] if c in ntp.columns],
        )
    if oehha_csv:
        oehha = load_optional_table(oehha_csv)
        df = left_enrich_by_names(
            df,
            oehha,
            seed_name_columns=["display_label", "preferred_query_name", "casrn", "iarc_casrn"],
            ref_name_columns=["Chemical", "chemical", "Name", "name", "CAS No.", "CAS", "casrn"],
            output_prefix="oehha",
            keep_columns=[c for c in ["Chemical", "Name", "CAS No.", "CAS", "Toxicity Endpoint", "Listing Mechanism"] if c in oehha.columns],
        )

    if ctd_diseases:
        df = summarize_ctd_diseases(df, ctd_diseases)
    if ctd_genes:
        df = summarize_ctd_genes(df, ctd_genes)
    return df


def build_dataset(
    seed_df: pd.DataFrame,
    cfg: PipelineConfig,
    cache_dir: str | Path,
    iarc_csv: str | Path | None = None,
    ntp_csv: str | Path | None = None,
    oehha_csv: str | Path | None = None,
    ctd_diseases: str | Path | None = None,
    ctd_genes: str | Path | None = None,
    use_pubchem: bool = True,
    use_classyfire: bool = True,
    classyfire_cache_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    normalized = normalize_seed(seed_df)
    provenance = [provenance_for_seed(row) for _, row in normalized.iterrows()]
    if use_pubchem:
        enriched, pubchem_prov = enrich_pubchem(normalized, cfg, cache_dir=cache_dir)
        provenance.extend(pubchem_prov)
    else:
        enriched = normalized

    if use_classyfire:
        cf_cache = (
            Path(classyfire_cache_dir)
            if classyfire_cache_dir
            else Path(cache_dir).parent / "classyfire"
        )
        enriched, classyfire_prov = enrich_classyfire(enriched, cfg, cache_dir=cf_cache)
        provenance.extend(classyfire_prov)
    else:
        empty = pd.DataFrame(
            [_empty_classyfire_record("disabled") for _ in range(len(enriched))]
        )
        enriched = pd.concat([enriched.reset_index(drop=True), empty], axis=1)

    enriched = finalize_carcinogen_class(enriched, cfg)
    if "norm_carcinogen_class" not in enriched.columns:
        enriched["norm_carcinogen_class"] = enriched.get("carcinogen_class_normalized", "")
    enriched = enrich_local_references(enriched, iarc_csv=iarc_csv, ntp_csv=ntp_csv, oehha_csv=oehha_csv, ctd_diseases=ctd_diseases, ctd_genes=ctd_genes)
    return enriched, provenance
