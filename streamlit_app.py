from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from carcinogen_harmonizer.app_utils import (
    DEFAULT_RESULT_COLUMNS,
    filter_dataset,
    iarc_group_summary,
    read_static_dataset,
    sorted_unique_values,
)
from carcinogen_harmonizer.pipeline import PipelineRunOptions, run_pipeline
from carcinogen_harmonizer.sources.iarc import DEFAULT_IARC_URL


ROOT = Path(__file__).resolve().parent
DEFAULT_STATIC_DATA = ROOT / "outputs" / "carcinogens_enriched.csv"
DEFAULT_SEED = ROOT / "data" / "seed.csv"
DEFAULT_OUT = ROOT / "outputs"
DEFAULT_CACHE = ROOT / "cache"
DEFAULT_CONFIG = ROOT / "config.yaml"
IARC_RUN_GROUPS = ["Group 1", "Group 2A", "Group 2B", "Group 3"]


st.set_page_config(
    page_title="Carcinogen Harmonizer",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_data(show_spinner=False)
def _load_static_dataset(path: str, mtime: float) -> pd.DataFrame:
    del mtime
    return read_static_dataset(path)


def _path_mtime(path: str) -> float:
    p = Path(path)
    return p.stat().st_mtime if p.exists() else 0.0


def _optional_path(value: str) -> str | None:
    value = value.strip()
    return value or None


def _selected_display_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in DEFAULT_RESULT_COLUMNS if col in df.columns]


def render_static_search() -> None:
    st.header("Search Static Data")
    static_path = st.text_input("Static dataset CSV", value=str(DEFAULT_STATIC_DATA))
    if st.button("Reload static data"):
        st.cache_data.clear()

    try:
        df = _load_static_dataset(static_path, _path_mtime(static_path))
    except FileNotFoundError as exc:
        st.warning(str(exc))
        return

    metric_cols = st.columns(4)
    metric_cols[0].metric("Rows", f"{len(df):,}")
    metric_cols[1].metric("IARC groups", f"{len(sorted_unique_values(df, 'iarc_group_normalized')):,}")
    metric_cols[2].metric("Classes", f"{len(sorted_unique_values(df, 'carcinogen_class_normalized')):,}")
    if "needs_manual_review" in df.columns:
        needs_review = df["needs_manual_review"].astype(str).str.lower().isin(["true", "1", "yes"]).sum()
        metric_cols[3].metric("Needs review", f"{int(needs_review):,}")
    else:
        metric_cols[3].metric("Needs review", "n/a")

    group_summary = iarc_group_summary(df)
    with st.expander("IARC group list", expanded=True):
        st.dataframe(group_summary, use_container_width=True, hide_index=True)

    filters = st.container(border=True)
    query = filters.text_input(
        "Search",
        placeholder="Name, alias, CAS, PubChem CID, InChIKey, IARC agent...",
    )
    left, mid, right = filters.columns(3)
    groups = left.multiselect("IARC groups", sorted_unique_values(df, "iarc_group_normalized"))
    classes = mid.multiselect("Carcinogen classes", sorted_unique_values(df, "carcinogen_class_normalized"))
    entity_types = right.multiselect("Entity types", sorted_unique_values(df, "agent_entity_type"))
    manual_review = filters.segmented_control(
        "Review status",
        ["All", "Needs review", "No review flag"],
        default="All",
    )

    result = filter_dataset(
        df,
        query,
        groups=groups,
        classes=classes,
        entity_types=entity_types,
        manual_review=manual_review,
    )
    st.subheader(f"Results: {len(result):,}")
    display_cols = _selected_display_columns(result)
    st.dataframe(result[display_cols] if display_cols else result, use_container_width=True, hide_index=True)

    csv_bytes = result.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download filtered CSV",
        data=csv_bytes,
        file_name="carcinogen_search_results.csv",
        mime="text/csv",
        disabled=result.empty,
    )

    if not result.empty:
        label_col = "display_label" if "display_label" in result.columns else result.columns[0]
        id_col = "carcinogen_id" if "carcinogen_id" in result.columns else result.columns[0]
        options = result.index.tolist()
        selected = st.selectbox(
            "Record details",
            options=options,
            format_func=lambda idx: f"{result.loc[idx, id_col]} | {result.loc[idx, label_col]}",
        )
        detail = result.loc[selected].to_dict()
        details_left, details_right = st.columns([2, 1])
        detail_fields = [
            "carcinogen_id",
            "display_label",
            "preferred_query_name",
            "pubchem_lookup_name",
            "pubchem_synonyms",
            "iarc_group_normalized",
            "carcinogen_class_normalized",
            "agent_entity_type",
            "pubchem_cid",
            "casrn",
            "inchikey",
            "iupac_name",
            "confidence_score",
            "needs_manual_review",
        ]
        details_left.json({k: detail.get(k, "") for k in detail_fields if k in detail})
        links: list[tuple[str, str]] = []
        if detail.get("pubchem_source_url"):
            links.append(("PubChem", str(detail["pubchem_source_url"])))
        if detail.get("classyfire_reference_url"):
            links.append(("ClassyFire/ChemOnt", str(detail["classyfire_reference_url"])))
        if detail.get("iarc_source_url"):
            links.append(("IARC source", str(detail["iarc_source_url"])))
        if links:
            details_right.write("External records")
            for label, url in links:
                details_right.link_button(label, url)


def render_pipeline_runner() -> None:
    st.header("Run Pipeline")
    run_mode = st.radio(
        "Run scope",
        [
            "Normalize seed only",
            "Build without online enrichment",
            "Build with PubChem",
            "Full build",
        ],
        horizontal=True,
        index=3,
    )

    basic = st.container(border=True)
    seed = basic.text_input("Seed CSV", value=str(DEFAULT_SEED))
    out = basic.text_input("Output directory", value=str(DEFAULT_OUT))
    cache = basic.text_input("Cache directory", value=str(DEFAULT_CACHE))
    config = basic.text_input("Config YAML", value=str(DEFAULT_CONFIG))

    iarc = st.container(border=True)
    expand_iarc = iarc.checkbox("Expand from IARC classifications", value=True)
    iarc_source = iarc.text_input("IARC source", value=DEFAULT_IARC_URL)
    iarc_groups = iarc.multiselect("IARC groups to include", IARC_RUN_GROUPS, default=IARC_RUN_GROUPS[:3])
    include_group3 = "Group 3" in iarc_groups

    optional = st.expander("Optional local reference tables")
    iarc_csv = optional.text_input("Local IARC reference CSV/XLSX/TSV", value="")
    ntp_csv = optional.text_input("NTP reference CSV", value="")
    oehha_csv = optional.text_input("OEHHA reference CSV", value="")
    ctd_diseases = optional.text_input("CTD diseases TSV/TSV.GZ", value="")
    ctd_genes = optional.text_input("CTD genes TSV/TSV.GZ", value="")

    use_pubchem = run_mode in {"Build with PubChem", "Full build"}
    use_classyfire = run_mode == "Full build"
    normalize_only = run_mode == "Normalize seed only"

    st.caption(
        f"Effective run: PubChem={'on' if use_pubchem else 'off'}, "
        f"ClassyFire={'on' if use_classyfire else 'off'}, "
        f"IARC expansion={'on' if expand_iarc else 'off'}"
    )

    if st.button("Run pipeline", type="primary"):
        with st.status("Running pipeline", expanded=True) as status:
            st.write("Reading seed and configuration")
            options = PipelineRunOptions(
                seed=seed,
                out=out,
                config=config,
                cache=cache,
                use_pubchem=use_pubchem,
                use_classyfire=use_classyfire,
                iarc_csv=_optional_path(iarc_csv),
                ntp_csv=_optional_path(ntp_csv),
                oehha_csv=_optional_path(oehha_csv),
                ctd_diseases=_optional_path(ctd_diseases),
                ctd_genes=_optional_path(ctd_genes),
                expand_iarc=expand_iarc,
                iarc_source=iarc_source,
                include_group3=include_group3,
                iarc_groups=iarc_groups,
                normalize_only=normalize_only,
            )
            result = run_pipeline(options)
            status.update(label="Pipeline complete", state="complete")

        summary_cols = st.columns(4)
        summary_cols[0].metric("Rows", f"{len(result.df):,}")
        summary_cols[1].metric("IARC records loaded", f"{result.iarc_records_loaded:,}")
        summary_cols[2].metric("Manual review", f"{result.report.get('manual_review_count', 0):,}")
        summary_cols[3].metric("Output", str(result.output_dir))
        st.dataframe(result.df[_selected_display_columns(result.df)], use_container_width=True, hide_index=True)


def main() -> None:
    st.title("Carcinogen Harmonizer")
    st.caption("Search generated carcinogen data or run the harmonization pipeline from the same UI.")
    tab_search, tab_run = st.tabs(["Static Search", "Run Pipeline"])
    with tab_search:
        render_static_search()
    with tab_run:
        render_pipeline_runner()


if __name__ == "__main__":
    main()
