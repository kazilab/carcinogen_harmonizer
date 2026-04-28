from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from .config import PipelineConfig, load_config
from .enrich import build_dataset
from .io import (
    read_seed_csv,
    write_csv,
    write_dataframe_json,
    write_json,
    write_jsonl,
)
from .quality import add_quality_flags, graph_tables, make_qa_report
from .sources.iarc import DEFAULT_IARC_URL, load_iarc_agents, merge_seed_with_iarc_expansion


@dataclass
class PipelineRunOptions:
    seed: str | Path
    out: str | Path = "outputs"
    config: str | Path = "config.yaml"
    cache: str | Path = "cache"
    use_pubchem: bool = True
    use_classyfire: bool = True
    classyfire_cache: str | Path | None = None
    iarc_csv: str | Path | None = None
    ntp_csv: str | Path | None = None
    oehha_csv: str | Path | None = None
    ctd_diseases: str | Path | None = None
    ctd_genes: str | Path | None = None
    expand_iarc: bool = False
    iarc_source: str | Path = DEFAULT_IARC_URL
    include_group3: bool = False
    iarc_groups: list[str] = field(default_factory=lambda: ["Group 1", "Group 2A", "Group 2B"])
    output_format: str = "both"
    write_outputs: bool = True
    normalize_only: bool = False


@dataclass
class PipelineRunResult:
    df: pd.DataFrame
    provenance: list[dict[str, Any]]
    report: dict[str, Any]
    graph: dict[str, pd.DataFrame]
    iarc_skipped: pd.DataFrame | None
    original_seed_count: int
    iarc_records_loaded: int
    output_dir: Path

    @property
    def appended_iarc_count(self) -> int:
        return max(0, len(self.df) - self.original_seed_count)


def run_pipeline(options: PipelineRunOptions, cfg: PipelineConfig | None = None) -> PipelineRunResult:
    out_dir = Path(options.out)
    cache_dir = Path(options.cache)
    if options.write_outputs:
        out_dir.mkdir(parents=True, exist_ok=True)
        cache_dir.mkdir(parents=True, exist_ok=True)

    cfg = cfg or load_config(options.config)
    seed = read_seed_csv(options.seed)
    original_seed_count = len(seed)
    iarc_skipped = None
    iarc_records_loaded = 0

    if options.expand_iarc:
        iarc_agents = load_iarc_agents(
            options.iarc_source,
            include_groups=options.iarc_groups,
            include_group3=options.include_group3,
        )
        iarc_records_loaded = len(iarc_agents)
        seed, iarc_skipped = merge_seed_with_iarc_expansion(seed, iarc_agents)

    if options.normalize_only:
        from .normalize import normalize_seed, provenance_for_seed

        df = normalize_seed(seed)
        provenance = [provenance_for_seed(row) for _, row in df.iterrows()]
    else:
        df, provenance = build_dataset(
            seed,
            cfg=cfg,
            cache_dir=cache_dir / "pubchem",
            iarc_csv=options.iarc_csv,
            ntp_csv=options.ntp_csv,
            oehha_csv=options.oehha_csv,
            ctd_diseases=options.ctd_diseases,
            ctd_genes=options.ctd_genes,
            use_pubchem=options.use_pubchem,
            use_classyfire=options.use_classyfire,
            classyfire_cache_dir=Path(options.classyfire_cache)
            if options.classyfire_cache
            else cache_dir / "classyfire",
        )

    df = add_quality_flags(df, cfg)
    report = make_qa_report(df)
    graph = graph_tables(df)
    result = PipelineRunResult(
        df=df,
        provenance=provenance,
        report=report,
        graph=graph,
        iarc_skipped=iarc_skipped,
        original_seed_count=original_seed_count,
        iarc_records_loaded=iarc_records_loaded,
        output_dir=out_dir,
    )

    if options.write_outputs:
        write_pipeline_outputs(result, options.output_format)

    return result


def write_pipeline_outputs(result: PipelineRunResult, output_format: str = "both") -> None:
    out_dir = result.output_dir
    df = result.df
    write_csv(df, out_dir / "carcinogens_enriched.csv")
    if output_format in {"json", "both"}:
        write_dataframe_json(df, out_dir / "carcinogens_enriched.json")
    if "needs_manual_review" in df.columns:
        write_csv(df.loc[df["needs_manual_review"]].copy(), out_dir / "manual_review_queue.csv")
    if result.iarc_skipped is not None:
        write_csv(result.iarc_skipped, out_dir / "iarc_expansion_skipped_duplicates.csv")
    write_jsonl(result.provenance, out_dir / "source_provenance.jsonl")
    write_json(result.report, out_dir / "qa_report.json")
    for name, table in result.graph.items():
        write_csv(table, out_dir / f"{name}.csv")
        if output_format in {"json", "both"}:
            write_dataframe_json(table, out_dir / f"{name}.json")
