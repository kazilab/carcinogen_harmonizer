"""Public package interface for carcinogen-harmonizer."""
from __future__ import annotations

__version__ = "0.4.0"

from .config import load_config
from .enrich import build_dataset, enrich_classyfire, finalize_carcinogen_class
from .io import read_seed_csv, write_csv, write_dataframe_json, write_json, write_jsonl
from .normalize import clean_text, normalize_seed, provenance_for_seed
from .pipeline import PipelineRunOptions, PipelineRunResult, run_pipeline
from .quality import add_quality_flags, graph_tables, make_qa_report
from .sources.classyfire import ClassyFireClient, chemont_reference_url, chemont_to_class_bucket
from .sources.iarc import DEFAULT_IARC_URL, load_iarc_agents, merge_seed_with_iarc_expansion

__all__ = [
    "__version__",
    "DEFAULT_IARC_URL",
    "ClassyFireClient",
    "PipelineRunOptions",
    "PipelineRunResult",
    "add_quality_flags",
    "build_dataset",
    "clean_text",
    "chemont_reference_url",
    "chemont_to_class_bucket",
    "enrich_classyfire",
    "finalize_carcinogen_class",
    "graph_tables",
    "load_config",
    "load_iarc_agents",
    "make_qa_report",
    "merge_seed_with_iarc_expansion",
    "normalize_seed",
    "provenance_for_seed",
    "read_seed_csv",
    "run_pipeline",
    "write_csv",
    "write_dataframe_json",
    "write_json",
    "write_jsonl",
]
