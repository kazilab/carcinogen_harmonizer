from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

@dataclass
class PubChemConfig:
    enabled: bool = True
    rate_limit_seconds: float = 0.2
    timeout_seconds: int = 20
    max_synonyms: int = 25
    properties: list[str] = field(default_factory=lambda: [
        "MolecularFormula",
        "MolecularWeight",
        "CanonicalSMILES",
        "IsomericSMILES",
        "InChI",
        "InChIKey",
        "IUPACName",
        "XLogP",
        "TPSA",
    ])


@dataclass
class ClassyFireConfig:
    """Settings for the ClassyFire (ChemOnt) chemical taxonomy lookups.

    The Wishart Lab server returns HTTP 429 if pinged too fast. Defaults
    are tuned to keep a one-shot run for ~600 IARC chemical-like agents
    reliable: ~1.5s between requests with exponential backoff on 429.
    All settled responses (``matched`` / ``no_classyfire_match``) are
    cached on disk; transient ``api_error`` results are not cached so
    they are retried on the next run.
    """

    enabled: bool = True
    rate_limit_seconds: float = 1.5
    timeout_seconds: int = 30
    max_retries: int = 4
    backoff_base_seconds: float = 4.0
    upgrade_unknown_only: bool = True


@dataclass
class QualityConfig:
    required_identifier_any: list[str] = field(default_factory=lambda: [
        "pubchem_cid",
        "inchikey",
        "casrn",
        "iarc_casrn",
        "dtxsid",
        "chebi_id",
    ])
    controlled_iarc_groups: list[str] = field(default_factory=lambda: [
        "Group 1",
        "Group 2A",
        "Group 2B",
        "Group 3",
        "Not classified",
    ])


@dataclass
class PipelineConfig:
    pubchem: PubChemConfig = field(default_factory=PubChemConfig)
    classyfire: ClassyFireConfig = field(default_factory=ClassyFireConfig)
    quality: QualityConfig = field(default_factory=QualityConfig)
    raw: dict[str, Any] = field(default_factory=dict)


def load_config(path: str | Path | None = None) -> PipelineConfig:
    data: dict[str, Any] = {}
    if path:
        p = Path(path)
        if p.exists():
            try:
                import yaml
            except ModuleNotFoundError as exc:
                raise RuntimeError(
                    "Missing dependency 'pyyaml'. Install with: pip install pyyaml"
                ) from exc
            data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pubchem = PubChemConfig(**(data.get("pubchem") or {}))
    classyfire = ClassyFireConfig(**(data.get("classyfire") or {}))
    quality = QualityConfig(**(data.get("quality") or {}))
    return PipelineConfig(pubchem=pubchem, classyfire=classyfire, quality=quality, raw=data)
