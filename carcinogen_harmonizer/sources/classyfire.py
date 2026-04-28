"""ClassyFire (ChemOnt) structure-based chemical taxonomy lookups.

ClassyFire is a free Wishart-Lab service that returns a hierarchical
chemical class for any chemical structure. It is used here as a
"validated chemical class" source for rows where the seed CSV did
not provide one and the simple regex fallback in
``normalize.infer_carcinogen_class`` returned ``Unknown``.

Documentation: http://classyfire.wishartlab.com/

We only call the ``GET /entities/{inchikey}.json`` endpoint, which
returns pre-computed classifications for any structure that has
already been classified. This is sufficient for nearly every IARC
agent and Group 1/2A/2B chemical, which are well-known compounds
already present in the public ClassyFire cache.

If a row has no InChIKey (e.g. mixtures, exposures, biological
agents, radiation), the call is skipped.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import quote
import json
import time

import requests


CLASSYFIRE_BASE = "http://classyfire.wishartlab.com"
CHEMONT_BROWSER = "http://classyfire.wishartlab.com/tax_nodes"


@dataclass
class ClassyFireResult:
    inchikey: str = ""
    smiles: str = ""
    kingdom: str = ""
    kingdom_id: str = ""
    superclass: str = ""
    superclass_id: str = ""
    chemical_class: str = ""
    chemical_class_id: str = ""
    subclass: str = ""
    subclass_id: str = ""
    direct_parent: str = ""
    direct_parent_id: str = ""
    alternative_parents: list[str] = field(default_factory=list)
    description: str = ""
    match_status: str = "not_queried"
    error: str = ""
    source_url: str = ""

    def to_record(self) -> dict[str, Any]:
        return asdict(self)


class ClassyFireClient:
    """Thin, file-cached HTTP client for ClassyFire entity lookups.

    The Wishart Lab service applies a fairly aggressive rate limit. The
    client therefore:

    - sleeps ``rate_limit_seconds`` between requests
    - retries on HTTP 429 (Too Many Requests) with exponential backoff
    - **does not cache** transient API errors so they can be retried on
      the next pipeline run

    Successful matches and confirmed 404s (no ChemOnt entry) are
    cached on disk for fast re-runs.
    """

    def __init__(
        self,
        cache_dir: str | Path = "cache/classyfire",
        timeout: int = 30,
        rate_limit_seconds: float = 1.5,
        max_retries: int = 4,
        backoff_base_seconds: float = 4.0,
    ) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.rate_limit_seconds = rate_limit_seconds
        self.max_retries = max_retries
        self.backoff_base_seconds = backoff_base_seconds
        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": "carcinogen-pipeline/0.3",
                "Accept": "application/json",
            }
        )

    def lookup_inchikey(self, inchikey: str) -> ClassyFireResult:
        inchikey = str(inchikey or "").strip()
        if not inchikey:
            return ClassyFireResult(match_status="empty_query")
        cache_path = self._cache_path(inchikey)
        if cache_path.exists():
            try:
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
                cached_result = ClassyFireResult(**cached)
                # Treat previously cached transient errors as not-cached
                # so we retry them on this run.
                if cached_result.match_status != "api_error":
                    return cached_result
            except Exception:
                pass
        result = ClassyFireResult(inchikey=inchikey)
        url = f"{CLASSYFIRE_BASE}/entities/{quote(inchikey)}.json"
        result.source_url = url
        try:
            self._fetch_with_retry(url, result)
        except requests.HTTPError as e:
            result.match_status = "api_error"
            status = e.response.status_code if e.response is not None else ""
            result.error = f"HTTPError: {status} {e}"
        except Exception as e:
            result.match_status = "api_error"
            result.error = f"{type(e).__name__}: {e}"
        # Only cache settled outcomes (matched / no_classyfire_match).
        if result.match_status in {"matched", "no_classyfire_match"}:
            self._write_cache(cache_path, result)
        return result

    def _fetch_with_retry(self, url: str, result: ClassyFireResult) -> None:
        attempt = 0
        while True:
            time.sleep(self.rate_limit_seconds)
            resp = self.session.get(url, timeout=self.timeout)
            if resp.status_code == 404:
                result.match_status = "no_classyfire_match"
                return
            if resp.status_code == 429 and attempt < self.max_retries:
                retry_after = self._retry_after_seconds(resp, attempt)
                time.sleep(retry_after)
                attempt += 1
                continue
            resp.raise_for_status()
            data = resp.json()
            self._populate(result, data)
            result.match_status = "matched"
            return

    def _retry_after_seconds(self, response: requests.Response, attempt: int) -> float:
        header = response.headers.get("Retry-After")
        if header:
            try:
                return max(float(header), self.backoff_base_seconds)
            except ValueError:
                pass
        return self.backoff_base_seconds * (2**attempt)

    @staticmethod
    def _populate(result: ClassyFireResult, data: dict[str, Any]) -> None:
        result.smiles = str(data.get("smiles") or "")
        node_map = [
            ("kingdom", "kingdom", "kingdom_id"),
            ("superclass", "superclass", "superclass_id"),
            ("class", "chemical_class", "chemical_class_id"),
            ("subclass", "subclass", "subclass_id"),
            ("direct_parent", "direct_parent", "direct_parent_id"),
        ]
        for json_key, name_attr, id_attr in node_map:
            node = data.get(json_key) or {}
            if isinstance(node, dict):
                setattr(result, name_attr, str(node.get("name", "") or ""))
                setattr(result, id_attr, str(node.get("chemont_id", "") or ""))
        alt_parents = data.get("alternative_parents") or []
        result.alternative_parents = [
            str((entry or {}).get("name", ""))
            for entry in alt_parents
            if isinstance(entry, dict) and entry.get("name")
        ]
        result.description = str(data.get("description") or "")

    def _cache_path(self, inchikey: str) -> Path:
        safe = str(inchikey).replace("/", "_")
        return self.cache_dir / f"{safe}.json"

    def _write_cache(self, path: Path, result: ClassyFireResult) -> None:
        path.write_text(
            json.dumps(result.to_record(), ensure_ascii=False, indent=2, sort_keys=True),
            encoding="utf-8",
        )


# ---------------------------------------------------------------------------
# Mapping helpers: ChemOnt class → project class buckets used by
# ``carcinogen_class_normalized``. Order matters; the first matching rule
# wins. Patterns are case-insensitive substring matches against the
# ChemOnt direct_parent / class / superclass strings, in that order.
# ---------------------------------------------------------------------------

CHEMONT_CLASS_RULES: list[tuple[str, str]] = [
    ("polychlorinated dibenzodioxin", "Dioxin"),
    ("dibenzodioxin", "Dioxin"),
    ("dibenzofuran", "Dioxin"),
    ("polychlorinated biphenyl", "PCB"),
    ("chlorinated biphenyl", "PCB"),
    ("polycyclic aromatic hydrocarbon", "PAH"),
    ("aromatic homopolycyclic", "PAH"),
    ("aflatoxin", "Mycotoxin"),
    ("mycotoxin", "Mycotoxin"),
    ("nitrosamine", "Nitrosamine"),
    ("n-nitroso", "Nitrosamine"),
    ("aldehyde", "Aldehyde"),
    ("aniline", "Aromatic_Amine"),
    ("aromatic amine", "Aromatic_Amine"),
    ("aminobiphenyl", "Aromatic_Amine"),
    ("benzidine", "Aromatic_Amine"),
    ("estrogen", "Estrogen"),
    ("estradiol", "Estrogen"),
    ("androgen", "Androgen"),
    ("organochlor", "Organochlorine"),
    ("organochloride", "Organochlorine"),
    ("organohalide", "Organochlorine"),
    ("imidazopyridine", "HCA"),
    ("imidazoquinoline", "HCA"),
    ("imidazoquinoxaline", "HCA"),
    ("alkylating", "Alkylating"),
    ("nitrosourea", "Alkylating"),
    ("epoxide", "Alkylating"),
    ("aziridine", "Alkylating"),
    ("ethylene oxide", "Alkylating"),
    ("nitrogen mustard", "Alkylating"),
    ("sulfur mustard", "Alkylating"),
    ("haloalkane", "Chlorinated_Solvent"),
    ("chlorinated solvent", "Chlorinated_Solvent"),
    ("benzene", "Solvent"),
    ("vinyl halide", "Solvent"),
    ("alcohol", "Alcohol"),
    ("organometal", "Heavy_Metal"),
    ("metalloid", "Heavy_Metal"),
]


def chemont_to_class_bucket(result: ClassyFireResult | dict[str, Any] | None) -> str:
    """Map a ClassyFire result to one of the project's class buckets.

    Searches ``direct_parent``, ``subclass``, ``class``, ``superclass``,
    and the full ``alternative_parents`` list. ChemOnt frequently
    categorizes molecules by skeleton (e.g. "Biphenyls and derivatives")
    while the functional-group classifications appear in the
    ``alternative_parents`` list (e.g. "Primary aromatic amines"), so
    both must be considered.

    Returns an empty string if no rule matched.
    """
    if result is None:
        return ""
    if isinstance(result, ClassyFireResult):
        candidates = [
            result.direct_parent,
            result.subclass,
            result.chemical_class,
            result.superclass,
            *(result.alternative_parents or []),
        ]
    else:
        alt_parents = result.get("alternative_parents") or []
        if isinstance(alt_parents, str):
            alt_parents = [p.strip() for p in alt_parents.split(";") if p.strip()]
        candidates = [
            str(result.get("direct_parent", "") or ""),
            str(result.get("subclass", "") or ""),
            str(result.get("chemical_class", "") or result.get("class", "") or ""),
            str(result.get("superclass", "") or ""),
            *[str(p or "") for p in alt_parents],
        ]
    haystack = " | ".join(c.lower() for c in candidates if c)
    if not haystack:
        return ""
    for needle, bucket in CHEMONT_CLASS_RULES:
        if needle in haystack:
            return bucket
    return ""


def chemont_reference_url(chemont_id: str) -> str:
    chemont_id = str(chemont_id or "").strip()
    if not chemont_id:
        return ""
    return f"{CHEMONT_BROWSER}/{chemont_id}"
