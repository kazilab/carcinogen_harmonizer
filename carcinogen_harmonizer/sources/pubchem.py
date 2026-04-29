from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote
import hashlib
import json
import re
import time

import requests

try:
    from ..normalize import clean_text
except ImportError:
    from normalize import clean_text


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CAS_RE = re.compile(r"\b[1-9][0-9]{1,6}-\d{2}-\d\b")


@dataclass
class PubChemResult:
    query: str
    pubchem_cid: str = ""
    molecular_formula: str = ""
    molecular_weight: str = ""
    canonical_smiles: str = ""
    isomeric_smiles: str = ""
    inchi: str = ""
    inchikey: str = ""
    iupac_name: str = ""
    xlogp: str = ""
    tpsa: str = ""
    synonyms: list[str] | None = None
    casrn: str = ""
    match_status: str = "not_queried"
    error: str = ""
    source_url: str = ""

    def to_record(self) -> dict[str, Any]:
        data = asdict(self)
        data["synonyms"] = self.synonyms or []
        return data


class PubChemClient:
    def __init__(self, cache_dir: str | Path = "cache/pubchem", timeout: int = 20, rate_limit_seconds: float = 0.2):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.rate_limit_seconds = rate_limit_seconds
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "carcinogen-pipeline/0.1"})

    def enrich_name(self, name: str, properties: list[str], max_synonyms: int = 25) -> PubChemResult:
        name = clean_text(name)
        if not name:
            return PubChemResult(query=name, match_status="empty_query")
        cache_path = self._cache_path(name)
        if cache_path.exists():
            try:
                return PubChemResult(**json.loads(cache_path.read_text(encoding="utf-8")))
            except Exception:
                pass

        result = PubChemResult(query=name)
        try:
            cid = self._name_to_cid(name)
            if not cid:
                result.match_status = "no_pubchem_match"
                self._write_cache(cache_path, result)
                return result
            result.pubchem_cid = str(cid)
            result.source_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
            props = self._compound_properties(cid, properties)
            result.molecular_formula = str(props.get("MolecularFormula", "") or "")
            result.molecular_weight = str(props.get("MolecularWeight", "") or "")
            result.canonical_smiles = str(props.get("CanonicalSMILES", "") or "")
            result.isomeric_smiles = str(props.get("IsomericSMILES", "") or "")
            result.inchi = str(props.get("InChI", "") or "")
            result.inchikey = str(props.get("InChIKey", "") or "")
            result.iupac_name = str(props.get("IUPACName", "") or "")
            result.xlogp = str(props.get("XLogP", "") or "")
            result.tpsa = str(props.get("TPSA", "") or "")
            synonyms = self._compound_synonyms(cid)[:max_synonyms]
            result.synonyms = synonyms
            result.casrn = pick_casrn(synonyms)
            result.match_status = "matched"
        except requests.HTTPError as e:
            status = e.response.status_code if e.response else ""
            if status == 404 and not result.pubchem_cid:
                result.match_status = "no_pubchem_match"
            else:
                result.match_status = "api_error"
                result.error = f"HTTPError: {status} {str(e)}"
        except Exception as e:
            result.match_status = "api_error"
            result.error = f"{type(e).__name__}: {e}"
        if result.match_status in {"matched", "no_pubchem_match"}:
            self._write_cache(cache_path, result)
        return result

    def _name_to_cid(self, name: str) -> str:
        url = f"{PUBCHEM_BASE}/compound/name/{quote(name)}/cids/JSON"
        try:
            data = self._get_json(url)
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                return ""
            raise
        ids = data.get("IdentifierList", {}).get("CID", [])
        return str(ids[0]) if ids else ""

    def _compound_properties(self, cid: str, properties: list[str]) -> dict[str, Any]:
        prop_string = ",".join(properties)
        url = f"{PUBCHEM_BASE}/compound/cid/{quote(str(cid))}/property/{prop_string}/JSON"
        data = self._get_json(url)
        props = data.get("PropertyTable", {}).get("Properties", [])
        return props[0] if props else {}

    def _compound_synonyms(self, cid: str) -> list[str]:
        url = f"{PUBCHEM_BASE}/compound/cid/{quote(str(cid))}/synonyms/JSON"
        data = self._get_json(url)
        infos = data.get("InformationList", {}).get("Information", [])
        if not infos:
            return []
        return [str(x) for x in infos[0].get("Synonym", [])]

    def _get_json(self, url: str) -> dict[str, Any]:
        time.sleep(self.rate_limit_seconds)
        resp = self.session.get(url, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    def _cache_path(self, name: str) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "empty"
        if len(safe) > 120:
            digest = hashlib.sha1(name.encode("utf-8")).hexdigest()[:12]
            safe = f"{safe[:107]}_{digest}"
        return self.cache_dir / f"{safe}.json"

    def _write_cache(self, path: Path, result: PubChemResult) -> None:
        path.write_text(json.dumps(result.to_record(), ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def pick_casrn(synonyms: list[str] | None) -> str:
    if not synonyms:
        return ""
    for syn in synonyms:
        match = CAS_RE.search(str(syn))
        if match:
            return match.group(0)
    return ""
