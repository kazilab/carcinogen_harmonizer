#!/usr/bin/env python
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import date
import json
from pathlib import Path
import sys
from typing import Any

try:
    from carcinogen_harmonizer.sources.pubchem import PubChemClient
except ModuleNotFoundError:
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))
    from carcinogen_harmonizer.sources.pubchem import PubChemClient


DEFAULT_ALIAS_MAP_PATH = (
    Path(__file__).resolve().parents[1]
    / "carcinogen_harmonizer"
    / "data"
    / "alias_map.json"
)
DEFAULT_CACHE_DIR = Path(__file__).resolve().parents[1] / "cache" / "pubchem"
PUBCHEM_PROPERTIES = [
    "MolecularFormula",
    "MolecularWeight",
    "CanonicalSMILES",
    "IsomericSMILES",
    "InChI",
    "InChIKey",
    "IUPACName",
]


def load_alias_map(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def format_alias_map_json(data: dict[str, Any]) -> str:
    return json.dumps(data, ensure_ascii=False, indent=2) + "\n"


def refresh_alias_data(
    data: dict[str, Any],
    client: PubChemClient,
    *,
    retrieved_on: str,
    max_synonyms: int = 25,
) -> tuple[dict[str, Any], list[str]]:
    refreshed = deepcopy(data)
    metadata = refreshed.setdefault("metadata", {})
    metadata.update(
        {
            "generated_from": "PubChem PUG REST compound name lookups",
            "primary_source": "PubChem",
            "primary_source_url": "https://pubchem.ncbi.nlm.nih.gov/",
            "retrieved_on": retrieved_on,
            "source_policy": (
                "Entries must be backed by a recognized external chemical database record; "
                "local seed files are not used as authority for this map."
            ),
        }
    )

    aliases = refreshed.get("aliases")
    if not isinstance(aliases, dict):
        raise ValueError("alias_map.json must contain an object at key 'aliases'")

    errors: list[str] = []
    sorted_aliases: dict[str, Any] = {}
    for alias in sorted(aliases):
        raw_entry = aliases[alias]
        if isinstance(raw_entry, str):
            entry: dict[str, Any] = {"preferred_query_name": raw_entry}
        elif isinstance(raw_entry, dict):
            entry = dict(raw_entry)
        else:
            errors.append(f"{alias}: entry must be an object or string")
            sorted_aliases[alias] = raw_entry
            continue

        preferred_query_name = str(entry.get("preferred_query_name", "")).strip()
        if not preferred_query_name:
            errors.append(f"{alias}: missing preferred_query_name")
            sorted_aliases[alias] = entry
            continue

        result = client.enrich_name(
            preferred_query_name,
            properties=PUBCHEM_PROPERTIES,
            max_synonyms=max_synonyms,
        )
        if result.match_status != "matched" or not result.pubchem_cid:
            detail = f"{result.match_status}"
            if result.error:
                detail = f"{detail}: {result.error}"
            errors.append(f"{alias}: PubChem lookup for {preferred_query_name!r} failed ({detail})")
            sorted_aliases[alias] = entry
            continue

        cid = str(result.pubchem_cid)
        entry["source"] = {
            "database": "PubChem",
            "record_type": "Compound",
            "cid": cid,
            "url": result.source_url or f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        }
        sorted_aliases[alias] = entry

    refreshed["aliases"] = sorted_aliases
    return refreshed, errors


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Refresh carcinogen alias-map PubChem CID/URL provenance."
    )
    parser.add_argument(
        "--alias-map",
        type=Path,
        default=DEFAULT_ALIAS_MAP_PATH,
        help=f"Alias map JSON to refresh (default: {DEFAULT_ALIAS_MAP_PATH})",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help=f"PubChem cache directory (default: {DEFAULT_CACHE_DIR})",
    )
    parser.add_argument("--timeout", type=int, default=20, help="PubChem request timeout in seconds")
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=0.2,
        help="Seconds to sleep between PubChem requests",
    )
    parser.add_argument(
        "--max-synonyms",
        type=int,
        default=25,
        help="Maximum synonyms to request from PubChem while refreshing",
    )
    parser.add_argument(
        "--retrieved-on",
        default="",
        help="Override metadata.retrieved_on, in YYYY-MM-DD form",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate that the alias map is current without writing changes",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    data = load_alias_map(args.alias_map)
    existing_date = str(data.get("metadata", {}).get("retrieved_on", "")).strip()
    retrieved_on = args.retrieved_on or (existing_date if args.check else date.today().isoformat())

    client = PubChemClient(
        cache_dir=args.cache_dir,
        timeout=args.timeout,
        rate_limit_seconds=args.rate_limit,
    )
    refreshed, errors = refresh_alias_data(
        data,
        client,
        retrieved_on=retrieved_on,
        max_synonyms=args.max_synonyms,
    )
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1

    old_text = format_alias_map_json(data)
    new_text = format_alias_map_json(refreshed)
    if args.check:
        if old_text != new_text:
            print(f"{args.alias_map} is stale; run tools/refresh_alias_map.py", file=sys.stderr)
            return 1
        print(f"{args.alias_map} is current")
        return 0

    args.alias_map.write_text(new_text, encoding="utf-8")
    print(f"Refreshed {len(refreshed['aliases'])} aliases in {args.alias_map}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
