from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace

import pandas as pd
import requests

from carcinogen_harmonizer.config import PipelineConfig
from carcinogen_harmonizer.enrich import enrich_pubchem
from carcinogen_harmonizer.sources import pubchem
from carcinogen_harmonizer.sources.pubchem import PubChemClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.refresh_alias_map import refresh_alias_data  # noqa: E402


class _Response404:
    status_code = 404

    def raise_for_status(self) -> None:
        raise requests.HTTPError(response=self)


class _Response500:
    status_code = 500

    def raise_for_status(self) -> None:
        raise requests.HTTPError(response=self)


def test_pubchem_404_is_reported_as_no_match(monkeypatch, tmp_path):
    client = PubChemClient(cache_dir=tmp_path)

    def fake_get_json(self, url):  # noqa: ANN001
        raise requests.HTTPError(response=_Response404())

    monkeypatch.setattr(PubChemClient, "_get_json", fake_get_json)
    result = client.enrich_name("Nonexistent compound", properties=["InChIKey"])
    assert result.match_status == "no_pubchem_match"
    assert result.error == ""


def test_pubchem_server_failure_remains_api_error(monkeypatch, tmp_path):
    client = PubChemClient(cache_dir=tmp_path)

    def fake_get_json(self, url):  # noqa: ANN001
        raise requests.HTTPError(response=_Response500())

    monkeypatch.setattr(PubChemClient, "_get_json", fake_get_json)
    result = client.enrich_name("Nonexistent compound", properties=["InChIKey"])
    assert result.match_status == "api_error"
    assert "HTTPError" in result.error


def test_enrich_pubchem_exposes_lookup_name_and_synonyms(monkeypatch, tmp_path):
    class FakePubChemClient:
        def __init__(self, *args, **kwargs):  # noqa: ANN002, ANN003
            pass

        def enrich_name(self, name, properties, max_synonyms=25):  # noqa: ANN001
            return pubchem.PubChemResult(
                query=name,
                pubchem_cid="2336",
                iupac_name="benzo[a]pyrene",
                synonyms=["Benzo[a]pyrene", "BaP"],
                match_status="matched",
                source_url="https://pubchem.ncbi.nlm.nih.gov/compound/2336",
            )

    monkeypatch.setattr("carcinogen_harmonizer.enrich.PubChemClient", FakePubChemClient)
    df = pd.DataFrame(
        [
            {
                "carcinogen_id": "BaP",
                "preferred_query_name": "Benzo[a]pyrene",
                "display_label": "BaP",
                "node_label": "BaP",
                "is_chemical_like": True,
            }
        ]
    )

    out, provenance = enrich_pubchem(df, PipelineConfig(), cache_dir=tmp_path)

    assert out.loc[0, "pubchem_lookup_name"] == "Benzo[a]pyrene"
    assert out.loc[0, "pubchem_synonyms"] == "Benzo[a]pyrene; BaP"
    assert out.loc[0, "pubchem_cid"] == "2336"
    assert provenance[0]["query"] == "Benzo[a]pyrene"


def test_refresh_alias_data_updates_pubchem_provenance():
    class FakePubChemClient:
        def enrich_name(self, name, properties, max_synonyms=25):  # noqa: ANN001
            assert name == "Benzo[a]pyrene"
            return SimpleNamespace(
                match_status="matched",
                pubchem_cid="2336",
                source_url="https://pubchem.ncbi.nlm.nih.gov/compound/2336",
                error="",
            )

    data = {
        "metadata": {"retrieved_on": "2020-01-01"},
        "aliases": {
            "BaP": {
                "preferred_query_name": "Benzo[a]pyrene",
                "source": {
                    "database": "PubChem",
                    "record_type": "Compound",
                    "cid": "old",
                    "url": "https://example.test/old",
                },
            }
        },
    }

    refreshed, errors = refresh_alias_data(
        data,
        FakePubChemClient(),
        retrieved_on="2026-04-28",
    )

    assert errors == []
    assert refreshed["metadata"]["retrieved_on"] == "2026-04-28"
    assert refreshed["aliases"]["BaP"]["source"]["cid"] == "2336"
    assert refreshed["aliases"]["BaP"]["source"]["url"] == "https://pubchem.ncbi.nlm.nih.gov/compound/2336"


def test_refresh_alias_data_reports_pubchem_failures():
    class FakePubChemClient:
        def enrich_name(self, name, properties, max_synonyms=25):  # noqa: ANN001
            return SimpleNamespace(
                match_status="no_pubchem_match",
                pubchem_cid="",
                source_url="",
                error="",
            )

    refreshed, errors = refresh_alias_data(
        {"aliases": {"BAD": {"preferred_query_name": "missing"}}},
        FakePubChemClient(),
        retrieved_on="2026-04-28",
    )

    assert refreshed["aliases"]["BAD"]["preferred_query_name"] == "missing"
    assert errors == ["BAD: PubChem lookup for 'missing' failed (no_pubchem_match)"]
