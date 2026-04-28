from __future__ import annotations

from importlib import resources
import json

import pandas as pd

from carcinogen_harmonizer import build_dataset, graph_tables, load_config, make_qa_report, normalize_seed
from carcinogen_harmonizer.normalize import ALIAS_MAP


def test_normalize_seed_cleans_labels_and_adds_id_aliases():
    seed = pd.DataFrame(
        [
            {
                "node_id": "CUSTOM_1",
                "node_label": "1-<i>tert</i>-Butoxypropan-2-ol",
                "node_type": "Carcinogen",
                "carcinogen_class": "",
                "iarc": "Group 1",
                "Status": "1",
            },
            {
                "node_id": "CUSTOM_2",
                "node_label": "1-Naphthylthiourea (ANTU)",
                "node_type": "Carcinogen",
                "carcinogen_class": "",
                "iarc": "Group 2B",
                "Status": "1",
            },
        ]
    )

    out = normalize_seed(seed)
    assert out.loc[0, "display_label"] == "1-tert-Butoxypropan-2-ol"
    assert out.loc[0, "preferred_query_name"] == "1-tert-Butoxypropan-2-ol"
    assert out.loc[1, "preferred_query_name"] == "1-Naphthylthiourea"
    assert out.loc[0, "standard_node_id"] == "CUSTOM_1"
    assert out.loc[0, "srandard_node_id"] == "CUSTOM_1"


def test_alias_map_is_loaded_from_pubchem_sourced_json():
    text = (
        resources.files("carcinogen_harmonizer")
        .joinpath("data", "alias_map.json")
        .read_text(encoding="utf-8")
    )
    data = json.loads(text)

    assert data["metadata"]["primary_source"] == "PubChem"
    assert data["metadata"]["source_policy"].find("local seed files are not used") >= 0
    assert ALIAS_MAP["NNK"] == "4-(methylnitrosamino)-1-(3-pyridyl)-1-butanone"
    for entry in data["aliases"].values():
        assert entry["source"]["database"] == "PubChem"
        assert entry["source"]["cid"]
        assert entry["source"]["url"].startswith("https://pubchem.ncbi.nlm.nih.gov/compound/")


def test_build_dataset_exposes_normalized_aliases_and_graph_schema(tmp_path):
    seed = pd.DataFrame(
        [
            {
                "node_id": "TCDD",
                "node_label": "2,3,7,8-tetrachlorodibenzo-p-dioxin",
                "node_type": "Carcinogen",
                "carcinogen_class": "Dioxin",
                "iarc": "Group 1",
                "Status": "1",
            }
        ]
    )
    cfg = load_config("config.yaml")

    df, _ = build_dataset(
        seed,
        cfg=cfg,
        cache_dir=tmp_path / "cache",
        use_pubchem=False,
        use_classyfire=False,
    )

    assert "standard_node_id" in df.columns
    assert "srandard_node_id" in df.columns
    assert "norm_carcinogen_class" in df.columns
    assert df.loc[0, "norm_carcinogen_class"] == df.loc[0, "carcinogen_class_normalized"]

    graph = graph_tables(df)
    nodes = graph["nodes_carcinogens"]
    edges = graph["edges_carcinogen_facts"]

    assert {"standard_node_id", "srandard_node_id", "node_id", "node_label", "node_type"}.issubset(nodes.columns)
    assert {"source_standard_node_id", "source_srandard_node_id", "target_standard_node_id", "target_srandard_node_id"}.issubset(edges.columns)

    report = make_qa_report(df)
    assert "identifier_coverage" in report
    assert report["identifier_coverage"]["standard_node_id"]["coverage"] == 1.0
    assert report["standard_node_id_duplicate_count"] == 0
