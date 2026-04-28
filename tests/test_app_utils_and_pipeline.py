from __future__ import annotations

import pandas as pd

from carcinogen_harmonizer.app_utils import filter_dataset, iarc_group_summary
from carcinogen_harmonizer.pipeline import PipelineRunOptions, run_pipeline


def test_filter_dataset_searches_aliases_and_filters_groups():
    df = pd.DataFrame(
        [
            {
                "carcinogen_id": "BaP",
                "display_label": "Benzo[a]pyrene",
                "pubchem_synonyms": "Benzo[a]pyrene; BaP",
                "iarc_group_normalized": "Group 1",
                "carcinogen_class_normalized": "PAH",
                "agent_entity_type": "Chemical",
                "needs_manual_review": "False",
            },
            {
                "carcinogen_id": "TCE",
                "display_label": "Trichloroethylene",
                "pubchem_synonyms": "",
                "iarc_group_normalized": "Group 2A",
                "carcinogen_class_normalized": "Chlorinated_Solvent",
                "agent_entity_type": "Chemical",
                "needs_manual_review": "True",
            },
        ]
    )

    out = filter_dataset(df, "bap", groups=["Group 1"], manual_review="No review flag")

    assert len(out) == 1
    assert out.loc[0, "carcinogen_id"] == "BaP"


def test_iarc_group_summary_uses_known_group_order():
    df = pd.DataFrame({"iarc_group_normalized": ["Group 2B", "Group 1", "Group 2B"]})

    summary = iarc_group_summary(df)

    assert summary["iarc_group_normalized"].tolist() == ["Group 1", "Group 2B"]
    assert summary["count"].tolist() == [1, 2]


def test_run_pipeline_normalize_only_writes_standard_outputs(tmp_path):
    seed = tmp_path / "seed.csv"
    seed.write_text(
        "\n".join(
            [
                "node_id,node_label,node_type,carcinogen_class,iarc,Status",
                "TCDD,\"2,3,7,8-tetrachlorodibenzo-p-dioxin\",Carcinogen,Dioxin,Group 1,1",
            ]
        ),
        encoding="utf-8",
    )

    result = run_pipeline(
        PipelineRunOptions(
            seed=seed,
            out=tmp_path / "out",
            cache=tmp_path / "cache",
            normalize_only=True,
        )
    )

    assert len(result.df) == 1
    assert result.report["row_count"] == 1
    assert (tmp_path / "out" / "carcinogens_enriched.csv").exists()
