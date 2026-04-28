from __future__ import annotations


def test_public_package_exports():
    import carcinogen_harmonizer as ch

    assert ch.__version__ == "0.4.0"
    assert callable(ch.build_dataset)
    assert callable(ch.load_config)
    assert callable(ch.graph_tables)
    assert callable(ch.run_pipeline)
    assert ch.DEFAULT_IARC_URL.startswith("https://")
