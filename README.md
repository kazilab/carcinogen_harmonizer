# Carcinogen Harmonizer

<!-- PyPI version badge -->
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://carcinogen-harmonizer.streamlit.app)
<!--[![PyPI version](https://img.shields.io/pypi/v/ExposoGraph.svg)](https://pypi.org/project/ExposoGraph/)-->
<!--[![Documentation Status](https://readthedocs.org/projects/ExposoGraph/badge/?version=latest)](https://ExposoGraph.readthedocs.io/en/latest/?badge=latest)-->
<!--[![ResearchSquare](https://img.shields.io/badge/ResearchSquare-rs--9202489%2Fv1-00A0E0.svg)](https://www.researchsquare.com/article/rs-9202489/v1)-->
<!--[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.03.22.713456-b31b1b.svg)](https://doi.org/10.64898/2026.03.22.713456)-->
<!-- PyPI version badge -->
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-kazilab%2Fcarcinogen_harmonizer-181717?logo=github&logoColor=white)](https://github.com/kazilab/carcinogen_harmonizer)
[![@KaziLab.se](https://img.shields.io/website?url=https://www.kazilab.se/)](https://www.kazilab.se/)
<!-- PyPI version badge -->


Build a normalized carcinogen dataset from a curated seed CSV, optionally expand from IARC classifications, enrich with PubChem and ClassyFire (ChemOnt), and export graph-ready outputs for downstream analysis.

## Features

- Seed CSV normalization and QA flags
- Optional IARC expansion (`Group 1`, `2A`, `2B`, optional `3`)
- Optional PubChem enrichment (CID, SMILES, InChI, InChIKey, CAS, properties)
- **Optional ClassyFire / ChemOnt enrichment** (structure-based, validated chemical class with citable ChemOnt IDs)
- Optional local-reference enrichment (IARC/NTP/OEHHA)
- Optional CTD disease and gene association summaries
- `carcinogen_class_normalized` with explicit `carcinogen_class_source` provenance (`seed_csv`, `rules`, `classyfire`, `iarc_entity_type`, or `unknown`)
- Output in CSV and/or JSON (plus provenance JSONL and QA report JSON)
- Graph-ready nodes/edges tables (with `HAS_CHEMONT_CLASS` / `HAS_CHEMONT_SUPERCLASS` edges)

## Requirements

- Python 3.10+
- Packages:
  - `pandas`
  - `requests`
  - `pyyaml`
  - `lxml` (recommended for `pandas.read_html`)
  - `openpyxl` (if reading `.xlsx` IARC/reference files)

Install example:

```bash
pip install carcinogen-harmonizer
```

Public Python API:

```python
from carcinogen_harmonizer import load_config, build_dataset, graph_tables
```

## Seed CSV Format

The seed file is the curated starting catalog. The pipeline uses it to anchor stable node IDs, preserve curated labels and classes, and attach provenance before optional IARC/PubChem/ClassyFire enrichment runs.
It is currently required even if you plan to expand from IARC.

Minimum required columns:

- `node_id`
- `node_label`
- `node_type`
- `carcinogen_class`
- `iarc`
- `Status`

Example row:

```csv
node_id,node_label,node_type,carcinogen_class,iarc,Status
TCDD,"2,3,7,8-tetrachlorodibenzo-p-dioxin",Carcinogen,Dioxin,Group 1,1
```

> Note: Keep CSV values quoted as needed if labels contain commas.

## Alias Map Source

Short carcinogen IDs and common abbreviations are normalized through the packaged JSON file at `carcinogen_harmonizer/data/alias_map.json`.
This file is not generated from the local seed CSV. Each entry must carry an external source record; the current map is backed by PubChem Compound CIDs and PubChem URLs.
The alias map is only bootstrap data for early matching and PubChem lookup. After a PubChem match, the pipeline also exposes `pubchem_lookup_name` and `pubchem_synonyms` so downstream consumers can use the current PubChem names and aliases without overwriting curated seed labels.

Refresh or validate the packaged map with:

```bash
python tools/refresh_alias_map.py
python tools/refresh_alias_map.py --check
```

## Quick Start

Run with seed file only:

```bash
carcinogen-harmonizer --seed data/seed.csv --out outputs
```

Run with IARC expansion and both CSV+JSON outputs:

```bash
carcinogen-harmonizer \
  --seed data/seed.csv \
  --expand-iarc \
  --iarc-source "https://monographs.iarc.who.int/list-of-classifications" \
  --output-format both \
  --out outputs
```

Run with local reference sources:

```bash
carcinogen-harmonizer \
  --seed data/seed.csv \
  --expand-iarc \
  --iarc-source data/iarc_export.xlsx \
  --ntp-csv data/ntp_roc.csv \
  --oehha-csv data/oehha_prop65.csv \
  --ctd-diseases data/CTD_chemicals_diseases.tsv.gz \
  --ctd-genes data/CTD_chem_gene_ixns.tsv.gz \
  --output-format both \
  --out outputs
```

## Streamlit App

Run the browser app from the repository root:

```bash
streamlit run streamlit_app.py
```

The app has two workflows:

- Static Search: loads `outputs/carcinogens_enriched.csv`, shows IARC group counts, and filters by text, IARC group, carcinogen class, entity type, and review status.
- Run Pipeline: runs the same pipeline as the CLI with selectable scopes: normalize seed only, build without online enrichment, build with PubChem, or full PubChem + ClassyFire enrichment.

If Streamlit is not installed, install the optional app dependency:

```bash
pip install ".[app]"
```

## CLI Options

- `--seed` (required): seed CSV path
- `--out`: output directory (default: `outputs`)
- `--config`: YAML config path (default: `config.yaml`)
- `--cache`: cache directory for PubChem and ClassyFire responses (default: `cache`)
- `--no-pubchem`: skip PubChem API enrichment
- `--no-classyfire`: skip ClassyFire / ChemOnt enrichment
- `--classyfire-cache`: override ClassyFire cache directory
- `--iarc-csv`: optional local IARC table for matching existing rows
- `--ntp-csv`: optional local NTP Report on Carcinogens table
- `--oehha-csv`: optional local OEHHA Prop 65 table
- `--ctd-diseases`: optional CTD chemical-disease file (`.tsv` / `.tsv.gz`)
- `--ctd-genes`: optional CTD chemical-gene file (`.tsv` / `.tsv.gz`)
- `--expand-iarc`: append IARC agents not already in seed
- `--iarc-source`: IARC URL or local CSV/XLSX/TSV source
- `--include-group3`: include IARC Group 3 when expanding
- `--iarc-groups`: comma-separated groups to include
- `--output-format`: `csv`, `json`, or `both` (default: `both`)

## Outputs

Always generated:

- `carcinogens_enriched.csv`
- `manual_review_queue.csv` (if rows need review)
- `iarc_expansion_skipped_duplicates.csv` (when `--expand-iarc` used)
- `source_provenance.jsonl`
- `qa_report.json` (includes identifier coverage and duplicate-ID summaries)
- `nodes_carcinogens.csv` (standard `standard_node_id` / `node_id` / `node_label` / `node_type` columns plus compatibility aliases)
- `edges_carcinogen_facts.csv`

Generated when `--output-format json` or `both`:

- `carcinogens_enriched.json`
- `nodes_carcinogens.json`
- `edges_carcinogen_facts.json`

## Configuration

Edit `config.yaml` to control:

- PubChem behavior (`enabled`, timeout, rate limit, properties)
- ClassyFire behavior (`enabled`, timeout, rate limit, `upgrade_unknown_only`)
- QA checks (required identifiers, controlled IARC groups)

Example `config.yaml`:

```yaml
pubchem:
  enabled: true
  rate_limit_seconds: 0.2
classyfire:
  enabled: true
  rate_limit_seconds: 1.0
  upgrade_unknown_only: true   # set to false to overwrite regex-derived classes
```

## Chemical Class Provenance

The pipeline assigns `carcinogen_class_normalized` using a layered strategy and records the winning layer in `carcinogen_class_source`:

| Source              | Meaning                                                                   |
|---------------------|---------------------------------------------------------------------------|
| `seed_csv`          | Curated value from your `data/seed.csv`                                   |
| `rules`             | Project regex rules in `carcinogen_harmonizer.normalize.CLASS_RULES`      |
| `classyfire`        | Mapped from ClassyFire / ChemOnt taxonomy (validated by structure)        |
| `iarc_entity_type`  | Non-chemical IARC agent (Mixture / Exposure / Radiation / BiologicalAgent)|
| `unknown`           | No source matched; row likely needs manual review                         |

Each row also carries the raw ChemOnt taxonomy in `chemont_kingdom`, `chemont_superclass`, `chemont_class`, `chemont_subclass`, `chemont_direct_parent` (plus the matching `*_id` ChemOnt identifiers) and a citable `classyfire_reference_url`. A boolean `classyfire_validated` is `True` when ChemOnt independently agrees with the assigned bucket.

Compatibility aliases:

- `norm_carcinogen_class` mirrors `carcinogen_class_normalized`
- `standard_node_id` mirrors the stable internal `carcinogen_id`
- `srandard_node_id` is retained as a typo-compatible alias of `standard_node_id`

## Trust Model

- `seed_csv` rows are curated inputs and should be treated as the primary source of truth for identifiers and manually assigned classes.
- `rules` rows are heuristic and should be reviewed if you need citable evidence.
- `classyfire` rows are structure-backed but still depend on remote taxonomy availability.
- `no_pubchem_match` means a name lookup did not resolve; it is not the same as an API failure.
- `api_error` means the network or remote service failed.
- `confidence_score` is a deterministic heuristic for triage, not a statistical probability.

## Notes on IARC Source

- The IARC classifications page may be JavaScript-rendered.
- The pipeline attempts:
  1. direct HTML table parsing,
  2. auto-discovery of downloadable CSV/XLSX/TSV links.
- If neither works, export IARC classifications locally and pass with `--iarc-source`.

## Typical Workflow for Graph/Flux

1. Run pipeline with `--output-format both`.
2. Load `nodes_carcinogens.*` and `edges_carcinogen_facts.*` into your graph layer.
3. Use `qa_report.json` and `manual_review_queue.csv` to filter or triage uncertain records.
4. Use `source_provenance.jsonl` for traceability and auditing.
