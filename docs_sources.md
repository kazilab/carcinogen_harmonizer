# Source notes

## IARC Monographs

Primary source for hazard classification. The pipeline supports IARC as either:

1. a local CSV/XLSX/TSV passed with `--iarc-source` and `--expand-iarc`, or
2. a URL passed with `--iarc-source`.

The IARC `list-of-classifications` page may be rendered by browser-side JavaScript, so local CSV/XLSX export is the preferred reproducible source.

Suggested fields:

- `Agent`
- `Group`
- `Volume`
- `Year`
- `CAS` or `CASRN`

The public CLI is `carcinogen-harmonizer` and the public Python import namespace is `carcinogen_harmonizer`.

## Package layout

Runtime code lives under the single public import namespace `carcinogen_harmonizer`.
The repository-root `data/` directory is intentionally kept outside the package because it is example/user input, not library code. Runtime reference files that must ship with the package live under `carcinogen_harmonizer/data/`.

## PubChem

Used only for chemical-like agents. Non-chemical IARC agents such as exposures, occupations, biological agents, and radiation types receive `pubchem_match_status = skipped_non_chemical_agent`.
Name lookups that do not resolve now return `pubchem_match_status = no_pubchem_match`; network or service failures remain `api_error`.

The packaged alias map (`carcinogen_harmonizer/data/alias_map.json`) is also PubChem-backed. It maps local short IDs and abbreviations to PubChem-resolvable query names and records the PubChem Compound CID/URL for each entry. Local seed files are not treated as authority for alias-map entries.
Use `python tools/refresh_alias_map.py` to refresh those PubChem CIDs/URLs, or `python tools/refresh_alias_map.py --check` to validate that the JSON is current. PubChem enrichment output also includes `pubchem_lookup_name` and `pubchem_synonyms`, which expose the matched lookup string and PubChem-returned aliases after the bootstrap lookup succeeds.

## Output schema

Graph-ready exports now include `standard_node_id`, typo-compatible `srandard_node_id`, and `norm_carcinogen_class` aliases for downstream consumers that expect explicit standardized names.

## Streamlit app

`streamlit_app.py` is a browser UI over the same pipeline functions used by the CLI. It can search existing static output files or run partial/full pipeline jobs and write the standard output schema.

## Seed CSV

The seed CSV is the curated starting catalog. It anchors stable IDs, labels, and provenance before optional enrichment or expansion. The current pipeline expects it as the entry point rather than bootstrapping from an empty workspace.

## CTD

Optional CTD files can enrich chemical-disease and chemical-gene associations with `--ctd-diseases` and `--ctd-genes`.

## Local regulatory tables

NTP, OEHHA, ECHA, EPA IRIS, and other regulatory lists can be loaded as local normalized tables. Add source-specific modules if their schemas are stable enough for direct ingestion.
