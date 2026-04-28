from __future__ import annotations

import argparse
from pathlib import Path

from . import (
    DEFAULT_IARC_URL,
    __version__,
)
from .pipeline import PipelineRunOptions, run_pipeline


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="carcinogen-harmonizer",
        description="Build a normalized/enriched carcinogen dataset.",
    )
    parser.add_argument(
        "--seed",
        required=True,
        help="Curated seed CSV that anchors stable IDs and provenance; must include node_id,node_label,node_type,carcinogen_class,iarc,Status",
    )
    parser.add_argument("--out", default="outputs", help="Output directory")
    parser.add_argument("--config", default="config.yaml", help="YAML config path")
    parser.add_argument("--cache", default="cache", help="Cache directory")
    parser.add_argument("--no-pubchem", action="store_true", help="Skip PubChem online enrichment")
    parser.add_argument(
        "--no-classyfire",
        action="store_true",
        help="Skip ClassyFire (ChemOnt) chemical-class enrichment",
    )
    parser.add_argument(
        "--classyfire-cache",
        default=None,
        help="Override ClassyFire cache directory (default: <cache>/classyfire)",
    )
    parser.add_argument("--iarc-csv", default=None, help="Optional local IARC reference CSV for matching existing rows")
    parser.add_argument("--ntp-csv", default=None, help="Optional local NTP Report on Carcinogens reference CSV")
    parser.add_argument("--oehha-csv", default=None, help="Optional local OEHHA Prop 65 reference CSV")
    parser.add_argument("--ctd-diseases", default=None, help="Optional CTD chemical-disease TSV/TSV.GZ")
    parser.add_argument("--ctd-genes", default=None, help="Optional CTD chemical-gene TSV/TSV.GZ")
    parser.add_argument("--expand-iarc", action="store_true", help="Append IARC agents that are not already in the seed CSV")
    parser.add_argument("--iarc-source", default=DEFAULT_IARC_URL, help="IARC source for --expand-iarc: local CSV/XLSX/TSV or the IARC classifications URL")
    parser.add_argument("--include-group3", action="store_true", help="Include IARC Group 3 records during --expand-iarc")
    parser.add_argument("--iarc-groups", default="Group 1,Group 2A,Group 2B", help="Comma-separated IARC groups to include for --expand-iarc")
    parser.add_argument(
        "--output-format",
        choices=["csv", "json", "both"],
        default="both",
        help="Primary table and graph output format(s)",
    )
    parser.add_argument("--version", action="version", version=f"carcinogen-harmonizer {__version__}")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    include_groups = [g.strip() for g in args.iarc_groups.split(",") if g.strip()]
    result = run_pipeline(
        PipelineRunOptions(
            seed=args.seed,
            out=args.out,
            config=args.config,
            cache=args.cache,
            use_pubchem=not args.no_pubchem,
            use_classyfire=not args.no_classyfire,
            classyfire_cache=Path(args.classyfire_cache) if args.classyfire_cache else None,
            iarc_csv=args.iarc_csv,
            ntp_csv=args.ntp_csv,
            oehha_csv=args.oehha_csv,
            ctd_diseases=args.ctd_diseases,
            ctd_genes=args.ctd_genes,
            expand_iarc=args.expand_iarc,
            iarc_source=args.iarc_source,
            include_group3=args.include_group3,
            iarc_groups=include_groups,
            output_format=args.output_format,
        )
    )
    if args.expand_iarc:
        print(
            f"Loaded {result.iarc_records_loaded} IARC records "
            f"and appended {result.appended_iarc_count} new rows"
        )
    print(f"Wrote {len(result.df)} rows to {result.output_dir / 'carcinogens_enriched.csv'}")
    print(f"Manual review count: {result.report['manual_review_count']}")


if __name__ == "__main__":
    main()
