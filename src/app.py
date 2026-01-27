import argparse
import sys
from pathlib import Path
from typing import Optional, Literal
import json

import numpy as np
import pandas as pd
from pyensembl import species

from flex_probe_pipeline import HumanBackgroundFlexProbeConfig, MouseBackgroundFlexProbeConfig, FlexProbeDesigner, \
    SnvProbeHelper, FlexProbeConfig, MskImpactSnvProbeHelper, ManeSelectSnvProbeHelper

# Set random seed for reproducibility
np.random.seed(42)


def main(
        config_file: Path,
        organism: Literal["human", "mouse"],
        technology: Literal["FlexV1", "FlexV2", "VisiumHD"],
        barcodes: int,
        output_format: Literal["csv", "tsv", "xlsx"],
        targets: Path,
        name: str,
        skip_errors: bool = False,
        msk: bool = False,
        mane: bool = False,
        ensembl_release: int = 111
):
    if not targets.exists() or not targets.is_file():
        print(f"Error: Targets file '{targets}' does not exist or is not a file.", file=sys.stderr)
        exit(1)

    # Load configuration
    if not config_file.exists():
        config_data = dict()
    else:
        with open(config_file, 'r') as f:
            config_data = json.load(f)

    if organism == "human":
        config, release = HumanBackgroundFlexProbeConfig(ensembl_release=ensembl_release, **config_data)
    elif organism == "mouse":
        config, release = MouseBackgroundFlexProbeConfig(ensembl_release=ensembl_release, **config_data)
    else:
        config = FlexProbeConfig(**config_data)
        release = None

    match technology:
        case "FlexV1":
            if organism == "human":
                probeset = "human_flex_v1"
            elif organism == "mouse":
                probeset = "mouse_flex_v1"
            else:
                raise ValueError(f"Unsupported organism for FlexV1: {organism}")
        case "FlexV2":
            if organism == "human":
                probeset = "human_flex_v2"
            elif organism == "mouse":
                probeset = "mouse_flex_v2"
            else:
                raise ValueError(f"Unsupported organism for FlexV2: {organism}")
        case "VisiumHD":
            if organism == "human":
                probeset = "human_visiumhd"
            elif organism == "mouse":
                probeset = "mouse_visiumhd"
            else:
                raise ValueError(f"Unsupported organism for VisiumHD: {organism}")

    # noinspection PyTypeChecker
    designer = FlexProbeDesigner(
        working_dir="./",
        reference_probe_set=probeset,
        config=config
    )

    if organism == "human":
        if msk:
            snv_probe_helper = MskImpactSnvProbeHelper(config, release)
        elif mane:
            snv_probe_helper = ManeSelectSnvProbeHelper(config, release)
        else:
            snv_probe_helper = SnvProbeHelper(config, release)
    else:
        snv_probe_helper = SnvProbeHelper(config, release)

    targets = pd.read_csv(targets, header=None)
    if targets.shape[0] == 0:
        print(f"Error: Targets file '{targets}' is empty.", file=sys.stderr)
        exit(1)
    if targets.shape[1] == 1:
        targets.columns = ['HGVSc']
        targets['Sequence'] = None
        targets['Gene'] = targets.HGVSc.str.split(" ").str[0]
    elif targets.shape[1] == 2:
        targets.columns = ['HGVSc', 'Sequence']
        targets['Gene'] = targets.HGVSc.str.split(" ").str[0]
    elif targets.shape[1] == 3:
        targets.columns = ['Gene', 'HGVSc', 'Sequence']
    else:
        print(f"Error: Targets file '{targets}' has too many columns (max 3).", file=sys.stderr)
        exit(1)

    # Run the main pipeline
    transcripts = dict()
    target_starts = []
    target_ends = []
    skipped = []
    for index, row in targets.iterrows():
        gene = row['Gene']
        hgvsc = row['HGVSc']
        sequence: Optional[str] = row['Sequence'] if 'Sequence' in row and pd.notna(row['Sequence']) else None

        is_zerobp = False
        if len(hgvsc.split(" ")) == 1 or '0bp' in hgvsc:
            is_zerobp = True

        if is_zerobp:
            hgvsc = f"{gene} 0bp"

        print(hgvsc)

        try:
            snv_start, snv_end, snv_action, snv_data, mutated_snv_length = snv_probe_helper.parse_snv_info(hgvsc.split(" ")[1])
            original_sequence, mutated_sequence, mutated_snv_start, mutated_snv_end = snv_probe_helper.get_gene_sequence(
                gene, snv_start, snv_end, snv_action, snv_data, sequence
            )
        except Exception as e:
            if skip_errors:
                print(f"Warning: Skipping '{hgvsc}': {e}", file=sys.stderr)
                skipped.append((hgvsc, str(e)))
                continue
            else:
                raise

        transcripts[hgvsc] = (original_sequence, mutated_sequence)
        target_starts.append((snv_start, mutated_snv_start))
        target_ends.append((snv_end, mutated_snv_end))

    if skipped:
        print(f"\nSkipped {len(skipped)} targets due to errors:", file=sys.stderr)
        for target, reason in skipped:
            print(f"  {target}: {reason}", file=sys.stderr)
        print(file=sys.stderr)

    probe_df = designer.generate_gapfilling_flex_probe_set_df(
        transcripts,
        target_starts,
        target_ends,
        expect_hits=[],
        n_probes=1,
        visium='visium' in technology,
        barcode=barcodes
    )

    if output_format == "csv":
        probe_df.to_csv(f"{name}.csv", index=False)
    elif output_format == "tsv":
        probe_df.to_csv(f"{name}.tsv", index=False, sep="\t")
    elif output_format == "xlsx":
        probe_df.to_excel(f"{name}.xlsx", index=False)

    print(f"Generated {len(probe_df)} probes and saved to {name}.{output_format}'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the Flex Probe Design pipeline with specified configuration."
    )
    parser.add_argument(
        "--config_file",
        type=Path,
        help="Path to the Flex Probe Design configuration file (JSON format).",
    )
    parser.add_argument(
        "--organism",
        type=str,
        help="Organism name for the analysis (e.g., 'human', 'mouse').",
        choices=["human", "mouse"],
        default="human",
    )
    parser.add_argument(  # Flex v1, Flex v2, or visiumhd
        "--technology",
        type=str,
        help="Technology used for the analysis (e.g., 'FlexV1', 'FlexV2', 'VisiumHD').",
        default="FlexV1",
        choices=["FlexV1", "FlexV2", "VisiumHD"],
    )
    parser.add_argument(
        "--barcodes",
        type=int,
        help="Number of barcodes used in the experiment, currently only supports Flex v1 multiplexing.",
        default=1
    )

    parser.add_argument(
        "--output_format",
        type=str,
        help="Output format for the designed probes",
        default="tsv",
        choices=["csv", "tsv", "xlsx"]
    )

    parser.add_argument(
        "--skip_errors",
        action="store_true",
        help="Skip targets that fail (e.g., gene symbol not found, invalid variant) instead of aborting."
    )

    parser.add_argument(
        "--msk",
        action="store_true",
        help="Use MSKCC's internal canonical transcript override when generating probes."
    )

    parser.add_argument(
        "--mane",
        action="store_true",
        help="Use MANE Select canonical transcripts (recommended for reproducibility)."
    )

    parser.add_argument(
        "--release",
        type=int,
        default=111,
        help="Ensembl release version to use (default: 111)."
    )

    parser.add_argument(
        "targets",
        type=Path,
        help="Path to a csv file containing a list of genes, one per line (no header)."
             " If a single column, it is assumed to be HGVSc-like notation, you may use either gene symbols or Ensembl transcript IDs. "
             " If two columns, the first is HGVSc-like notation and the second is full transcript sequence."
             " If three columns, first is gene name, second is HGVSc-like notation, third is full transcript sequence.",
    )

    parser.add_argument(
        "name",
        type=str,
        help="Base name for the output files (without extension).",
    )

    args = parser.parse_args()

    main(args.config_file, args.organism, args.technology, args.barcodes, args.output_format, args.targets, args.name, args.skip_errors, args.msk, args.mane, args.release)
