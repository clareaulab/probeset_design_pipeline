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


def collapse_nearby_targets(
        transcripts: dict,
        transcript_ids: dict,
        target_positions: dict,
        collapse_window: int,
) -> tuple[dict, dict, list, list]:
    """
    Merge targets that share the same original transcript sequence and whose HGVSc SNV positions
    span <= collapse_window bp into a single probe entry covering the combined region.

    The combined probe is designed on the wildtype (original) sequence so its gap spans all
    variant sites, allowing a single probe pair to detect any variant in the group.
    0bp (wildtype) entries are never collapsed with SNV entries.

    target_positions: dict mapping hgvsc key -> (orig_start, orig_end, mut_start, mut_end)
    """
    from collections import defaultdict

    entries = []
    for key in transcripts.keys():
        orig_seq, mut_seq = transcripts[key]
        orig_start, orig_end, mut_start, mut_end = target_positions[key]
        is_zerobp = orig_start is None
        entries.append({
            'key': key,
            'orig_seq': orig_seq,
            'mut_seq': mut_seq,
            'orig_start': orig_start,
            'orig_end': orig_end,
            'mut_start': mut_start,
            'mut_end': mut_end,
            'tid': transcript_ids.get(key),
            'is_zerobp': is_zerobp,
        })

    # Group by original sequence; 0bp entries form their own singleton groups
    groups_by_seq = defaultdict(list)
    zerobp_entries = []
    for e in entries:
        if e['is_zerobp']:
            zerobp_entries.append(e)
        else:
            groups_by_seq[e['orig_seq']].append(e)

    new_transcripts = {}
    new_transcript_ids = {}
    new_target_starts = []
    new_target_ends = []

    def _emit(e):
        new_transcripts[e['key']] = (e['orig_seq'], e['mut_seq'])
        new_transcript_ids[e['key']] = e['tid']
        new_target_starts.append((e['orig_start'], e['mut_start']))
        new_target_ends.append((e['orig_end'], e['mut_end']))

    # Process SNV groups
    for orig_seq, group in groups_by_seq.items():
        group.sort(key=lambda e: e['orig_start'])

        # Greedy merge: extend the current cluster while its span stays within collapse_window
        clusters = []
        current = [group[0]]
        for e in group[1:]:
            cur_min = min(x['orig_start'] for x in current)
            cur_max = max(x['orig_end'] for x in current)
            new_min = min(cur_min, e['orig_start'])
            new_max = max(cur_max, e['orig_end'])
            if new_max - new_min <= collapse_window:
                current.append(e)
            else:
                clusters.append(current)
                current = [e]
        clusters.append(current)

        for cluster in clusters:
            if len(cluster) == 1:
                _emit(cluster[0])
            else:
                gene_prefix = cluster[0]['key'].split(' ')[0]
                variants_part = '+'.join(
                    e['key'].split(' ', 1)[1] if ' ' in e['key'] else e['key']
                    for e in cluster
                )
                combined_key = f"{gene_prefix} {variants_part}"
                combined_start = min(e['orig_start'] for e in cluster)
                combined_end = max(e['orig_end'] for e in cluster)
                print(f"Collapsing {len(cluster)} targets into '{combined_key}' "
                      f"(HGVSc positions {combined_start}-{combined_end})")
                new_transcripts[combined_key] = (orig_seq, orig_seq)
                new_transcript_ids[combined_key] = cluster[0]['tid']
                new_target_starts.append((combined_start, combined_start))
                new_target_ends.append((combined_end, combined_end))

    # Re-append 0bp entries in original order
    for e in zerobp_entries:
        _emit(e)

    return new_transcripts, new_transcript_ids, new_target_starts, new_target_ends


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
        ensembl_release: int = 111,
        search_method: Literal["brute_force", "optimization"] = None,
        fast: bool = False,
        collapse_window: Optional[int] = None,
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

    # Use cached flex_probes.csv next to config if available, otherwise download
    flex_probes_cache = config_file.parent / "flex_probes.csv" if config_file else None
    if flex_probes_cache and flex_probes_cache.exists():
        print(f"Loading reference probes from '{flex_probes_cache}'...")
        reference_probes = pd.read_csv(flex_probes_cache, comment='#')
    else:
        reference_probes = probeset

    # noinspection PyTypeChecker
    designer = FlexProbeDesigner(
        working_dir="./",
        reference_probe_set=reference_probes,
        config=config,
        fast=fast
    )

    # Cache the reference probes next to the config file if not already cached
    if flex_probes_cache and not flex_probes_cache.exists():
        flex_probes_cache.parent.mkdir(parents=True, exist_ok=True)
        designer.reference_probes.to_csv(flex_probes_cache, index=False)
        print(f"Saved reference probes to '{flex_probes_cache}'.")

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
    transcript_ids = dict()  # Map hgvsc -> transcript_id
    target_positions = dict()  # Map hgvsc -> (orig_start, orig_end, mut_start, mut_end)
    target_starts = []
    target_ends = []
    skipped = []
    first = True
    for index, row in targets.iterrows():
        gene = row['Gene']
        hgvsc = row['HGVSc']
        sequence: Optional[str] = row['Sequence'] if 'Sequence' in row and pd.notna(row['Sequence']) else None

        # Normalize colon-separated format: ENST00000378444.4:c.4376A>G -> ENST00000378444.4 c.4376A>G
        if ' ' not in hgvsc and ':' in hgvsc:
            transcript_part, hgvsc_part = hgvsc.split(':', 1)
            gene = transcript_part
            hgvsc = f"{transcript_part} {hgvsc_part}"

        is_zerobp = False
        if len(hgvsc.split(" ")) == 1 or '0bp' in hgvsc:
            is_zerobp = True

        if is_zerobp:
            hgvsc = f"{gene} 0bp"

        try:
            snv_start, snv_end, snv_action, snv_data, mutated_snv_length = snv_probe_helper.parse_snv_info(hgvsc.split(" ")[-1])
            original_sequence, mutated_sequence, mutated_snv_start, mutated_snv_end = snv_probe_helper.get_gene_sequence(
                gene, snv_start, snv_end, snv_action, snv_data, sequence
            )

            # Get the transcript ID that was used
            used_transcript_id = snv_probe_helper.last_transcript_id

            # If input was an ENST, replace with gene name in the output name for readability
            if gene.startswith("ENST"):
                transcript = snv_probe_helper.genome.transcript_by_id(gene.split(".")[0])
                gene_name = transcript.gene_name
                hgvsc = hgvsc.replace(gene, gene_name, 1)
        except Exception as e:
            if first:
                first = False  # If there is a header, silently skip
                continue
            if skip_errors:
                print(f"Warning: Skipping '{hgvsc}': {e}", file=sys.stderr)
                skipped.append((hgvsc, str(e)))
                continue
            else:
                raise

        print("Parsed", hgvsc)

        transcripts[hgvsc] = (original_sequence, mutated_sequence)
        transcript_ids[hgvsc] = used_transcript_id
        target_starts.append((snv_start, mutated_snv_start))
        target_ends.append((snv_end, mutated_snv_end))
        target_positions[hgvsc] = (snv_start, snv_end, mutated_snv_start, mutated_snv_end)
        first = False

    if skipped:
        print(f"\nSkipped {len(skipped)} targets due to errors:", file=sys.stderr)
        for target, reason in skipped:
            print(f"  {target}: {reason}", file=sys.stderr)
        print(file=sys.stderr)

    if collapse_window is not None:
        transcripts, transcript_ids, target_starts, target_ends = collapse_nearby_targets(
            transcripts, transcript_ids, target_positions, collapse_window
        )

    probe_df = designer.generate_gapfilling_flex_probe_set_df(
        transcripts,
        target_starts,
        target_ends,
        expect_hits=None,
        n_probes=1,
        visium='visium' in technology,
        barcode=barcodes,
        search_method=search_method
    )

    # Add transcript_id column by mapping from name
    # Match by exact name, name without suffix, or name prefix
    probe_df['transcript_id'] = probe_df['name'].apply(
        lambda n: transcript_ids.get(n) or transcript_ids.get(n.rsplit('_', 1)[0]) or
                  next((tid for hgvsc, tid in transcript_ids.items() if n.startswith(hgvsc)), None)
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
        "--blast",
        action="store_true",
        help="Enable BLAST-based off-target filtering (slower but more accurate; requires BLAST+ in PATH)."
    )

    search_method_group = parser.add_mutually_exclusive_group()
    search_method_group.add_argument(
        "--brute-force",
        action="store_true",
        dest="brute_force",
        help="Force brute force (exhaustive) search for probe optimization."
    )
    search_method_group.add_argument(
        "--optimize",
        action="store_true",
        help="Force dual annealing optimization search for probe optimization."
    )

    parser.add_argument(
        "--collapse_window",
        type=int,
        default=None,
        metavar="N",
        help="Merge targets on the same transcript whose SNV positions span <= N bp into a single "
             "probe. The combined probe is designed on the wildtype sequence with its gap spanning "
             "all variant sites, so one probe pair detects any variant in the group. "
             "Disabled by default. A value around your max_bridge_length (default 10) is a "
             "reasonable starting point."
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

    # Determine search method from mutually exclusive flags
    if args.brute_force:
        search_method = "brute_force"
    elif args.optimize:
        search_method = "optimization"
    else:
        search_method = None

    main(args.config_file, args.organism, args.technology, args.barcodes, args.output_format, args.targets, args.name, args.skip_errors, args.msk, args.mane, args.release, search_method, fast=not args.blast, collapse_window=args.collapse_window)
