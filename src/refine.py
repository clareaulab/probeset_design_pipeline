import argparse
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Literal, Dict, Any, List, Tuple
from scipy.optimize import minimize
import time

# Import exactly as requested/shown in app.py
from flex_probe_pipeline import (
    HumanBackgroundFlexProbeConfig,
    MouseBackgroundFlexProbeConfig,
    FlexProbeDesigner,
    SnvProbeHelper,
    FlexProbeConfig,
    MskImpactSnvProbeHelper,
    ManeSelectSnvProbeHelper
)

# Set random seed
np.random.seed(42)


def cleanup_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes the 'name' column in the ground truth dataframe by removing
    common suffixes like _0, _1, _LHS, _RHS, etc.
    """
    if 'name' not in df.columns:
        return df

    df = df.copy()
    # Ensure column is string type before regex operations
    df['name'] = df['name'].astype(str)

    df['name'] = (
        df['name'].str.replace(r"(?<=.)_0 ", " ", regex=True)
        .str.replace(r"(?<=.)_1 ", " ", regex=True)
        .str.replace(r"_0$", " ", regex=True)
        .str.replace(r"_1$", " ", regex=True)
        .str.replace("_LHS", "", regex=False)
        .str.replace("_RHS", "", regex=False)
        .str.replace("_both", "", regex=False)
        .str.strip()
    )
    return df


def levenshtein_distance(s1: str, s2: str) -> int:
    """
    Calculates the Levenshtein distance between two strings.
    """
    if s1 == s2:
        return 0
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def get_designer_and_config(
        config_dict: Dict[str, Any],
        organism: str,
        technology: str,
        ensembl_release: int,
        probeset: str
) -> FlexProbeDesigner:
    """
    Instantiates the Designer with a dynamic configuration dictionary.
    """
    if organism == "human":
        # We discard the genome object here as we only need the config for the designer
        config, _ = HumanBackgroundFlexProbeConfig(ensembl_release=ensembl_release, **config_dict)
    elif organism == "mouse":
        config, _ = MouseBackgroundFlexProbeConfig(ensembl_release=ensembl_release, **config_dict)
    else:
        config = FlexProbeConfig(**config_dict)

    # Initialize designer with the new config
    designer = FlexProbeDesigner(
        working_dir="./",
        reference_probe_set=probeset,
        config=config
    )
    return designer


def calculate_loss(
        params: np.ndarray,
        param_keys: List[str],
        base_config: Dict[str, Any],
        ground_truth_df: pd.DataFrame,
        # Pipeline specific args passed via partial or closure
        transcripts: Dict,
        target_starts: List,
        target_ends: List,
        organism: str,
        technology: str,
        ensembl_release: int,
        probeset: str,
        visium: bool,
        barcodes: int
) -> float:
    """
    Objective function: Updates config, runs designer, calculates Levenshtein distance.
    """
    # Initialize static counters if they don't exist
    if not hasattr(calculate_loss, "iters"):
        calculate_loss.iters = 0
        calculate_loss.best_loss = float('inf')
        calculate_loss.start_time = time.time()

    calculate_loss.iters += 1
    iter_start = time.time()

    # 1. Update Configuration
    current_config = base_config.copy()
    param_print_str = []
    for key, val in zip(param_keys, params):
        # Constraints: penalties usually shouldn't be negative
        val = abs(float(val))
        current_config[key] = val
        param_print_str.append(f"{key}={val:.2f}")

    print(f"\n--- Iteration {calculate_loss.iters} ---")
    print(f"Testing Params: {', '.join(param_print_str)}")

    # 2. Run Designer (Probe Selection Step only)
    try:
        designer = get_designer_and_config(current_config, organism, technology, ensembl_release, probeset)

        # Note: app.py passes expect_hits=[] by default
        probe_df = designer.generate_gapfilling_flex_probe_set_df(
            transcripts,
            target_starts,
            target_ends,
            expect_hits=[],
            n_probes=1,
            visium=visium,
            barcode=barcodes
        )
    except Exception as e:
        print(f"!! Iteration {calculate_loss.iters} Failed (Pipeline Crash): {e}")
        return 1e6

    if probe_df.empty:
        print(f"!! Iteration {calculate_loss.iters} Failed (No probes generated)")
        return 1e6

    # 3. Compare with Ground Truth
    gen_lhs_col = 'lhs_probe' if 'lhs_probe' in probe_df.columns else 'lhs_sequence'
    gen_rhs_col = 'rhs_probe' if 'rhs_probe' in probe_df.columns else 'rhs_sequence'

    true_lhs_col = 'lhs_probe' if 'lhs_probe' in ground_truth_df.columns else 'lhs_sequence'
    true_rhs_col = 'rhs_probe' if 'rhs_probe' in ground_truth_df.columns else 'rhs_sequence'

    try:
        merged = pd.merge(
            probe_df,
            ground_truth_df,
            on='name',
            how='inner',
            suffixes=('_gen', '_true')
        )
    except KeyError:
        print("!! Merge Error: 'name' column missing.")
        return 1e6

    if merged.empty:
        print(f"!! Iteration {calculate_loss.iters} Failed (No matching targets in Ground Truth)")
        return 1e6

    # Calculate Distance
    total_dist = 0
    for _, row in merged.iterrows():
        l_gen_key = f"{gen_lhs_col}_gen" if gen_lhs_col == true_lhs_col else gen_lhs_col
        l_true_key = f"{true_lhs_col}_true" if gen_lhs_col == true_lhs_col else true_lhs_col
        r_gen_key = f"{gen_rhs_col}_gen" if gen_rhs_col == true_rhs_col else gen_rhs_col
        r_true_key = f"{true_rhs_col}_true" if gen_rhs_col == true_rhs_col else true_rhs_col

        lhs_gen = str(row.get(l_gen_key, ''))
        lhs_true = str(row.get(l_true_key, ''))
        d_lhs = levenshtein_distance(lhs_gen, lhs_true)

        rhs_gen = str(row.get(r_gen_key, ''))
        rhs_true = str(row.get(r_true_key, ''))
        d_rhs = levenshtein_distance(rhs_gen, rhs_true)

        total_dist += (d_lhs + d_rhs)

    # Normalize by number of targets
    loss = total_dist / len(merged)

    elapsed = time.time() - iter_start
    is_best = ""
    if loss < calculate_loss.best_loss:
        calculate_loss.best_loss = loss
        is_best = " [NEW BEST]"

    print(f"Result: Loss={loss:.4f}{is_best} (Time: {elapsed:.1f}s)")

    return loss


def main():
    parser = argparse.ArgumentParser(description="Optimize FlexProbeConfig parameters.")
    parser.add_argument("config_file", type=Path, help="Path to initial config.json")
    parser.add_argument("inputs", type=Path, help="Path to targets CSV")
    parser.add_argument("ground_truth", type=Path,
                        help="Path to ground truth (CSV/TSV with name, lhs_probe/sequence, rhs_probe/sequence)")
    parser.add_argument("output_config", type=Path, help="Path to save optimized config")

    # Standard app.py arguments
    parser.add_argument("--organism", choices=["human", "mouse"], default="human")
    parser.add_argument("--technology", choices=["FlexV1", "FlexV2", "VisiumHD"], default="FlexV1")
    parser.add_argument("--barcodes", type=int, default=1)
    parser.add_argument("--release", type=int, default=111)
    parser.add_argument("--msk", action="store_true")
    parser.add_argument("--mane", action="store_true")

    args = parser.parse_args()

    # 1. Load Initial Configuration
    if not args.config_file.exists():
        print("Config file not found.")
        sys.exit(1)

    with open(args.config_file, 'r') as f:
        base_config = json.load(f)

    # 2. Determine Probeset String (Logic from app.py)
    if args.technology == "FlexV1":
        probeset = f"{args.organism}_flex_v1"
    elif args.technology == "FlexV2":
        probeset = f"{args.organism}_flex_v2"
    elif args.technology == "VisiumHD":
        probeset = f"{args.organism}_visiumhd"
    else:
        probeset = "unknown"

    # 3. Initialize Helper and Parse Targets (PRE-OPTIMIZATION STEP)
    if args.organism == "human":
        temp_config, genome_obj = HumanBackgroundFlexProbeConfig(ensembl_release=args.release, **base_config)
    elif args.organism == "mouse":
        temp_config, genome_obj = MouseBackgroundFlexProbeConfig(ensembl_release=args.release, **base_config)
    else:
        temp_config = FlexProbeConfig(**base_config)
        genome_obj = None

    # Init Helper
    if args.organism == "human":
        if args.msk:
            snv_probe_helper = MskImpactSnvProbeHelper(temp_config, genome_obj)
        elif args.mane:
            snv_probe_helper = ManeSelectSnvProbeHelper(temp_config, genome_obj)
        else:
            snv_probe_helper = SnvProbeHelper(temp_config, genome_obj)
    else:
        snv_probe_helper = SnvProbeHelper(temp_config, genome_obj)

    # Parse Targets CSV
    targets_df = pd.read_csv(args.inputs, header=None)
    if targets_df.shape[1] == 1:
        targets_df.columns = ['HGVSc']
        targets_df['Sequence'] = None
        targets_df['Gene'] = targets_df.HGVSc.str.split(" ").str[0]
    elif targets_df.shape[1] == 2:
        targets_df.columns = ['HGVSc', 'Sequence']
        targets_df['Gene'] = targets_df.HGVSc.str.split(" ").str[0]
    elif targets_df.shape[1] == 3:
        targets_df.columns = ['Gene', 'HGVSc', 'Sequence']

    # Pre-calculate sequences
    transcripts = dict()
    target_starts = []
    target_ends = []

    print("Pre-calculating transcript sequences...")
    for _, row in targets_df.iterrows():
        gene = row['Gene']
        hgvsc = row['HGVSc']
        sequence = row['Sequence'] if 'Sequence' in row and pd.notna(row['Sequence']) else None

        is_zerobp = False
        if len(hgvsc.split(" ")) == 1 or '0bp' in hgvsc:
            is_zerobp = True
            hgvsc = f"{gene} 0bp"

        try:
            snv_start, snv_end, snv_action, snv_data, _ = snv_probe_helper.parse_snv_info(hgvsc.split(" ")[1])
            original_sequence, mutated_sequence, mutated_snv_start, mutated_snv_end = snv_probe_helper.get_gene_sequence(
                gene, snv_start, snv_end, snv_action, snv_data, sequence
            )

            if gene.startswith("ENST"):
                transcript = snv_probe_helper.genome.transcript_by_id(gene.split(".")[0])
                gene_name = transcript.gene_name
                hgvsc = hgvsc.replace(gene, gene_name, 1)

            transcripts[hgvsc] = (original_sequence, mutated_sequence)
            target_starts.append((snv_start, mutated_snv_start))
            target_ends.append((snv_end, mutated_snv_end))

        except Exception as e:
            print(f"Skipping {hgvsc}: {e}")
            continue

    print(f"Pre-calculation complete. {len(transcripts)} valid targets.")

    # 4. Load Ground Truth (Robust loading for TSV/CSV)
    gt_path = args.ground_truth
    try:
        if gt_path.suffix.lower() == '.tsv':
            gt_df = pd.read_csv(gt_path, sep='\t')
        elif gt_path.suffix.lower() in ['.xls', '.xlsx']:
            gt_df = pd.read_excel(gt_path)
        else:
            # Try comma, if it looks like one column, try tab
            gt_df = pd.read_csv(gt_path)
            if len(gt_df.columns) == 1:
                # Rewind and try tab
                gt_df = pd.read_csv(gt_path, sep='\t')
    except Exception as e:
        print(f"Error loading ground truth file: {e}")
        sys.exit(1)

    # Normalize column names (strip whitespace, lowercase)
    gt_df.columns = [str(c).lower().strip() for c in gt_df.columns]

    if 'name' not in gt_df.columns:
        print(f"Error: 'name' column not found in ground truth. Found: {list(gt_df.columns)}")
        sys.exit(1)

    # --- APPLY CLEANUP ---
    print("Normalizing ground truth names...")
    gt_df = cleanup_names(gt_df)

    print(f"Loaded {len(gt_df)} ground truth records.")

    # 5. Optimization Setup
    keys_to_optimize = [
        "lenient_gc_penalty",
        "last_base_penalty",
        "first_base_penalty",
        "suboptimal_gap_penalty",
        "invalid_gap_penalty",
        "distance_from_max_bridge_length_penalty",
        "unbalanced_gc_penalty",
        "homopolymer_penalty",
        "lhs_tm_penalty",
        "rhs_tm_penalty",
        "flex_overlap_penalty",
        "tandem_repeat_penalty"
    ]

    keys_to_optimize = [k for k in keys_to_optimize if k in base_config]
    initial_values = [base_config[k] for k in keys_to_optimize]

    print(f"Optimizing {len(keys_to_optimize)} parameters.")

    # 6. Run Optimization
    print(f"Optimizing {len(keys_to_optimize)} parameters.")

    res = minimize(
        calculate_loss,
        x0=np.array(initial_values),
        args=(
            keys_to_optimize,
            base_config,
            gt_df,
            transcripts,
            target_starts,
            target_ends,
            args.organism,
            args.technology,
            args.release,
            probeset,
            ('Visium' in args.technology),
            args.barcodes
        ),
        method='Nelder-Mead',
        # CHANGE HERE: Use maxfev to limit total pipeline runs, not just algorithm steps
        options={'maxfev': 50, 'disp': True}
    )

    print("\nOptimization Complete.")
    print(f"Final Loss: {res.fun}")
    print(f"Best Params: {res.x}")

    # 7. Save Result
    final_config = base_config.copy()
    for key, val in zip(keys_to_optimize, res.x):
        final_config[key] = abs(float(val))

    with open(args.output_config, 'w') as f:
        json.dump(final_config, f, indent=2)

    print(f"Saved optimized config to {args.output_config}")


if __name__ == "__main__":
    main()