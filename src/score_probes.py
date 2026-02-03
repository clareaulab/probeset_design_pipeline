import argparse
import json
import sys
from pathlib import Path
from Bio.Seq import reverse_complement

# Ensure we can import from the local directory
sys.path.append(str(Path(__file__).parent.parent))

try:
    from flex_probe_pipeline import FlexProbeConfig
    from utils import max_homopolymer_length, get_melting_temp, has_tandem_repeat
except ImportError:
    # Fallback if running from a different directory structure
    sys.path.append(".")
    from flex_probe_pipeline import FlexProbeConfig
    from utils import max_homopolymer_length, get_melting_temp, has_tandem_repeat


def score_probe_with_breakdown(
        config: FlexProbeConfig,
        lhs_gene_sequence: str,
        rhs_gene_sequence: str,
        target_gap_gene_sequence: str = None
) -> dict:
    """
    Replicates the scoring logic from FlexProbeDesigner.score_probe but returns a dictionary breakdown.
    Expects inputs in GENE SENSE (Transcript orientation).
    """
    penalties = {}
    score = 0.0

    # If target_gap_gene_sequence is None, it is treated as a standard adjacent probe (Non-Gapfill)
    is_gapfill = target_gap_gene_sequence is not None

    # 1. Convert to Probe Space
    lhs_probe_sequence = reverse_complement(lhs_gene_sequence)
    rhs_probe_sequence = reverse_complement(rhs_gene_sequence)

    target_gap_probe_sequence = ""
    if is_gapfill:
        target_gap_probe_sequence = reverse_complement(target_gap_gene_sequence)

    # 2. Probe Length Penalties
    len_penalty = 0.0
    len_penalty += config.probe_length_penalty * abs(len(lhs_probe_sequence) - config.lhs_probe_length)
    len_penalty += config.probe_length_penalty * abs(len(rhs_probe_sequence) - config.rhs_probe_length)
    len_penalty += config.probe_length_penalty * abs(len(lhs_probe_sequence) - len(rhs_probe_sequence))
    penalties['probe_length_penalty'] = len_penalty
    score += len_penalty

    # 3. Identity Check
    if lhs_probe_sequence == rhs_probe_sequence:
        penalties['identity_error'] = config.invalid_score
        score += config.invalid_score

    # 4. Junction Checks
    junction_penalty = 0.0

    if is_gapfill:
        # --- Gap Fill Logic ---

        # Check gap length constraints
        if len(target_gap_probe_sequence) < config.min_bridge_length or len(
                target_gap_probe_sequence) > config.max_bridge_length:
            penalties['invalid_bridge_length'] = config.invalid_score
            score += config.invalid_score

        # Junction: Last base of gap + First base of RHS probe
        # (LHS -> Gap -> RHS)
        if len(target_gap_probe_sequence) > 0 and len(rhs_probe_sequence) > 0:
            junction = target_gap_probe_sequence[-1].upper() + rhs_probe_sequence[0].upper()
        else:
            junction = "NN"
    else:
        # --- Adjacent/Non-Gapfill Logic ---

        # Junction: Last base of LHS probe + First base of RHS probe
        # (LHS -> RHS)
        if len(lhs_probe_sequence) > 0 and len(rhs_probe_sequence) > 0:
            junction = lhs_probe_sequence[-1].upper() + rhs_probe_sequence[0].upper()
        else:
            junction = "NN"

    if junction not in ["TA", "TC", "TG", "TT", "CT", "CA", "AT"]:
        if junction in ["GA", "GT"]:
            junction_penalty += config.suboptimal_gap_penalty
        elif junction in ["AA", "AC", "AG"]:
            junction_penalty += config.suboptimal_gap_penalty
        elif junction in ["CC", "GC", "CG", "GG"]:
            junction_penalty += config.invalid_gap_penalty
        elif is_gapfill and len(target_gap_probe_sequence) > 1:
            junction_penalty += config.invalid_gap_penalty
        else:
            junction_penalty += config.invalid_score

    penalties['junction_penalty'] = junction_penalty
    penalties['junction_bases'] = junction
    score += junction_penalty

    # 5. Low Complexity (Gene Space)
    complexity_penalty = 0.0
    complexity_penalty += config.low_complexity_penalty * sum([c.islower() for c in lhs_gene_sequence])
    complexity_penalty += config.low_complexity_penalty * sum([c.islower() for c in rhs_gene_sequence])
    penalties['low_complexity_penalty'] = complexity_penalty
    score += complexity_penalty

    # 6. Gap Length Bias
    gap_bias_penalty = 0.0
    if is_gapfill:
        if config.bias_longer_gaps:
            gap_bias_penalty += config.distance_from_max_bridge_length_penalty * abs(
                len(target_gap_gene_sequence) - config.max_bridge_length)
        elif len(target_gap_gene_sequence) < 2 or len(target_gap_gene_sequence) > 4:
            gap_bias_penalty += config.distance_from_max_bridge_length_penalty * abs(len(target_gap_gene_sequence) - 4)
    penalties['gap_bias_penalty'] = gap_bias_penalty
    score += gap_bias_penalty

    # 7. GC Content
    if len(lhs_probe_sequence) == 0:
        lhs_gc = 0.5
    else:
        lhs_gc = (lhs_probe_sequence.upper().count("G") + lhs_probe_sequence.upper().count("C")) / len(
            lhs_probe_sequence)

    if len(rhs_probe_sequence) == 0:
        rhs_gc = 0.5
    else:
        rhs_gc = (rhs_probe_sequence.upper().count("G") + rhs_probe_sequence.upper().count("C")) / len(
            rhs_probe_sequence)

    gc_penalty = 0.0
    if config.strict_gc_content:
        if lhs_gc < config.tx_min_gc or lhs_gc > config.tx_max_gc:
            gc_penalty = config.invalid_score
        if rhs_gc < config.tx_min_gc or rhs_gc > config.tx_max_gc:
            gc_penalty = config.invalid_score
    elif config.penalize_gc_content:
        if lhs_gc < config.tx_min_gc: gc_penalty += 100 * config.scaled_gc_penalty * (config.tx_min_gc - lhs_gc)
        if lhs_gc > config.tx_max_gc: gc_penalty += 100 * config.scaled_gc_penalty * (lhs_gc - config.tx_max_gc)
        if rhs_gc < config.tx_min_gc: gc_penalty += 100 * config.scaled_gc_penalty * (config.tx_min_gc - rhs_gc)
        if rhs_gc > config.tx_max_gc: gc_penalty += 100 * config.scaled_gc_penalty * (rhs_gc - config.tx_max_gc)
    else:
        if lhs_gc < config.tx_min_gc or lhs_gc > config.tx_max_gc: gc_penalty += config.lenient_gc_penalty
        if rhs_gc < config.tx_min_gc or rhs_gc > config.tx_max_gc: gc_penalty += config.lenient_gc_penalty

    penalties['gc_content_penalty'] = gc_penalty
    score += gc_penalty

    # 8. Unbalanced GC
    unbalanced_penalty = config.unbalanced_gc_penalty * (abs(lhs_gc - rhs_gc) * 100)
    penalties['unbalanced_gc_penalty'] = unbalanced_penalty
    score += unbalanced_penalty

    # 9. Homopolymer
    homopolymer_score = 0.0
    homopolymer_score += config.homopolymer_penalty * max_homopolymer_length(lhs_probe_sequence)
    homopolymer_score += config.homopolymer_penalty * max_homopolymer_length(rhs_probe_sequence)
    penalties['homopolymer_penalty'] = homopolymer_score
    score += homopolymer_score

    # 10. Melting Temp
    tm_score = 0.0
    lhs_tm = get_melting_temp(lhs_probe_sequence)
    rhs_tm = get_melting_temp(rhs_probe_sequence)
    tm_score += (config.lhs_tm_penalty * lhs_tm)
    tm_score += (config.rhs_tm_penalty * rhs_tm)
    penalties['tm_penalty'] = tm_score
    score += tm_score

    # 11. Tandem Repeats
    repeat_score = 0.0
    if not is_gapfill:
        repeat_score += (config.tandem_repeat_penalty * has_tandem_repeat(lhs_probe_sequence))
        repeat_score += (config.tandem_repeat_penalty * has_tandem_repeat(rhs_probe_sequence))
    penalties['tandem_repeat_penalty'] = repeat_score
    score += repeat_score

    penalties['final_score'] = score
    return penalties


def main():
    parser = argparse.ArgumentParser(description="Calculate detailed score breakdown for a specific probe pair.")
    parser.add_argument("config_file", type=Path, help="Path to config.json")
    parser.add_argument("full_sequence", type=str,
                        help="Full sequence containing LHS and RHS (Either Probe Sense OR Gene Sense).")
    parser.add_argument("lhs_probe", type=str, help="LHS Probe sequence (Probe Sense).")
    parser.add_argument("rhs_probe", type=str, help="RHS Probe sequence (Probe Sense).")

    args = parser.parse_args()

    # 1. Load Config
    with open(args.config_file, 'r') as f:
        config_data = json.load(f)
    config = FlexProbeConfig.from_dict(config_data)

    full_seq = args.full_sequence.strip()
    lhs_probe = args.lhs_probe.strip()
    rhs_probe = args.rhs_probe.strip()

    # 2. Detect Orientation
    lhs_in_full = full_seq.find(lhs_probe)
    rhs_in_full = full_seq.find(rhs_probe)

    # Check Gene Sense (Target RNA)
    rc_lhs = reverse_complement(lhs_probe)
    rc_rhs = reverse_complement(rhs_probe)

    rc_lhs_in_full = full_seq.find(rc_lhs)
    rc_rhs_in_full = full_seq.find(rc_rhs)

    gap_gene = None
    lhs_gene = None
    rhs_gene = None
    detected_mode = "Unknown"
    gap_len = 0

    if lhs_in_full != -1 and rhs_in_full != -1 and lhs_in_full < rhs_in_full:
        detected_mode = "Probe Sense (LHS -> RHS)"
        print(f"Orientation Detected: {detected_mode}")

        # Gap is between LHS end and RHS start
        gap_start = lhs_in_full + len(lhs_probe)
        gap_end = rhs_in_full
        gap_probe = full_seq[gap_start:gap_end]

        # Convert to Gene Sense for Scoring
        lhs_gene = reverse_complement(lhs_probe)
        rhs_gene = reverse_complement(rhs_probe)

        # Handle 0bp logic
        gap_len = len(gap_probe)
        if gap_len > 0:
            gap_gene = reverse_complement(gap_probe)
        else:
            gap_gene = None  # Signal to scorer that this is adjacent/non-gapfill

    elif rc_lhs_in_full != -1 and rc_rhs_in_full != -1 and rc_rhs_in_full < rc_lhs_in_full:
        detected_mode = "Gene Sense (RC(RHS) -> RC(LHS))"
        print(f"Orientation Detected: {detected_mode}")

        # In Gene Sense, RHS (binds 5') comes BEFORE LHS (binds 3')
        rhs_gene = rc_rhs
        lhs_gene = rc_lhs

        # Gap is between RHS end and LHS start
        gap_start = rc_rhs_in_full + len(rc_rhs)
        gap_end = rc_lhs_in_full
        gap_probe_gene = full_seq[gap_start:gap_end]

        # Handle 0bp logic
        gap_len = len(gap_probe_gene)
        if gap_len > 0:
            gap_gene = gap_probe_gene
        else:
            gap_gene = None  # Signal to scorer that this is adjacent/non-gapfill

    else:
        print("Error: Could not determine orientation.")
        sys.exit(1)

    # 3. Calculate Score Breakdown
    results = score_probe_with_breakdown(config, lhs_gene, rhs_gene, gap_gene)

    # 4. Report
    print("\n--- Scoring Breakdown ---")
    print(f"LHS Probe ({len(lhs_probe)}bp): {lhs_probe}")
    print(f"RHS Probe ({len(rhs_probe)}bp): {rhs_probe}")
    print(f"Gap: {gap_len}bp")
    print(f"Junction: {results.get('junction_bases', 'N/A')}")
    print("-" * 30)

    # Sort by penalty magnitude (descending)
    sorted_penalties = sorted(results.items(),
                              key=lambda item: abs(item[1]) if isinstance(item[1], (int, float)) else 0, reverse=True)

    for key, val in sorted_penalties:
        if key == 'final_score': continue
        if key == 'junction_bases': continue
        print(f"{key:<35}: {val:>15.4f}")

    print("-" * 30)
    print(f"{'FINAL SCORE':<35}: {results['final_score']:>15.4f}")


if __name__ == "__main__":
    main()