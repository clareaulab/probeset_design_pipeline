import argparse
import sys
sys.path.insert(0, 'src')
import json
from pathlib import Path

import pandas as pd
from flex_probe_pipeline import build_genome, FlexProbeDesigner, FlexProbeConfig, SnvProbeHelper, HumanBackgroundFlexProbeConfig


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def find_probe_position(probe: str, transcript: str) -> tuple[int, bool] | None:
    """Find the position of a probe in the transcript.

    Returns (position, is_reverse_complement) or None if not found.
    Position is 0-indexed start position.
    """
    # Try direct match first
    pos = transcript.find(probe)
    if pos != -1:
        return pos, False

    # Try reverse complement
    rc_probe = reverse_complement(probe)
    pos = transcript.find(rc_probe)
    if pos != -1:
        return pos, True

    return None


def visualize_transcript_alignment(
        transcript: str,
        new_lhs: str, new_rhs: str,
        old_lhs: str, old_rhs: str,
        context: int = 2
) -> str:
    """Visualize where both old and new probes align to the transcript."""
    lines = []

    # Find positions of all probes
    new_lhs_pos = find_probe_position(new_lhs, transcript)
    new_rhs_pos = find_probe_position(new_rhs, transcript)
    old_lhs_pos = find_probe_position(old_lhs, transcript)
    old_rhs_pos = find_probe_position(old_rhs, transcript)

    def format_pos(pos_result, probe, label):
        if pos_result is None:
            return f"  {label}: NOT FOUND in transcript"
        pos, is_rc = pos_result
        rc_marker = " (RC)" if is_rc else ""
        return f"  {label}: pos {pos}-{pos + len(probe) - 1}{rc_marker}"

    lines.append("Transcript positions:")
    lines.append(format_pos(new_lhs_pos, new_lhs, "New LHS"))
    lines.append(format_pos(new_rhs_pos, new_rhs, "New RHS"))
    lines.append(format_pos(old_lhs_pos, old_lhs, "Old LHS"))
    lines.append(format_pos(old_rhs_pos, old_rhs, "Old RHS"))

    # If we found positions, show the relevant transcript region
    positions = []
    for pos_result, probe in [(new_lhs_pos, new_lhs), (new_rhs_pos, new_rhs),
                               (old_lhs_pos, old_lhs), (old_rhs_pos, old_rhs)]:
        if pos_result is not None:
            pos, _ = pos_result
            positions.append((pos, pos + len(probe)))

    if positions:
        min_pos = max(0, min(p[0] for p in positions) - context)
        max_pos = min(len(transcript), max(p[1] for p in positions) + context)

        lines.append(f"\nTranscript region [{min_pos}:{max_pos}]:")
        lines.append(f"  {transcript[min_pos:max_pos]}")

        # Create position markers
        def make_marker(pos_result, probe, char):
            if pos_result is None:
                return None
            pos, is_rc = pos_result
            marker = [' '] * (max_pos - min_pos)
            display_probe = reverse_complement(probe) if is_rc else probe
            for i, base in enumerate(display_probe):
                if min_pos <= pos + i < max_pos:
                    marker[pos + i - min_pos] = char
            return ''.join(marker)

        new_lhs_marker = make_marker(new_lhs_pos, new_lhs, 'N')
        new_rhs_marker = make_marker(new_rhs_pos, new_rhs, 'n')
        old_lhs_marker = make_marker(old_lhs_pos, old_lhs, 'O')
        old_rhs_marker = make_marker(old_rhs_pos, old_rhs, 'o')

        if new_lhs_marker:
            lines.append(f"  {new_lhs_marker}  <- New LHS (N)")
        if new_rhs_marker:
            lines.append(f"  {new_rhs_marker}  <- New RHS (n)")
        if old_lhs_marker:
            lines.append(f"  {old_lhs_marker}  <- Old LHS (O)")
        if old_rhs_marker:
            lines.append(f"  {old_rhs_marker}  <- Old RHS (o)")

    return '\n'.join(lines)


def find_best_alignment(a: str, b: str) -> tuple[int, int]:
    """Find the optimal shift and corresponding distance between two sequences.

    Returns (best_shift, best_dist) where:
    - positive shift means 'a' starts after 'b'
    - negative shift means 'b' starts after 'a'
    """
    max_shift = abs(len(a) - len(b)) + min(len(a), len(b))
    best_shift = 0
    best_dist = None
    for shift in range(-max_shift, max_shift + 1):
        if shift >= 0:
            a_slice = a[shift:]
            b_slice = b[:len(a_slice)]
            shift_penalty = shift
        else:
            b_slice = b[-shift:]
            a_slice = a[:len(b_slice)]
            shift_penalty = -shift
        mismatches = sum(1 for x, y in zip(a_slice, b_slice) if x != y)
        dist = shift_penalty + mismatches
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_shift = shift
    return best_shift, best_dist if best_dist is not None else 0


def sliding_distance(a: str, b: str) -> int:
    _, dist = find_best_alignment(a, b)
    return dist


def align_sequences(a: str, b: str) -> str:
    """Create a visual alignment showing matches between two sequences."""
    if a == b:
        return f"{a}\n{'|' * len(a)}\n{b}"

    best_shift, _ = find_best_alignment(a, b)

    if best_shift >= 0:
        # 'a' starts after 'b' by best_shift positions
        aligned_a = ' ' * best_shift + a
        aligned_b = b + ' ' * max(0, best_shift + len(a) - len(b))
    else:
        # 'b' starts after 'a' by -best_shift positions
        aligned_a = a + ' ' * max(0, -best_shift + len(b) - len(a))
        aligned_b = ' ' * (-best_shift) + b

    # Pad to same length
    max_len = max(len(aligned_a), len(aligned_b))
    aligned_a = aligned_a.ljust(max_len)
    aligned_b = aligned_b.ljust(max_len)

    # Build match line: '|' for match, ' ' for mismatch/gap
    match_line = ''.join(
        '|' if (ca != ' ' and cb != ' ' and ca == cb) else ' '
        for ca, cb in zip(aligned_a, aligned_b)
    )

    return f"{aligned_a}\n{match_line}\n{aligned_b}"


def create_designer_with_config(config_path: str) -> FlexProbeDesigner:
    """Create a FlexProbeDesigner with config from JSON file."""
    with open(config_path, 'r') as f:
        config_data = json.load(f)

    config = FlexProbeConfig(**config_data)
    designer = FlexProbeDesigner(
        working_dir="./",
        reference_probe_set="human_flex_v1",
        config=config
    )
    return designer


def score_probe_from_sequences(
        designer: FlexProbeDesigner,
        lhs_probe: str,
        rhs_probe: str,
        gap_gene_sequence: str = None,
        target_gap_start: int = None,
        target_gap_end: int = None
) -> float:
    """Score a probe using the designer's scoring function.

    Takes probe sequences (not gene sequences) and converts them.
    """
    # Convert probe sequences to gene sequences (reverse complement)
    lhs_gene_seq = reverse_complement(lhs_probe)
    rhs_gene_seq = reverse_complement(rhs_probe)

    return designer.score_probe(
        expect_hits=1,
        lhs_gene_sequence=lhs_gene_seq,
        rhs_gene_sequence=rhs_gene_seq,
        target_gap_gene_sequence=gap_gene_sequence,
        gap_target_start_position=target_gap_start,
        gap_target_end_position=target_gap_end
    )


def get_detailed_probe_score(
        designer: FlexProbeDesigner,
        lhs_probe: str,
        rhs_probe: str,
        gap_gene_sequence: str = None,
        target_gap_start: int = None,
        target_gap_end: int = None
) -> dict:
    """Get detailed scoring breakdown for a probe.

    Returns a dictionary with individual score components.
    """
    lhs_gene_seq = reverse_complement(lhs_probe)
    rhs_gene_seq = reverse_complement(rhs_probe)

    config = designer.config

    # Calculate individual components
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq

    details = {}

    # GC content
    lhs_gc = (lhs_probe.count('G') + lhs_probe.count('C')) / len(lhs_probe)
    rhs_gc = (rhs_probe.count('G') + rhs_probe.count('C')) / len(rhs_probe)
    details['lhs_gc'] = lhs_gc
    details['rhs_gc'] = rhs_gc
    details['gc_diff'] = abs(lhs_gc - rhs_gc)

    # Check GC bounds
    details['lhs_gc_valid'] = config.tx_min_gc <= lhs_gc <= config.tx_max_gc
    details['rhs_gc_valid'] = config.tx_min_gc <= rhs_gc <= config.tx_max_gc

    # Melting temperature
    try:
        lhs_tm = mt.Tm_NN(Seq(lhs_probe))
        rhs_tm = mt.Tm_NN(Seq(rhs_probe))
        details['lhs_tm'] = lhs_tm
        details['rhs_tm'] = rhs_tm
    except:
        details['lhs_tm'] = None
        details['rhs_tm'] = None

    # Gap length
    if gap_gene_sequence:
        details['gap_length'] = len(gap_gene_sequence)
        details['gap_valid'] = config.min_bridge_length <= len(gap_gene_sequence) <= config.max_bridge_length
    else:
        details['gap_length'] = 0
        details['gap_valid'] = True

    # Homopolymer check
    def has_homopolymer(seq, length=4):
        for base in 'ATGC':
            if base * length in seq:
                return True
        return False

    details['lhs_homopolymer'] = has_homopolymer(lhs_probe)
    details['rhs_homopolymer'] = has_homopolymer(rhs_probe)

    # Junction chemistry (check last base of LHS and first base of RHS)
    details['lhs_last_base'] = lhs_probe[-1]
    details['rhs_first_base'] = rhs_probe[0]

    # Total score
    details['total_score'] = score_probe_from_sequences(
        designer, lhs_probe, rhs_probe, gap_gene_sequence, target_gap_start, target_gap_end
    )

    return details


def score_old_probe_on_new_transcript(
        designer: FlexProbeDesigner,
        old_lhs: str,
        old_rhs: str,
        new_transcript: str,
        variant_pos: int = None
) -> dict:
    """Score the old probe sequences on the new transcript.

    This helps understand why the old probes might not have been selected.
    """
    result = {
        'old_lhs': old_lhs,
        'old_rhs': old_rhs,
        'found_in_transcript': False,
        'issues': []
    }

    # Find old probes in the new transcript
    old_lhs_pos = find_probe_position(old_lhs, new_transcript)
    old_rhs_pos = find_probe_position(old_rhs, new_transcript)

    if old_lhs_pos is None:
        result['issues'].append(f"Old LHS probe NOT FOUND in transcript")
    else:
        result['old_lhs_pos'] = old_lhs_pos[0]
        result['old_lhs_rc'] = old_lhs_pos[1]

    if old_rhs_pos is None:
        result['issues'].append(f"Old RHS probe NOT FOUND in transcript")
    else:
        result['old_rhs_pos'] = old_rhs_pos[0]
        result['old_rhs_rc'] = old_rhs_pos[1]

    if old_lhs_pos and old_rhs_pos:
        result['found_in_transcript'] = True

        # Calculate gap between probes
        lhs_end = old_lhs_pos[0] + len(old_lhs)
        rhs_start = old_rhs_pos[0]

        # Probes are on reverse complement, so positions are reversed
        if old_lhs_pos[1]:  # RC orientation
            gap_start = old_rhs_pos[0] + len(old_rhs)
            gap_end = old_lhs_pos[0]
        else:
            gap_start = lhs_end
            gap_end = rhs_start

        result['calculated_gap_start'] = gap_start
        result['calculated_gap_end'] = gap_end
        result['calculated_gap_length'] = abs(gap_end - gap_start)

        if gap_start < gap_end:
            gap_seq = new_transcript[gap_start:gap_end]
        else:
            gap_seq = new_transcript[gap_end:gap_start]
        result['calculated_gap_seq'] = gap_seq

        # Score the old probes
        try:
            result['old_probe_details'] = get_detailed_probe_score(
                designer, old_lhs, old_rhs, gap_seq,
                min(gap_start, gap_end), max(gap_start, gap_end)
            )
        except Exception as e:
            result['issues'].append(f"Could not score old probe: {e}")

        # Check if gap is valid
        config = designer.config
        if result['calculated_gap_length'] < config.min_bridge_length:
            result['issues'].append(f"Gap too short: {result['calculated_gap_length']} < {config.min_bridge_length}")
        if result['calculated_gap_length'] > config.max_bridge_length:
            result['issues'].append(f"Gap too long: {result['calculated_gap_length']} > {config.max_bridge_length}")

    return result


def compare_probes(
        new_probe_lhs: str,
        new_probe_rhs: str,
        old_probe_lhs: str,
        old_probe_rhs: str
):
    if new_probe_lhs == old_probe_lhs:
        lhs_dist = 0
    else:
        lhs_dist = sliding_distance(new_probe_lhs, old_probe_lhs)

    if new_probe_rhs == old_probe_rhs:
        rhs_dist = 0
    else:
        rhs_dist = sliding_distance(new_probe_rhs, old_probe_rhs)

    # Format visual alignment
    lhs_alignment = align_sequences(new_probe_lhs, old_probe_lhs)
    rhs_alignment = align_sequences(new_probe_rhs, old_probe_rhs)

    return lhs_dist, rhs_dist, lhs_alignment, rhs_alignment


def find_best_transcript_for_probes(
        gene_name: str,
        old_lhs: str,
        old_rhs: str,
        new_lhs: str,
        new_rhs: str,
        new_transcript: str,
        genome=None,
        ensembl_release: int = 111
) -> list[tuple[str, int, int, int]]:
    """Find which transcript(s) would produce probes overlapping with new probes.

    Compares old probe positions on candidate transcripts with new probe positions
    on the new transcript to find the best matching transcript.

    Returns list of (transcript_id, lhs_pos, rhs_pos, overlap_score) sorted by best overlap.
    """
    if genome is None:
        genome = build_genome(ensembl_release)

    # Find where new probes are on the new transcript
    new_lhs_pos = find_probe_position(new_lhs, new_transcript)
    new_rhs_pos = find_probe_position(new_rhs, new_transcript)

    new_lhs_start = new_lhs_pos[0] if new_lhs_pos else -1
    new_rhs_start = new_rhs_pos[0] if new_rhs_pos else -1

    # Get all transcripts for the gene
    try:
        genes = genome.genes_by_name(gene_name)
    except ValueError:
        print(f"Gene {gene_name} not found")
        return []

    results = []
    for gene in genes:
        for transcript in gene.transcripts:
            try:
                seq = transcript.coding_sequence
                if seq is None:
                    continue
            except (ValueError, AttributeError):
                # Some transcripts don't have coding sequences
                continue

            old_lhs_pos = find_probe_position(old_lhs, seq)
            old_rhs_pos = find_probe_position(old_rhs, seq)

            old_lhs_start = old_lhs_pos[0] if old_lhs_pos else -1
            old_rhs_start = old_rhs_pos[0] if old_rhs_pos else -1

            # Calculate overlap score (lower = better)
            # Check if old probe positions overlap with new probe positions
            if old_lhs_start >= 0 and new_lhs_start >= 0:
                lhs_overlap = abs(old_lhs_start - new_lhs_start)
            else:
                lhs_overlap = 9999

            if old_rhs_start >= 0 and new_rhs_start >= 0:
                rhs_overlap = abs(old_rhs_start - new_rhs_start)
            else:
                rhs_overlap = 9999

            # Total overlap score
            overlap_score = lhs_overlap + rhs_overlap

            # Found score: 0 = both found, 1 = one found, 2 = neither found
            found_score = (0 if old_lhs_pos else 1) + (0 if old_rhs_pos else 1)

            results.append((
                transcript.transcript_id,
                old_lhs_start,
                old_rhs_start,
                overlap_score,
                found_score,
                transcript.support_level or 99,
                len(seq)
            ))

    # Sort by: found_score, overlap_score, support level, length
    results.sort(key=lambda x: (x[4], x[3], x[5], -x[6]))

    print(f"\nTranscripts for {gene_name} (new LHS@{new_lhs_start}, new RHS@{new_rhs_start}):")
    print(f"{'Transcript ID':<20} {'Old LHS':>10} {'Old RHS':>10} {'Overlap':>10} {'Support':>8}")
    print("-" * 65)
    for tid, lhs, rhs, overlap, found, support, length in results[:10]:  # Show top 10
        lhs_str = str(lhs) if lhs >= 0 else "N/A"
        rhs_str = str(rhs) if rhs >= 0 else "N/A"
        overlap_str = str(overlap) if overlap < 9999 else "N/A"
        print(f"{tid:<20} {lhs_str:>10} {rhs_str:>10} {overlap_str:>10} {support:>8}")

    return [(tid, lhs, rhs, overlap) for tid, lhs, rhs, overlap, _, _, _ in results]


def cleanup_names(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
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


def test_probe_results(
        old_probes_file: str,
        new_probes_file: str,
        config_file: str = None,
        output_file: str = "probe_discrepancies.txt"
):
    """Compare two probe files and report discrepancies.

    Args:
        old_probes_file: Path to the original/reference probe file (TSV)
        new_probes_file: Path to the new probe file to compare (TSV)
        config_file: Optional path to config JSON for probe scoring
        output_file: Path for output discrepancies file
    """
    orig_df = pd.read_table(old_probes_file)
    orig_df = cleanup_names(orig_df)
    new_df = pd.read_table(new_probes_file)

    # Join on cleaned names
    merged_df = pd.merge(orig_df, new_df, on='name', suffixes=('_old', '_new'))
    # Print any missing probes before comparing
    missing_probes = set(orig_df['name']) - set(new_df['name'])
    if missing_probes:
        print(f"Missing probes in new set: {missing_probes}")

    # Create designer with config for scoring if provided
    designer = create_designer_with_config(config_file) if config_file else None

    exact = 0
    discrepancies = []

    for _, row in merged_df.iterrows():
        lhs_dist, rhs_dist, lhs_alignment, rhs_alignment = compare_probes(
            row['lhs_probe_new'], row['rhs_probe_new'], row['lhs_probe_old'], row['rhs_probe_old']
        )
        if lhs_dist != 0 or rhs_dist != 0:
            output_lines = []
            output_lines.append(f"Discrepancy found for probe '{row['name']}':")

            # Score both probes if designer is available
            if designer:
                # New probe score - use gap_gene_sequence from the new probe data
                gap_gene_seq = row.get('gap_gene_sequence') if pd.notna(row.get('gap_gene_sequence')) else None
                target_start = int(row['target_gap_start']) if pd.notna(row.get('target_gap_start')) else None
                target_end = int(row['target_gap_end']) if pd.notna(row.get('target_gap_end')) else None

                new_score = score_probe_from_sequences(
                    designer,
                    row['lhs_probe_new'], row['rhs_probe_new'],
                    gap_gene_seq, target_start, target_end
                )

                # Score old probe on the NEW transcript to understand why it wasn't selected
                transcript = row.get('transcript') if pd.notna(row.get('transcript')) else None
                if transcript:
                    old_probe_analysis = score_old_probe_on_new_transcript(
                        designer,
                        row['lhs_probe_old'], row['rhs_probe_old'],
                        transcript,
                        target_start
                    )

                    if old_probe_analysis['found_in_transcript']:
                        old_details = old_probe_analysis.get('old_probe_details', {})
                        old_score = old_details.get('total_score', float('inf'))
                    else:
                        old_score = float('inf')

                    score_diff = new_score - old_score
                    if old_score == float('inf'):
                        score_better = "OLD NOT VALID ON NEW TRANSCRIPT"
                    else:
                        score_better = "NEW BETTER" if score_diff < 0 else ("OLD BETTER" if score_diff > 0 else "EQUAL")
                    output_lines.append(f"Scores: New={new_score:.2f}, Old={old_score:.2f} ({score_better}, diff={score_diff:.2f})")

                    # Show issues with old probe
                    if old_probe_analysis['issues']:
                        output_lines.append("Old probe issues:")
                        for issue in old_probe_analysis['issues']:
                            output_lines.append(f"  - {issue}")

                    # Show detailed comparison if both are valid
                    if old_probe_analysis['found_in_transcript'] and 'old_probe_details' in old_probe_analysis:
                        old_details = old_probe_analysis['old_probe_details']
                        new_details = get_detailed_probe_score(
                            designer, row['lhs_probe_new'], row['rhs_probe_new'],
                            gap_gene_seq, target_start, target_end
                        )

                        output_lines.append(f"GC content - New LHS: {new_details['lhs_gc']:.2f}, Old LHS: {old_details['lhs_gc']:.2f}")
                        output_lines.append(f"GC content - New RHS: {new_details['rhs_gc']:.2f}, Old RHS: {old_details['rhs_gc']:.2f}")
                        output_lines.append(f"Gap length - New: {new_details['gap_length']}, Old: {old_details['gap_length']}")

                        if not old_details.get('lhs_gc_valid') or not old_details.get('rhs_gc_valid'):
                            output_lines.append(f"  ** Old probe GC out of bounds [{designer.config.tx_min_gc}-{designer.config.tx_max_gc}]")
                        if not old_details.get('gap_valid'):
                            output_lines.append(f"  ** Old probe gap invalid [{designer.config.min_bridge_length}-{designer.config.max_bridge_length}]")
                else:
                    # No transcript available, use simple scoring
                    old_score = score_probe_from_sequences(
                        designer,
                        row['lhs_probe_old'], row['rhs_probe_old'],
                        gap_gene_seq, target_start, target_end
                    )

                    score_diff = new_score - old_score
                    score_better = "NEW BETTER" if score_diff < 0 else ("OLD BETTER" if score_diff > 0 else "EQUAL")
                    output_lines.append(f"Scores: New={new_score:.2f}, Old={old_score:.2f} ({score_better}, diff={score_diff:.2f})")

            output_lines.append(f"LHS distance: {lhs_dist}")
            output_lines.append(lhs_alignment)
            output_lines.append(f"RHS distance: {rhs_dist}")
            output_lines.append(rhs_alignment)

            # Show transcript alignment if transcript is available
            if 'transcript' in row and pd.notna(row['transcript']):
                output_lines.append("")
                output_lines.append(visualize_transcript_alignment(
                    row['transcript'],
                    row['lhs_probe_new'], row['rhs_probe_new'],
                    row['lhs_probe_old'], row['rhs_probe_old']
                ))
            output_lines.append("-" * 40)

            discrepancy_text = '\n'.join(output_lines)
            print(discrepancy_text)
            discrepancies.append(discrepancy_text)
        else:
            exact += 1

    print(f"Probes: {exact}/{len(merged_df)} probes match exactly.")

    # Write discrepancies to file
    if discrepancies or missing_probes:
        with open(output_file, "w") as f:
            f.write("> ")
            f.write(" ".join(sys.argv))
            if missing_probes:
                f.write(f"Missing probes in new set: {missing_probes}\n")
                f.write("=" * 60 + "\n\n")
            f.write('\n'.join(discrepancies))
            f.write(f"\n\nSummary: {exact}/{len(merged_df)} probes match exactly.\n")
            f.write(f"Discrepancies: {len(discrepancies)}\n")
            f.write(f"Missing: {len(missing_probes)}\n")
        print(f"\nDiscrepancies written to {output_file}")

    return merged_df, discrepancies


def find_transcripts_for_discrepancies(
        old_probes_file: str,
        new_probes_file: str,
        ensembl_release: int = 111
):
    """Analyze discrepancies and find which transcripts contain the old probes.

    Args:
        old_probes_file: Path to the original/reference probe file (TSV)
        new_probes_file: Path to the new probe file to compare (TSV)
        ensembl_release: Ensembl release version to use (default: 111)
    """
    orig_df = pd.read_table(old_probes_file)
    orig_df = cleanup_names(orig_df)
    new_df = pd.read_table(new_probes_file)

    merged_df = pd.merge(orig_df, new_df, on='name', suffixes=('_old', '_new'))

    # Build genome once
    genome = build_genome(ensembl_release)

    suggestions = []
    for _, row in merged_df.iterrows():
        lhs_dist, rhs_dist, _, _ = compare_probes(
            row['lhs_probe_new'], row['rhs_probe_new'],
            row['lhs_probe_old'], row['rhs_probe_old']
        )
        if lhs_dist != 0 or rhs_dist != 0:
            # Extract gene name from the probe name
            gene_name = row['name'].split()[0]

            print(f"\n{'=' * 60}")
            print(f"Analyzing: {row['name']} (LHS dist: {lhs_dist}, RHS dist: {rhs_dist})")

            results = find_best_transcript_for_probes(
                gene_name,
                row['lhs_probe_old'], row['rhs_probe_old'],
                row['lhs_probe_new'], row['rhs_probe_new'],
                row['transcript'],
                genome,
                ensembl_release
            )

            # results: (transcript_id, lhs_pos, rhs_pos, overlap_score)
            if results and results[0][1] >= 0 and results[0][2] >= 0:
                best_tid, _, _, overlap = results[0]
                # Extract variant from name (e.g., "EGFR c.2582T>G" -> "c.2582T>G")
                parts = row['name'].split(maxsplit=1)
                variant = parts[1] if len(parts) > 1 else "0bp"
                suggestions.append((row['name'], f"{best_tid} {variant}", overlap))
                print(f"\nSuggested input: {best_tid} {variant} (overlap: {overlap})")

    if suggestions:
        print(f"\n{'=' * 60}")
        print("SUGGESTED INPUTS (sorted by overlap, lower=better):")
        print("=" * 60)
        # Sort by overlap score
        suggestions.sort(key=lambda x: x[2])
        for orig_name, suggested, overlap in suggestions:
            print(f"{orig_name:40} -> {suggested:40} (overlap: {overlap})")


def test_all_transcripts_for_gene(
        gene_name: str,
        variant: str,
        old_lhs: str,
        old_rhs: str,
        config_path: str = None,
        genome=None,
        ensembl_release: int = 111
) -> list[tuple[str, str, str, int]]:
    """Test all transcripts for a gene and find which produces probes closest to old probes.

    Returns list of (transcript_id, generated_lhs, generated_rhs, distance) sorted by distance.
    """
    if genome is None:
        genome = build_genome(ensembl_release)

    # Load config
    if config_path:
        with open(config_path, 'r') as f:
            config_data = json.load(f)
    else:
        config_data = {}
    config, _ = HumanBackgroundFlexProbeConfig(ensembl_release=ensembl_release, **config_data)

    # Create probe helper
    snv_helper = SnvProbeHelper(config, genome)

    # Create designer
    designer = FlexProbeDesigner(
        working_dir="./",
        reference_probe_set="human_flex_v1",
        config=config
    )

    # Get all transcripts for the gene
    try:
        genes = genome.genes_by_name(gene_name)
    except ValueError:
        print(f"Gene {gene_name} not found")
        return []

    results = []
    for gene in genes:
        for transcript in gene.transcripts:
            try:
                seq = transcript.coding_sequence
                if seq is None:
                    continue
            except (ValueError, AttributeError):
                continue

            tid = transcript.transcript_id

            try:
                # Parse variant info
                snv_start, snv_end, snv_action, snv_data, _ = snv_helper.parse_snv_info(variant)

                # Get sequences for this transcript
                original_seq, mutated_seq, mut_start, mut_end = snv_helper.get_gene_sequence(
                    tid, snv_start, snv_end, snv_action, snv_data, None
                )

                # Generate probes
                transcripts_dict = {f"{gene_name} {variant}": (original_seq, mutated_seq)}
                target_starts = [(snv_start, mut_start)]
                target_ends = [(snv_end, mut_end)]

                probe_df = designer.generate_gapfilling_flex_probe_set_df(
                    transcripts_dict,
                    target_starts,
                    target_ends,
                    expect_hits=[1],
                    n_probes=1,
                    visium=False,
                    barcode=1
                )

                if len(probe_df) > 0:
                    new_lhs = probe_df.iloc[0]['lhs_probe']
                    new_rhs = probe_df.iloc[0]['rhs_probe']

                    # Calculate distance to old probes
                    lhs_dist = sliding_distance(new_lhs, old_lhs)
                    rhs_dist = sliding_distance(new_rhs, old_rhs)
                    total_dist = lhs_dist + rhs_dist

                    results.append((tid, new_lhs, new_rhs, total_dist, lhs_dist, rhs_dist))

            except Exception as e:
                # Skip transcripts that fail
                continue

    # Sort by total distance
    results.sort(key=lambda x: x[3])
    return results


def find_best_transcripts_for_all_discrepancies(
        old_probes_file: str,
        new_probes_file: str,
        config_file: str = None,
        ensembl_release: int = 111
):
    """Test all possible transcripts for each discrepancy to find the best match.

    Args:
        old_probes_file: Path to the original/reference probe file (TSV)
        new_probes_file: Path to the new probe file to compare (TSV)
        config_file: Optional path to config JSON for probe generation
        ensembl_release: Ensembl release version to use (default: 111)
    """
    orig_df = pd.read_table(old_probes_file)
    orig_df = cleanup_names(orig_df)
    new_df = pd.read_table(new_probes_file)

    merged_df = pd.merge(orig_df, new_df, on='name', suffixes=('_old', '_new'))

    # Build genome once
    genome = build_genome(ensembl_release)

    all_suggestions = []

    for _, row in merged_df.iterrows():
        lhs_dist, rhs_dist, _, _ = compare_probes(
            row['lhs_probe_new'], row['rhs_probe_new'],
            row['lhs_probe_old'], row['rhs_probe_old']
        )

        if lhs_dist != 0 or rhs_dist != 0:
            gene_name = row['name'].split()[0]
            parts = row['name'].split(maxsplit=1)
            variant = parts[1] if len(parts) > 1 else "0bp"

            print(f"\n{'=' * 70}")
            print(f"Testing all transcripts for: {row['name']}")
            print(f"Current distance: LHS={lhs_dist}, RHS={rhs_dist}, Total={lhs_dist + rhs_dist}")
            print(f"Old probes: LHS={row['lhs_probe_old']}, RHS={row['rhs_probe_old']}")
            print("-" * 70)

            results = test_all_transcripts_for_gene(
                gene_name, variant,
                row['lhs_probe_old'], row['rhs_probe_old'],
                config_path=config_file,
                genome=genome,
                ensembl_release=ensembl_release
            )

            if results:
                print(f"\n{'Transcript':<20} {'LHS Dist':>10} {'RHS Dist':>10} {'Total':>10} {'Match?':>10}")
                print("-" * 70)

                for tid, new_lhs, new_rhs, total, ld, rd in results[:10]:
                    match = "EXACT" if total == 0 else ""
                    print(f"{tid:<20} {ld:>10} {rd:>10} {total:>10} {match:>10}")

                best = results[0]
                if best[3] == 0:
                    print(f"\n*** EXACT MATCH FOUND: {best[0]} ***")
                elif best[3] < lhs_dist + rhs_dist:
                    print(f"\n*** BETTER MATCH: {best[0]} (dist={best[3]} vs current={lhs_dist + rhs_dist}) ***")

                all_suggestions.append((
                    row['name'],
                    best[0],  # transcript_id
                    variant,
                    best[3],  # total distance
                    lhs_dist + rhs_dist  # current distance
                ))
            else:
                print("No valid transcripts found")

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY - BEST TRANSCRIPT FOR EACH DISCREPANCY")
    print("=" * 70)
    print(f"{'Probe Name':<35} {'Best Transcript':<20} {'New Dist':>10} {'Old Dist':>10}")
    print("-" * 70)

    exact_matches = 0
    improvements = 0
    for name, tid, variant, new_dist, old_dist in sorted(all_suggestions, key=lambda x: x[3]):
        status = ""
        if new_dist == 0:
            status = "EXACT"
            exact_matches += 1
        elif new_dist < old_dist:
            status = "BETTER"
            improvements += 1
        print(f"{name:<35} {tid:<20} {new_dist:>10} {old_dist:>10} {status}")

    print(f"\nExact matches found: {exact_matches}")
    print(f"Improvements found: {improvements}")

    # Write suggestions to file
    with open("best_transcripts.txt", "w") as f:
        f.write("# Best transcript IDs for probe discrepancies\n")
        f.write("# Format: TRANSCRIPT_ID VARIANT\n\n")
        for name, tid, variant, new_dist, old_dist in sorted(all_suggestions, key=lambda x: x[3]):
            if new_dist <= old_dist:  # Only include if same or better
                f.write(f"{tid} {variant}\n")

    print(f"\nBest transcripts written to best_transcripts.txt")


def main(
        old_probes_file: str,
        new_probes_file: str,
        config_file: str = None,
        output_file: str = "probe_discrepancies.txt",
        find_transcripts: bool = False,
        find_best: bool = False,
        ensembl_release: int = 111
):
    """Compare probe files and optionally find best matching transcripts.

    Args:
        old_probes_file: Path to the original/reference probe file (TSV)
        new_probes_file: Path to the new probe file to compare (TSV)
        config_file: Optional path to config JSON for probe scoring
        output_file: Path for output discrepancies file
        find_transcripts: If True, analyze discrepancies to find matching transcripts
        find_best: If True, test ALL transcripts to find the closest match
        ensembl_release: Ensembl release version to use (default: 111)
    """
    test_probe_results(old_probes_file, new_probes_file, config_file, output_file)

    if find_transcripts:
        find_transcripts_for_discrepancies(old_probes_file, new_probes_file, ensembl_release)

    if find_best:
        find_best_transcripts_for_all_discrepancies(old_probes_file, new_probes_file, config_file, ensembl_release)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compare two probe files and report discrepancies."
    )
    parser.add_argument(
        "old_probes",
        type=str,
        help="Path to the original/reference probe file (TSV)"
    )
    parser.add_argument(
        "new_probes",
        type=str,
        help="Path to the new probe file to compare (TSV)"
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to config JSON for probe scoring"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="probe_discrepancies.txt",
        help="Path for output discrepancies file (default: probe_discrepancies.txt)"
    )
    parser.add_argument(
        "--find-transcripts",
        action="store_true",
        help="Analyze discrepancies to find which transcripts contain old probes"
    )
    parser.add_argument(
        "--find-best",
        action="store_true",
        help="Test ALL transcripts to find the closest match for each discrepancy"
    )
    parser.add_argument(
        "--release",
        type=int,
        default=111,
        help="Ensembl release version to use (default: 111)"
    )

    args = parser.parse_args()

    main(
        args.old_probes,
        args.new_probes,
        config_file=args.config,
        output_file=args.output,
        find_transcripts=args.find_transcripts,
        find_best=args.find_best,
        ensembl_release=args.release
    )
