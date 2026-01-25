import dataclasses
import functools
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, Literal
import os
import pickle
from pathlib import Path

os.environ["PYENSEMBL_CACHE_DIR"] = "./ensembl_cache"

import numpy as np
import pandas as pd
from Bio.Seq import reverse_complement
from pyensembl.species import Species, human, mouse
from pyensembl import EnsemblRelease, Transcript
from scipy.optimize import dual_annealing, brute, Bounds
from tqdm.auto import tqdm

from utils import max_homopolymer_length, blast_search, transcriptome, rank_and_filter_transcripts, get_melting_temp, \
    find_overlaps_with_flex, check_overlap_by_alignment, has_tandem_repeat, fetch_human_flex_v2_probeset, \
    fetch_mouse_visiumhd_probeset, fetch_mouse_flex_v1_probeset, fetch_mouse_flex_v2_probeset, \
    fetch_human_visiumhd_probeset, fetch_human_flex_v1_probeset

"""
Note: LHS refers to the left-hand side of the gene sequence, and RHS refers to the right-hand side of the gene sequence.
      Unless otherwise specified.
"""

# TODO: Allow for configuring adapters
# Probe ordering constants
rhs_order_prefix = "/5Phos/"
rhs_order_probe_barcode_bridge = "ACGCGGTTAGCACGTANN"
rhs_order_suffix = "CGGTCCTAGCAA"
rhs_order_suffix_visium = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

lhs_order_prefix = "CCTTGGCACCCGAGAATTCCA"
lhs_truseq_prefix = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
# From table 2: https://kb.10xgenomics.com/hc/en-us/articles/17623693026445-How-do-I-design-custom-probes-for-a-Single-Cell-Gene-Expression-Flex-i-e-Fixed-RNA-Profiling-for-multiplexed-samples-experiment
order_barcode_seqs = [
    "ACTTTAGG",
    "AACGGGAA",
    "AGTAGGCT",
    "ATGTTGAC",
    "ACAGACCT",
    "ATCCCAAC",
    "AAGTAGAG",
    "AGCTGTGA",
    "ACAGTCTG",
    "AGTGAGTG",
    "AGAGGCAA",
    "ACTACTCA",
    "ATACGTCA",
    "ATCATGTG",
    "AACGCCGA",
    "ATTCGGTT"
]


@dataclass(slots=True, kw_only=True)
class FlexProbeConfig:
    """
    A dataclass to store the configuration of the Flex Probe pipeline.
    """
    # The transcriptomic background to use for the analysis
    background: dict[str, str] = dataclasses.field(default_factory=dict)

    # Designer settings
    add_probes_to_blast: bool = True  # Do we add accepted probes to our database for off target hits?

    # Bridge settings
    max_bridge_length: int = 10
    min_bridge_length: int = 1

    # Scoring settings
    flexible_probe_length_range: int = 0  # If zero, we will not allow for flexible bridge lengths If non-zero we will test max_bridge_length +/- flexible_bridge_range
    low_complexity_penalty: float = 100.  # Penalize low complexity regions for probes
    bias_longer_gaps: bool = False  # If True, we will bias towards larger gaps
    distance_from_max_bridge_length_penalty: float = 5.  # Bias for larger gaps
    last_base_penalty: float = 1e6  # Bias away from the snv being at the end of the gap
    first_base_penalty: float = 250.  # Bias away from the snv being at the start of the gap
    strict_gc_content: bool = True  # If True, we will only allow probes with GC content between tx_min_gc and tx_max_gc
    penalize_gc_content: bool = False # if True, we will penalize probes based on how much they diverge from GC range (scaled_gc_penalty)
    lenient_gc_penalty: float = 100.  # Penalize GC content outside the bounds
    scaled_gc_penalty: float = 1. # Penalize for every percent outside gc range if penalize_gc_content
    unbalanced_gc_penalty: float = 1.  # Penalize probes with dissimilar GC content
    homopolymer_penalty: float = 2.  # Penalize probes with large homopolymers
    lhs_tm_penalty: float = 1. # Penalize probes with low melting temperature
    rhs_tm_penalty: float = -1. # Penalize probes with high rhs melting temperature
    offtarget_penalty_base: float = 150.  # Penalize probes with off-target hits by exponentiating this value
    offtarget_must_be_same_gene: bool = True  # If true, only penalize off target hits based on the hits that appear on both the LHS and RHS
    existing_probe_penalty: float = 100.  # Penalty per overlapping position with an existing probe
    max_probe_overlap: int = 0  # Maximum allowed overlap (in bp) between designed probes. Set to 0 to disallow any overlap.
    require_transcriptome_hit: bool = False  # If True, require probes to have a hit back to the transcriptome
    flex_overlap_penalty: float = 1e7  # Penalty for overlapping with a 10x flex probe
    tandem_repeat_penalty: float = 1e5 # Penalty for probe have >= 4 repeats of a sequence of at least 3bp
    exclude_probes: list[tuple[str, str]] = dataclasses.field(default_factory=list)  # Exclude probes that match these sequences
    # Constants
    lhs_probe_length: int = 25
    rhs_probe_length: int = 25
    probe_length_penalty: float = 0  # Penalize per bp that the probe is off from the standard probe length if using flexible probe lengths
    # Invalid scores should always yield no probe
    invalid_score: float = 1e9
    # If we can't match a T, use a fallback junction with this penalty This should always be worse than an optimal junction
    suboptimal_gap_penalty: float = 10
    # If we cannot get a T or a suboptimal junction, use this penalty This should always be worse than a suboptimal junction
    invalid_gap_penalty: float = 16
    tx_min_gc: float = 0.2
    tx_max_gc: float = 0.8
    evalue_cutoff: float = 1.0  # The evalue cutoff for blast hits
    flex_0bp_genotyping: bool = False  # If True, we will make 0bp genotyping probes
    variant_side: str = None  # Desired probe with the variant for flex_0bp_genotyping (lhs or rhs)

    # Parse from json dictionary
    @staticmethod
    def from_dict(config_dict: dict) -> 'FlexProbeConfig':
        """
        Create a FlexProbeConfig from a dictionary.
        :param config_dict: The dictionary to create the FlexProbeConfig from.
        :return: The FlexProbeConfig.
        """
        return FlexProbeConfig(**config_dict)

def build_genome(ensembl_release: int, species: Species = human) -> EnsemblRelease:
    """
    Build a genome dictionary from Ensembl.
    :param ensembl_release: The Ensembl release to use.
    :param species: The species to use.
    :return: The genome dictionary.
    """
    ensembl_genome = EnsemblRelease(release=ensembl_release, species=species)
    ensembl_genome.download(overwrite=False)
    ensembl_genome.index(overwrite=False)
    return ensembl_genome


def GenomeBackgroundFlexProbeConfig(*, ensembl_release: int = 111, species: Species = human, **kwargs) -> tuple[FlexProbeConfig, EnsemblRelease]:
    """
    Build a Flex Probe configuration for a genome as background.
    :param ensembl_release: The Ensembl release to use for the background. By default 111
    :param species: The species to use for the background. By default, humans.
    :param kwargs: The keyword arguments to pass to the FlexProbeConfig.
    :return: The final FlexProbeConfig.
    """
    # Build the initial config
    config = FlexProbeConfig(**kwargs)
    background = dict()
    if len(config.background) != 0:  # Add the human genome as background to existing background
        background = dict(config.background)
    else:
        background = dict()

    # Add the genome as background
    ensembl_genome = build_genome(ensembl_release, species)
    background |= transcriptome(ensembl_genome, filter_incomplete=False)  # Don't filter out incomplete transcripts
    config.background = background
    return config, ensembl_genome


def HumanBackgroundFlexProbeConfig(ensembl_release: int = 111, **kwargs) -> tuple[FlexProbeConfig, EnsemblRelease]:
    """
    Build a Flex Probe configuration for the human genome as background.
    """
    return GenomeBackgroundFlexProbeConfig(ensembl_release=ensembl_release, species=human, **kwargs)


def MouseBackgroundFlexProbeConfig(ensembl_release: int = 111, **kwargs) -> tuple[FlexProbeConfig, EnsemblRelease]:
    """
    Build a Flex Probe configuration for the mouse genome as background.
    """
    return GenomeBackgroundFlexProbeConfig(ensembl_release=ensembl_release, species=mouse, **kwargs)


def MskImpactSnvProbeHelper(config: FlexProbeConfig, genome: EnsemblRelease = None) -> 'SnvProbeHelper':
    """
    Create the SnvProbeHelper with isoform overrides from MSK IMPACT.
    """
    overrides = pd.read_csv(
        "https://raw.githubusercontent.com/genome-nexus/genome-nexus-importer/4c66ae73e1243ee473ce52322e6fe9cd37e753e3/data/common_input/isoform_overrides_oncokb_grch38.txt",
        sep="\t")
    overrides = {row['hugo_symbol']: row['enst_id'] for _, row in overrides.iterrows()}
    # Replace some outdated gene names
    return SnvProbeHelper(
        config=config,
        transcript_overrides=overrides,
        genome=genome,
    )


def load_mane_select_transcripts() -> dict[str, str]:
    """
    Load MANE Select transcript mappings from NCBI.
    MANE (Matched Annotation from NCBI and Ensembl) provides agreed-upon canonical transcripts.

    Returns a dict mapping gene symbol -> Ensembl transcript ID (without version).
    """
    mane_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz"
    try:
        mane_df = pd.read_csv(mane_url, sep="\t", compression="gzip")
        # Filter to MANE Select (not MANE Plus Clinical)
        mane_select = mane_df[mane_df['MANE_status'] == 'MANE Select']
        # Extract gene symbol and Ensembl transcript ID (remove version)
        overrides = {}
        for _, row in mane_select.iterrows():
            gene = row['symbol']
            enst = row['Ensembl_nuc'].split('.')[0]  # Remove version number
            overrides[gene] = enst
        return overrides
    except Exception as e:
        print(f"Warning: Could not load MANE Select data: {e}")
        return {}


def ManeSelectSnvProbeHelper(config: FlexProbeConfig, genome: EnsemblRelease = None) -> 'SnvProbeHelper':
    """
    Create the SnvProbeHelper with isoform overrides from MANE Select.
    MANE Select provides the agreed-upon canonical transcript between RefSeq and Ensembl.
    """
    overrides = load_mane_select_transcripts()
    return SnvProbeHelper(
        config=config,
        transcript_overrides=overrides,
        genome=genome,
    )


class SnvProbeHelper:
    """
    This class aims to assist in designing probes for specific SNVs.
    """

    def __init__(self,
                 config: FlexProbeConfig,
                 genome: EnsemblRelease = None,
                 transcript_overrides = None):
        """
        Initialize the SNV probe helper.
        :param config: The FlexProbeConfig to use for the SNV probe design.
        :param genome: The genome to use for the SNV probe design. If not specified, we will use the human genome.
        :param transcript_overrides: A dictionary of transcript overrides to use for the SNV probe design.
        """
        # if genome is None:
        #     genome = build_genome(111, human)
        self.genome = genome
        if transcript_overrides is None:
            self.overrides = dict()
        else:
            self.overrides = transcript_overrides
        self.config = config

    def create_probe_args(self,
                          gene_name: str,
                          dna_snv_string: str) -> tuple[tuple[str, int, Optional[int]], tuple[str, int, Optional[int]]]:
        """
        Generate arguments for creating SNV gap-filling probes.
        :param gene_name: The name of the gene.
        :param dna_snv_string: If specified, the DNA SNV to target. Ex. "c.1234G>T". If none, gets arguments without SNV.
        :return: Return a tuple of tuples. A tuple for the original and the mutated sequence which describes:
            The transcript to target, the SNV start, and the SNV end.
        """
        # hg38_location_string = hg38_location_string.replace("chr", "")
        # chromosome, hg38_start, hg38_end = parse_chromosome_location(hg38_location_string)
        snv_start, snv_end, snv_action, snv_data, mutated_snv_length = self.parse_snv_info(dna_snv_string)

        try:
            original_sequence, mutated_sequence, mutated_snv_start, mutated_snv_end = self.get_gene_sequence(
                gene_name, snv_start, snv_end, snv_action, snv_data)
        except ValueError as e:
            print(e)
            return None, None

        return (original_sequence, snv_start, snv_end), (mutated_sequence, mutated_snv_start, mutated_snv_end)

    def modify_gene_sequence(self,
                               transcript: Transcript | str,
                               snv_start: int = None,
                               snv_end: int = None,
                               snv_action: Literal["mutation", "inversion", "deletion", "insertion", "deletion-insertion", "duplication"] = None,
                               snv_data: str | int | None = None,
                               validate: bool = True
                               ) -> list[tuple[str, str, int, int]] | None:
        """
        Modify a gene sequence based on an SNV.
        :param transcript: The sequence to modify.
        :param snv_start: The starting position of the SNV.
        :param snv_end: The ending position of the SNV.
        :param snv_action: The action/type of the SNV.
        :param snv_data: The data for the SNV, if any.
        :param validate: If True, we will validate the SNV data against the sequence.
        :return: A list of tuples: The original sequence, the modified sequence, the start of the SNV, and the end of the SNV. Or None if the SNV is invalid.
            This list will either contain more than one element if either the LHS or RHS span an exon-exon junction.
        """
        if isinstance(transcript, str):
            original_sequence = transcript  # String, so don't check exons
        # elif snv_start is not None:  # We are targetting an snv, so we should check for exon junctions
        #     exon_lhs_check_range = self.config.lhs_probe_length + self.config.max_bridge_length + self.config.flexible_probe_length_range
        #     exon_rhs_check_range = self.config.rhs_probe_length + self.config.max_bridge_length + self.config.flexible_probe_length_range
        #     lhs_spans_junction = False
        #     rhs_spans_junction = False
        #
        #     # First get the exon that the SNV is in
        #     snv_exon = None
        #     snv_exon_i = None
        #     for i, exon in enumerate(transcript.exons):
        #         if exon.start <= snv_start <= exon.end:
        #             snv_exon = exon
        #             snv_exon_i = i
        #             break
        #
        #     if snv_exon is None:
        #         raise ValueError("SNV does not fall within an exon.")
        #
        #     # Check if the LHS spans an exon junction
        #     if snv_start - exon_lhs_check_range < snv_exon.start:
        #         lhs_spans_junction = True
        #     # Check if the RHS spans an exon junction
        #     if snv_end + exon_rhs_check_range > snv_exon.end:
        #         rhs_spans_junction = True
        #
        #     if (not lhs_spans_junction) and (not rhs_spans_junction):  # Same exon, so just continue the pipeline
        #         original_sequence = transcript.coding_sequence
        #     else:  # We have to generate sequence variants for each intron combination
        #         original_results = self.modify_gene_sequence(transcript.coding_sequence, snv_start, snv_end, snv_action, snv_data, validate=validate)
        #         if original_results is None:
        #             raise ValueError("Invalid SNV for " + transcript.gene_name + " at " + str(snv_start))
        #         exons = transcript.exon_intervals
        #         # Now inject the LHS intron sequence
        #         if lhs_spans_junction:
        #             lhs_intron_sequence = ""
        #             for i, exon in enumerate(exons):
        #                 if i == snv_exon_i - 1:
        #                     # lhs_intron_sequence += transcript.coding_sequence[exon[0]:exon[1] + exon_lhs_check_range]
        #                     genomic_start, genomic_end = transcript.offset_range(exon[0], exons[i + 1][0])
        #                     lhs_intron_sequence += transcript.gene
        #                 else:
        #                     lhs_intron_sequence += transcript.coding_sequence[exon[0]:exon[1]]
        else:
            original_sequence = transcript.coding_sequence

        if snv_action == '0bp':
            mutated_sequence = original_sequence
            mutated_snv_start = None
            mutated_snv_end = None
            return original_sequence, mutated_sequence, mutated_snv_start, mutated_snv_end

        if snv_start is not None:
            # Convert to 1 indexed
            snv_start_idx = snv_start - 1
            snv_end_idx = snv_end
            # print(snv_start_idx, snv_end_idx)

            if len(original_sequence) < snv_end_idx + self.config.rhs_probe_length:  # and validate:  <- Needed to prevent string index errors
                return None

            if snv_action == "mutation":
                from_nucleotide, to_nucleotide = snv_data.split(">")
                # assert from_nucleotide == original_sequence[snv_start], "SNV data does not match the original sequence."
                if from_nucleotide != original_sequence[snv_start_idx] and to_nucleotide != '*' and validate:
                    # print("SNV data does not match the original sequence.", from_nucleotide, to_nucleotide, original_sequence[snv_start_idx-1:snv_start_idx+2])
                    return None
                mutated_sequence = original_sequence[:snv_start_idx] + to_nucleotide + original_sequence[snv_end_idx:]
                mutated_snv_start = snv_start
                mutated_snv_end = snv_end
            elif snv_action == "inversion":
                conversion = snv_data.replace("inv", "")
                inversion_target_sequence = original_sequence[snv_start_idx:snv_end_idx]
                inverted_sequence = reverse_complement(inversion_target_sequence)
                # assert inverted_sequence == conversion, "SNV data does not match the inverted sequence."
                if inverted_sequence != conversion and validate:
                    # print("SNV data does not match the inverted sequence.")
                    return None
                mutated_sequence = original_sequence[:snv_start_idx] + inverted_sequence + original_sequence[snv_end_idx:]
                mutated_snv_start = snv_start
                mutated_snv_end = snv_end
            elif snv_action == "deletion":
                mutated_sequence = original_sequence[:snv_start_idx] + original_sequence[snv_end_idx:]
                mutated_snv_start = snv_start
                mutated_snv_end = snv_start - 1 ## in the mutated sequence, this could be a 0bp probe; target_length should be 0
            elif snv_action == "insertion":
                inserted_sequence = snv_data.replace("ins", "")
                mutated_sequence = original_sequence[:snv_start_idx] + inserted_sequence + original_sequence[snv_start_idx:]
                mutated_snv_start = snv_start
                mutated_snv_end = snv_end + len(inserted_sequence)
            elif snv_action == "deletion-insertion":
                inserted_sequence = snv_data.replace("delins", "")
                mutated_sequence = original_sequence[:snv_start_idx] + inserted_sequence + original_sequence[snv_end_idx:]
                mutated_snv_start = snv_start
                mutated_snv_end = snv_end + len(inserted_sequence) - (snv_end - snv_start) - 1
            elif snv_action == "duplication":
                duplicated = snv_data.replace("dup", "")
                # assert duplicated == original_sequence[snv_start:snv_end + 1], "SNV data does not match the original sequence."
                # Note that the duplication index is to the left of the duplicated nucleotide
                if not original_sequence[snv_start_idx:snv_end_idx].startswith(duplicated) and validate:
                    # print("SNV data does not match the original sequence.")
                    return None
                if duplicated == '':
                    # If the duplicated sequence is empty, we just insert the original sequence
                    duplicated = original_sequence[snv_start_idx:snv_end_idx]
                mutated_sequence = original_sequence[:snv_start_idx] + duplicated + original_sequence[snv_start_idx:]
                mutated_snv_start = snv_start
                mutated_snv_end = snv_end + len(duplicated)
            elif snv_action == "undefined":
                # We don't know what it is but we have a range to consider
                mutated_sequence = original_sequence
                mutated_snv_start = snv_start
                mutated_snv_end = snv_end
            else:
                return None
        else:
            mutated_sequence = original_sequence
            mutated_snv_start = snv_start
            mutated_snv_end = snv_end

        original_sequence = [char for char in original_sequence if char.isalpha()]
        original_sequence = "".join(original_sequence)
        mutated_sequence = [char for char in mutated_sequence if char.isalpha()]
        mutated_sequence = "".join(mutated_sequence)
        return original_sequence, mutated_sequence, mutated_snv_start, mutated_snv_end

    def get_gene_sequence(self,
                          gene_name: str,
                          snv_start: int = None,
                          snv_end: int = None,
                          snv_action: Literal["mutation", "inversion", "deletion", "insertion", "deletion-insertion", "duplication"] = None,
                          snv_data: str | int | None = None,
                          transcript_sequence: str | None = None,
                          validate: bool = False) -> list[tuple[str, str, int, int]]:
        """
        Get the sequence of a gene from the genome.
        :param gene_name: The name of the gene.
        :param snv_start: The starting position of the SNV. If None, the original sequence is returned.
        :param snv_end: The ending position of the SNV. If None, the original sequence is returned.
        :param snv_action: The action/type of the SNV. If None, the original sequence is returned.
        :param snv_data: The data for the SNV, if any. Example: the base that was inserted. If None, the original sequence is returned.
        :param transcript_sequence: The sequence of the transcript. If None, the sequence is fetched from the genome.
        :param validate: If True, we will validate the SNV data against the sequence.
        :return: Original gene sequence, the mutated gene sequence, and the region of interest within the mutated sequence.
        """
        if transcript_sequence is not None:
            results = self.modify_gene_sequence(transcript_sequence, snv_start, snv_end, snv_action, snv_data, validate=validate)
            if results is None:
                raise ValueError("Invalid SNV for " + gene_name + " at " + str(snv_start))
            return results

        if gene_name.startswith("ENST"):  # ENSEMBL ID
            # FIXME: Currently doesn't accept versions so we are assuming that hte latest versoin was provided
            # https://github.com/openvax/pyensembl/issues/242
            transcript = self.genome.transcript_by_id(gene_name.split(".")[0])
            results = self.modify_gene_sequence(transcript, snv_start, snv_end, snv_action, snv_data, validate=validate)
            if results is not None:
                return results
        else:
            gene_objs = self.genome.genes_by_name(gene_name)
            if len(gene_objs) != 1:
                print(f"WARNING: Found {len(gene_objs)} genes for {gene_name}")
            transcript_candidates = []
            for gene_obj in gene_objs:
                transcript_candidates.extend(rank_and_filter_transcripts(gene_obj.transcripts))
            # transcript_candidates = self.ensembl_genome.transcripts_at_locus(contig=chromosome, position=hg38_start, end=hg38_end)

            if len(transcript_candidates) == 0:
                raise ValueError("No transcripts found for the specified gene.")

            transcript = None
            # Prioritize msk canonical transcript, else use support level and length
            for candidate in sorted(transcript_candidates, key=lambda x: (1e8 if x.transcript_id in self.overrides.values() else 0, -(x.support_level or 0), ((x.biotype == 'protein_coding') + x.complete), x.length), reverse=True):  # Sort by support level and length
                # Skip transcripts with missing or malformed coding_sequence
                # Some pyensembl transcripts incorrectly include UTR in coding_sequence
                if candidate.coding_sequence is None or not candidate.coding_sequence.startswith('ATG'):
                    continue
                transcript = candidate

                results = self.modify_gene_sequence(transcript, snv_start, snv_end, snv_action, snv_data, validate=validate)
                if results is None:
                    continue
                else:
                    return results

        raise ValueError(f"No correct protein coding transcripts found for {gene_name} with SNV at {snv_start}.")

    def parse_snv_info(self, dna_snv_string: str) -> tuple[int, int, str, str, int]:
        """
        Parse the SNV information from the input strings.
        :param dna_snv_string: The DNA SNV string (ex. c.1849G>T).
        :return: The SNV start, SNV end, action, and data.
        """

        if dna_snv_string == '0bp':
            start = None
            end = None
            action = "0bp"
            data = "undefined"
            snv_length_change = 0
            return start, end, action, data, snv_length_change

        assert dna_snv_string.startswith("c."), "Invalid DNA SNV string."
        dna_snv_string = dna_snv_string[2:]
        # Classify the action first to determine how to parse the data
        if '>' in dna_snv_string:
            action = "mutation"
            split = dna_snv_string.split(">")
            start = int(split[0][:-1])  # Remove the base at the end
            end = start
            data = f"{split[0][-1]}>{split[1]}"  # recreate the part that was split
            snv_length_change = 0  # No change in region length
        elif 'delins' in dna_snv_string:
            action = "deletion-insertion"
            split = dna_snv_string.split("delins")
            dna_range = split[0].split("_")
            if len(dna_range) == 1:
                start = int(dna_range[0][:-1])
                end = start
            else:
                start = int(dna_range[0])
                end = int(dna_range[1])
            data = f"delins{split[1]}"
            snv_length_change = len(split[1]) - (end - start + 1)
        elif 'del' in dna_snv_string:
            action = "deletion"
            split = dna_snv_string.split("del")
            dna_range = split[0].split("_")
            if len(dna_range) == 1:
                start = int(dna_range[0])
                end = start
            else:
                start = int(dna_range[0])
                end = int(dna_range[1])
            data = "del"
            snv_length_change = -(end - start + 1)
        elif 'ins' in dna_snv_string:
            action = "insertion"
            split = dna_snv_string.split("ins")
            dna_range = split[0].split("_")
            if len(dna_range) == 1:
                start = int(dna_range[0])
                end = start
            else:
                start = int(dna_range[0])
                end = int(dna_range[1])
            data = f"ins{split[1]}"
            snv_length_change = len(split[1])
        elif 'dup' in dna_snv_string:
            action = "duplication"
            split = dna_snv_string.split("dup")
            dna_range = split[0].split("_")
            if len(dna_range) == 1:
                start = int(dna_range[0])
                end = start
            else:
                start = int(dna_range[0])
                end = int(dna_range[1])
            data = f"dup{split[1]}"
            snv_length_change = (end - start + 1)  # We are duplicating the region
        elif 'inv' in dna_snv_string:
            action = "inversion"
            split = dna_snv_string.split("inv")
            dna_range = split[0].split("_")
            if len(dna_range) == 1:
                start = int(dna_range[0])
                end = start
            else:
                start = int(dna_range[0])
                end = int(dna_range[1])
            data = f"inv{split[1]}"
            snv_length_change = 0
        else:
            print(f"WARNING: Unrecognized SNV action: {dna_snv_string}, assuming genomic range.")
            split = dna_snv_string.split("_")
            if len(split) == 1:
                start = int(split[0])
                end = start
            else:
                start = int(split[0])
                end = int(split[1])
            data = "undefined"
            action = "undefined"
            snv_length_change = 0

        return start, end, action, data, snv_length_change


@dataclass(slots=True, unsafe_hash=True)
class Probe:
    """
    A dataclass to store the information of a probe.
    Note that in this instance, we are considering LHS/RHS with respect to the PROBE, not the gene.
    """
    transcript_name: str  # Name
    # Probe-sense sequences
    rhs_probe: str
    lhs_probe: str
    # Gene-sense positions
    rhs_gene_start: int
    lhs_gene_start: int
    # Gene-sense sequences
    rhs_gene_sequence: str
    lhs_gene_sequence: str
    # GC content
    rhs_GC: float
    lhs_GC: float
    # BLAST hits
    expected_blast_hits: int
    total_rhs_blast_hits: int
    total_lhs_blast_hits: int
    rhs_blast_hits: list[str]
    lhs_blast_hits: list[str]
    # Overall score
    score: float
    # Transcript sequence
    transcript_sequence: str

    # Only for gap-filling probes
    gap_length: Optional[int] = None
    gap_probe_sequence: Optional[str] = None
    gap_gene_sequence: Optional[str] = None
    # Position of gap within the transcript
    target_start_gap: Optional[int] = None
    target_end_gap: Optional[int] = None
    # Additional info if available
    original_transcript_sequence: Optional[str] = None
    original_target_start_gap: Optional[int] = None
    original_target_end_gap: Optional[int] = None
    nonmutated_rhs_probe: Optional[int] = None
    nonmutated_lhs_probe: Optional[int] = None

    @property
    def is_gapfill(self) -> bool:
        return self.gap_length is not None

    # IDT compatible order format
    def rhs_probe_order(self, barcode: int = 0, visium: bool = False) -> str:
        if barcode > 0 and visium:
            raise ValueError("Visium probes do not support multiple barcodes.")
        if visium:
            return f"{rhs_order_prefix}{self.rhs_probe}{rhs_order_suffix_visium}"

        return f"{rhs_order_prefix}{self.rhs_probe}{rhs_order_probe_barcode_bridge}{order_barcode_seqs[barcode]}{rhs_order_suffix}"

    def lhs_probe_order(self, truseq: bool = False) -> str:
        return f"{lhs_truseq_prefix if truseq else lhs_order_prefix}{self.lhs_probe}"

    def probe_order(self, barcodes: int = 1, visium: bool = False) -> list[tuple[str, str]]:
        """
        Get the IDT compatible order for the probes.
        :param barcodes: The number of barcodes to order.
        :param visium: If True, we will use the Visium probe ordering.
        :return: The IDT compatible order (tuples of rhs and lhs probes).
        """
        idt_order = []
        for i in range(barcodes):
            idt_order.append((self.rhs_probe_order(i, visium=visium), self.lhs_probe_order()))
        return idt_order


class FlexProbeDesigner:
    """
    A class to configure and run the Flex Probe pipeline.
    """

    def __init__(self,
                 working_dir: str,
                 reference_probe_set: Literal['human_flex_v1', 'human_flex_v2', 'human_visiumhd',
                    'mouse_flex_v1', 'mouse_flex_v2', "mouse_visiumhd"] | pd.DataFrame,
                 config: FlexProbeConfig):
        self.working_dir = Path(working_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)
        self.config = config
        self.blast_db = self.working_dir / "blast_db" / "background"

        # We will cache the blast hits to avoid redundant calls
        self.blast_hits = functools.lru_cache(maxsize=None)(self.blast_hits)
        # We will also cache scores to avoid redundant calculations
        self.score_probe = functools.lru_cache(maxsize=None)(self.score_probe)

        if isinstance(reference_probe_set, pd.DataFrame):
            self.reference_probes = reference_probe_set.copy()
        else:
            match reference_probe_set:
                case 'human_flex_v2':
                    func = fetch_human_flex_v2_probeset
                case 'human_flex_v1':
                    func = fetch_human_flex_v1_probeset
                case 'human_visiumhd':
                    func = fetch_human_visiumhd_probeset
                case 'mouse_flex_v2':
                    func = fetch_mouse_flex_v2_probeset
                case 'mouse_flex_v1':
                    func = fetch_mouse_flex_v1_probeset
                case 'mouse_visiumhd':
                    func = fetch_mouse_visiumhd_probeset
                case _:
                    raise ValueError("Invalid reference probe set.")
            self.reference_probes = func()

        #print("Building the BLAST database...")
        #self.build_blast_db()

    def build_blast_db(self):
        """
        Build the blast database for the background.
        """
        if self.blast_db.parent.exists():
            shutil.rmtree(self.blast_db.parent)
        self.blast_db.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.NamedTemporaryFile(mode="w+", delete=True) as f:
            for name, sequence in self.config.background.items():
                f.write(f">{name}\n{sequence}\n")
            f.seek(0)

            subprocess.run(["makeblastdb", "-in", f.name, "-dbtype", "nucl", "-out", str(self.blast_db)], stdout=subprocess.DEVNULL)

    def append_to_background(self, seqs: dict[str, str]):
        """
        Append sequences to the background.
        :param seqs: The sequences to append.
        """
        if len(seqs) == 0:
            return

        # Add the sequences to the background
        for name, sequence in seqs.items():
            self.config.background[name] = sequence

        # Rebuild the database
        #self.build_blast_db()

    def blast_hits(self, sequence: str) -> tuple[int, list[str]]:
        """
        Blast a sequence against the background.
        :param sequence: The sequence to blast.
        :return: The number of hits and the names of the hits.
        """
        if len(sequence) == 0:
            return 0, list()

        blast_out = blast_search(
            blast_db=self.blast_db,
            sequence=sequence,
            evalue=self.config.evalue_cutoff
        )
        if blast_out is None or len(blast_out.alignments) == 0:
            return 0, list()

        # WARNING: We assume that all sequences have a non-unique identifier as the first word in the sequence title
        # I.e. gene name followed by a space followed by the unique transcript id

        # Build a tree of hits
        hits = defaultdict(list)
        for aln in blast_out.alignments:
            name = aln.hit_def.split(" ")[0]
            # 10X recommends at least 5 mismatches for a probe to not hybridize efficiently, so select alignments with < 5 mismatches to be "hits"
            matches = sum([hsp.identities for hsp in aln.hsps])
            if matches >= len(sequence) - 5:
                hits[name].append(aln)

        n_hits = len(hits)
        return n_hits, list(hits.keys())

    def score_probe(self,
                    expect_hits: int,
                    lhs_gene_sequence: str,
                    rhs_gene_sequence: str,
                    target_gap_gene_sequence: Optional[str] = None,
                    original_gap_gene_sequence: Optional[str] = None,
                    gap_target_start_position: Optional[int] = None,
                    gap_target_end_position: Optional[int] = None,
                    flex_overlap: Optional[int] = 0) -> float:
        """
        The scoring function for a probe based on heuristic rules. The lower the score, the better the probe.
        :param expect_hits: The number of base transcriptome hits to expect.
        :param lhs_gene_sequence: The left-hand side sequence of the gene.
        :param rhs_gene_sequence: The right-hand side sequence of the gene.
        :param target_gap_gene_sequence: The target gap sequence between the LHS and RHS gene sequences. If None, we assume this
            is not a gap-filling probe.
        :param gap_target_start_position: The start position of the targeted variant within the gap sequence.
        :param gap_target_end_position: The end position of the targeted variant within the gap sequence.
        :return: The score of the probe.
        """
        if target_gap_gene_sequence is not None:
            is_gapfill = True
        else:
            is_gapfill = False

        # Find the sequence in probe-space
        lhs_probe_sequence = reverse_complement(lhs_gene_sequence)
        rhs_probe_sequence = reverse_complement(rhs_gene_sequence)
        if is_gapfill:
            target_gap_probe_sequence = reverse_complement(target_gap_gene_sequence)
            if len(target_gap_probe_sequence) < self.config.min_bridge_length or len(target_gap_probe_sequence) > self.config.max_bridge_length:
                return self.config.invalid_score

        score = 0.0
        # Penalize the probe by how off the probe length is from optimal
        score += self.config.probe_length_penalty * abs(len(lhs_probe_sequence) - self.config.lhs_probe_length)
        score += self.config.probe_length_penalty * abs(len(rhs_probe_sequence) - self.config.rhs_probe_length)
        # Additionally penalize by imbalance in probe length
        score += self.config.probe_length_penalty * abs(len(lhs_probe_sequence) - len(rhs_probe_sequence))

        # Next, we should make sure the probes are unique to each other
        if lhs_probe_sequence == rhs_probe_sequence:
            return self.config.invalid_score

        # Note, the gap filling chemistry requires that the last base of the filled gap is T
        if is_gapfill:
            if len(rhs_probe_sequence) == 0:
                junction = ""
            else:
                junction = target_gap_probe_sequence[-1].upper() + rhs_probe_sequence[0].upper()
        else:
            # Non gapfill requires that the last base of the RHS is T
            junction = lhs_probe_sequence[-1].upper() + rhs_probe_sequence[0].upper()
        if junction not in ["TA", "TC", "TG", "TT", "CT", "CA", "AT"]:  # Empirically good junctions
            # Potential fallback junctions
            if junction in ["GA", "GT"]:  # Somewhat good junctions
                score += self.config.suboptimal_gap_penalty
            elif junction in ["AA", "AC", "AG"]:  # Other junctions are not good but not the worst.
                score += self.config.suboptimal_gap_penalty
            elif junction in ["CC", "GC", "CG", "GG"]:  # These were empirically the worst
                score += self.config.invalid_gap_penalty
            # If really bad junction, we need a small gap
            elif is_gapfill and len(target_gap_probe_sequence) > 1:
                score += self.config.invalid_gap_penalty
            else:
                return self.config.invalid_score

        # Penalize targetting annotated low complexity regions
        score += self.config.low_complexity_penalty * sum([c.islower() for c in lhs_gene_sequence])
        score += self.config.low_complexity_penalty * sum([c.islower() for c in rhs_gene_sequence])

        if is_gapfill:
            # Bias towards larger gaps
            if self.config.bias_longer_gaps:
                score += self.config.distance_from_max_bridge_length_penalty * abs(len(target_gap_gene_sequence) - self.config.max_bridge_length)
            # Bias toward 3-5bp
            elif len(target_gap_gene_sequence) < 2 or len(target_gap_gene_sequence) > 4: ## changed from 3 and 5
                score += self.config.distance_from_max_bridge_length_penalty * abs(len(target_gap_gene_sequence) - 4)
            # else:  # Bias towards middle size of gap
            #     score += self.config.distance_from_max_bridge_length_penalty * (len(target_gap_gene_sequence) - (self.config.max_bridge_length//2))
            # Bias away from the snv being at the end of the gap
            if gap_target_start_position is not None and (gap_target_start_position == 0):
                score += self.config.last_base_penalty
            if gap_target_end_position is not None and (gap_target_end_position == len(target_gap_gene_sequence)): ## this is actually the start position in probe space
                score += self.config.first_base_penalty


        # Next, test GC content
        if len(lhs_probe_sequence) == 0:
            lhs_gc = 0.5
        else:
            lhs_gc = (lhs_probe_sequence.upper().count("G") + lhs_probe_sequence.upper().count("C")) / len(lhs_probe_sequence)
        if len(rhs_probe_sequence) == 0:
            rhs_gc = 0.5
        else:
            rhs_gc = (rhs_probe_sequence.upper().count("G") + rhs_probe_sequence.upper().count("C")) / len(rhs_probe_sequence)

        # First check if the GC content is within the recommended bounds
        if self.config.strict_gc_content:  # We only allow probes with GC content between tx_min_gc and tx_max_gc
            if lhs_gc < self.config.tx_min_gc or lhs_gc > self.config.tx_max_gc:
                return self.config.invalid_score
            if rhs_gc < self.config.tx_min_gc or rhs_gc > self.config.tx_max_gc:
                return self.config.invalid_score
        elif self.config.penalize_gc_content:
            # Penalize GC content outside the bounds
            if lhs_gc < self.config.tx_min_gc:
                score += 100 * self.config.scaled_gc_penalty * (self.config.tx_min_gc - lhs_gc)
            if lhs_gc > self.config.tx_max_gc:
                score += 100 * self.config.scaled_gc_penalty * (lhs_gc - self.config.tx_max_gc)
            if rhs_gc < self.config.tx_min_gc:
                score += 100 * self.config.scaled_gc_penalty * (self.config.tx_min_gc - rhs_gc)
            if rhs_gc > self.config.tx_max_gc:
                score += 100 * self.config.scaled_gc_penalty * (rhs_gc - self.config.tx_max_gc)
        else:
            # Penalize GC content outside the bounds
            if lhs_gc < self.config.tx_min_gc:
                score += self.config.lenient_gc_penalty
            if lhs_gc > self.config.tx_max_gc:
                score += self.config.lenient_gc_penalty
            if rhs_gc < self.config.tx_min_gc:
                score += self.config.lenient_gc_penalty
            if rhs_gc > self.config.tx_max_gc:
                score += self.config.lenient_gc_penalty

        # Penalize probes with dissimilar GC content
        score += self.config.unbalanced_gc_penalty * (abs(lhs_gc - rhs_gc)*100)

        # Penalize probes with large homopolymers
        score += self.config.homopolymer_penalty * max_homopolymer_length(lhs_probe_sequence)
        score += self.config.homopolymer_penalty * max_homopolymer_length(rhs_probe_sequence)

        # Penalize probes by Tm
        lhs_tm = get_melting_temp(lhs_probe_sequence)
        rhs_tm = get_melting_temp(rhs_probe_sequence)
        score += (self.config.lhs_tm_penalty * lhs_tm)
        score += (self.config.rhs_tm_penalty * rhs_tm)

        score += (flex_overlap * self.config.flex_overlap_penalty)

        ## penalize probes with tandem repeats; only for 0bp probes that can move around
        if not is_gapfill:
            score += (self.config.tandem_repeat_penalty * has_tandem_repeat(lhs_probe_sequence))
            score += (self.config.tandem_repeat_penalty * has_tandem_repeat(rhs_probe_sequence))

        return score

    def _test_flex_probe(self,
                         args: np.ndarray,
                         transcript_name: str,
                         transcript_sequence: str,
                         original_transcript_sequence: str,
                         expect_hits: int,
                         is_gapfill: bool,
                         target_position: int,
                         target_length: int,
                         carryover_visits: dict[tuple[int, int, int, int], float],
                         existing_probes: dict[str, list[Probe]],
                         flex_probes: dict[str, list[Probe]],
                         original_target_start: int,
                         original_target_end: int,
                         predefined_lhs: list[str]
                         ) -> float:
        """Wrapper for dual annealing to test flex probes."""
        start_idx = int(args[0])
        if is_gapfill:
            has_length_adjustment = len(args) == 4
            gap_length = int(args[-1])
        else:
            has_length_adjustment = len(args) == 3
            gap_length = None
        if has_length_adjustment:
            lhs_probe_length_adjustment = int(args[1])
            rhs_probe_length_adjustment = int(args[2])
        else:
            lhs_probe_length_adjustment = 0
            rhs_probe_length_adjustment = 0
        if self.config.flex_0bp_genotyping and is_gapfill:
            shift = int(args[1])
        carryover_key = (start_idx, lhs_probe_length_adjustment, rhs_probe_length_adjustment, gap_length)
        score = 0.0

        lhs_length = self.config.lhs_probe_length + lhs_probe_length_adjustment
        rhs_length = self.config.rhs_probe_length + rhs_probe_length_adjustment

        if carryover_key in carryover_visits:
            return carryover_visits[carryover_key]

        if is_gapfill and not self.config.flex_0bp_genotyping:
            rhs_probe = transcript_sequence[start_idx:start_idx+rhs_length] #### these were reversed! now fixed
            gap_sequence = transcript_sequence[start_idx+rhs_length:start_idx+rhs_length+gap_length]
            lhs_probe = transcript_sequence[start_idx+rhs_length+gap_length:start_idx+rhs_length+gap_length+lhs_length] #### these were reversed!  now fixed
            # Find where in the gap sequence the target position is
            gap_target_start_position = (target_position - start_idx - rhs_length)
            gap_target_end_position = gap_target_start_position + target_length

            original_target_start_gap=(original_target_start - 1) + ((start_idx + len(rhs_probe)) - (start_idx + gap_target_start_position + len(rhs_probe))) if original_target_start is not None else None
            original_target_end_gap=(start_idx+rhs_length+gap_length) + (original_target_end - (start_idx + rhs_length + gap_target_end_position)) if original_target_end is not None else None
            original_gap_sequence = original_transcript_sequence[original_target_start_gap:original_target_end_gap]

            # Check for any inconsistencies with the gap fill
            if len(gap_sequence) != gap_length:
                carryover_visits[carryover_key] = self.config.invalid_score
                return self.config.invalid_score
            if len(gap_sequence) < self.config.min_bridge_length or len(gap_sequence) > self.config.max_bridge_length:
                carryover_visits[carryover_key] = self.config.invalid_score
                return self.config.invalid_score
            if target_length > 0:
                if gap_target_start_position < 0 or gap_target_start_position >= len(gap_sequence) or gap_target_end_position < 0 or gap_target_end_position > len(gap_sequence):
                    carryover_visits[carryover_key] = self.config.invalid_score
                    return self.config.invalid_score
            else: ## in this case gap_target_start_position == gap_target_end_position
                if gap_target_start_position < 0 or gap_target_start_position > len(gap_sequence):
                    return self.config.invalid_score
            args = (expect_hits, lhs_probe, rhs_probe, gap_sequence, original_gap_sequence, gap_target_start_position, gap_target_end_position)
        else:
            rhs_probe = transcript_sequence[start_idx:start_idx+rhs_length]#### these were reversed! now fixed
            lhs_probe = transcript_sequence[start_idx+rhs_length:start_idx+lhs_length+rhs_length] #### these were reversed! now fixed

            args = (expect_hits, lhs_probe, rhs_probe)

        if self.config.flex_0bp_genotyping and is_gapfill:
            original_gap_sequence = ''
            rhs_probe = original_transcript_sequence[start_idx + shift:start_idx+rhs_length+shift]
            lhs_probe = original_transcript_sequence[start_idx+rhs_length + shift:start_idx+lhs_length+rhs_length+shift]
            if original_transcript_sequence == transcript_sequence:
                return self.config.invalid_score ### don't make dual probes for a region

        # Filter by predefined lhs if available
        if predefined_lhs is not None:
            if reverse_complement(lhs_probe) != predefined_lhs:
                carryover_visits[carryover_key] = self.config.invalid_score
                return self.config.invalid_score
           
        # Inconsistencies with the probe
        if len(lhs_probe) != lhs_length:
            carryover_visits[carryover_key] = self.config.invalid_score
            return self.config.invalid_score
        if len(rhs_probe) != rhs_length:
            carryover_visits[carryover_key] = self.config.invalid_score
            return self.config.invalid_score

        if is_gapfill:
            if self.config.flex_0bp_genotyping:
                full_original_contig = rhs_probe + lhs_probe
            else:
                full_original_contig = rhs_probe + original_gap_sequence + lhs_probe

            if full_original_contig not in original_transcript_sequence:
                print('ERROR finding original contig for:',transcript_name, rhs_probe, original_gap_sequence, lhs_probe)

        # Check for overlap with existing probes
        for name in existing_probes.keys():
            for probe in existing_probes[name]:
                if probe.original_transcript_sequence == original_transcript_sequence:
                    existing_start = probe.rhs_gene_start
                    existing_original_gap_length = probe.original_target_end_gap - probe.original_target_start_gap if probe.is_gapfill else 0
                    existing_end = probe.rhs_gene_start + len(probe.lhs_gene_sequence) + len(probe.rhs_gene_sequence) + existing_original_gap_length
                    new_start = start_idx
                    new_original_gap_length = len(original_gap_sequence) if is_gapfill else 0
                    new_end = start_idx + len(lhs_probe) + len(rhs_probe) + new_original_gap_length
                    # Calculate overlap amount
                    if not ((new_end < existing_start) or (new_start > existing_end)):
                        overlap_amount = min(new_end, existing_end) - max(new_start, existing_start)
                        if overlap_amount > self.config.max_probe_overlap:
                            carryover_visits[carryover_key] = self.config.invalid_score
                            return self.config.invalid_score

        # Skip if manually excluded
        if (reverse_complement(lhs_probe),reverse_complement(rhs_probe)) in self.config.exclude_probes:
            carryover_visits[carryover_key] = self.config.invalid_score
            return self.config.invalid_score

        flex_overlap = 0
        # Check for overlap with flex probes
        for i, probe in flex_probes.items():
            if probe.original_transcript_sequence == original_transcript_sequence:
                existing_start = probe.rhs_gene_start
                existing_original_gap_length = 0
                existing_end = probe.rhs_gene_start + len(probe.lhs_gene_sequence) + len(probe.rhs_gene_sequence) + existing_original_gap_length
                new_start = start_idx
                new_original_gap_length = len(original_gap_sequence) if is_gapfill else 0
                new_end = start_idx + len(lhs_probe) + len(rhs_probe) + new_original_gap_length
                if not ((new_end <= existing_start) or (new_start >= existing_end)):
                    if (new_start < existing_end) & (new_end >= existing_end):
                        flex_overlap += (existing_end - new_start)
                    elif (existing_start < new_end) & (existing_end >= new_end):
                        flex_overlap += (new_end - existing_start)
                    else:
                        raise ValueError("Unexpected overlap condition.")

        if self.config.flex_0bp_genotyping and is_gapfill:
            score_mutated = self.score_probe(*args, flex_overlap=flex_overlap)
            args = list(args)
            mut_lhs_probe = args[1]
            mut_rhs_probe = args[2]
            args[1] = lhs_probe
            args[2] = rhs_probe
            args = tuple(args)
            score_nonmutated = self.score_probe(*args, flex_overlap=flex_overlap)
            score = np.mean([score_mutated, score_nonmutated])
            variant_side = None
            ### check if probes overlap with each other if they both can bind to the unexpected (wt/mut) sequence; this can happen with in/dels; allow if only option
            if (mut_lhs_probe in original_transcript_sequence) & (mut_rhs_probe in original_transcript_sequence):
                lhs_index = original_transcript_sequence.find(mut_lhs_probe)
                rhs_index = original_transcript_sequence.find(mut_rhs_probe)
                if lhs_index != -1 and rhs_index != -1:
                    lhs_end = lhs_index + len(mut_lhs_probe)
                    rhs_end = rhs_index + len(mut_rhs_probe)
                    # If the probes' regions overlap in the original sequence, mark the score as invalid
                    if not (lhs_end <= rhs_index or rhs_end <= lhs_index):
                        score = self.config.invalid_score - 1
            if (lhs_probe in transcript_sequence) & (rhs_probe in transcript_sequence):
                lhs_index = transcript_sequence.find(lhs_probe)
                rhs_index = transcript_sequence.find(rhs_probe)
                if lhs_index != -1 and rhs_index != -1:
                    lhs_end = lhs_index + len(lhs_probe)
                    rhs_end = rhs_index + len(rhs_probe)
                    # If the probes' regions overlap in the mutated sequence, mark the score as invalid
                    if not (lhs_end <= rhs_index or rhs_end <= lhs_index):
                        score = self.config.invalid_score - 1
            if (mut_lhs_probe == lhs_probe) & (mut_rhs_probe != rhs_probe): 
                variant_side = 'rhs'
                for empirical_offset in range(0, len(mut_rhs_probe)):
                    if (mut_rhs_probe[:-empirical_offset:] in original_transcript_sequence) & (rhs_probe[:-empirical_offset:] in transcript_sequence):
                        score += empirical_offset * 1e7 ### penalize if not on the end; this should supercede other scoring penalties
                        found = True
                        break
                if not found:
                    score = self.config.invalid_score 
            elif (mut_rhs_probe == rhs_probe) & (mut_lhs_probe != lhs_probe):
                variant_side = 'lhs'
                for empirical_offset in range(0, len(mut_lhs_probe)):
                    if (mut_lhs_probe[empirical_offset:] in original_transcript_sequence) & (lhs_probe[empirical_offset:] in transcript_sequence):
                        score += empirical_offset * 1e7 ### penalize if not on the end; this should supercede other scoring penalties
                        found = True
                        break
                if not found:
                    score = self.config.invalid_score 
            if self.config.variant_side != None:
                if variant_side != self.config.variant_side:
                    score = self.config.invalid_score
        else:
            score += self.score_probe(*args, flex_overlap=flex_overlap)
        if carryover_key not in carryover_visits:
            carryover_visits[carryover_key] = score
        return score

    def _simulated_annealing_probe_search(self,
                                          transcript_name: str,
                                          transcript_sequence: str,
                                          expect_hits: int = 1,
                                          target_start: int = None,
                                          target_end: int = None,
                                          original_transcript_sequence: str = None,
                                          original_target_start: int = None,
                                          original_target_end: int = None,
                                          n_probes: int = 1,
                                          max_iterations: int = 500,
                                          initial_temp: float = 500.,
                                          fast: bool = False,
                                          existing_probes: dict[str, list[Probe]] = None,
                                          flex_probes: dict[str, list[Probe]] = None,
                                          predefined_lhs: str = None
                                          ) -> list[Probe]:
        """
        Create probes for a transcript using simulated annealing.
        :param transcript_name: The name of the transcript.
        :param transcript_sequence: The sequence to target.
        :param expect_hits: The number of expected hits in the database.
        :param target_start: The start of the target region for gap fill. If None, no gap fill will be performed.
        :param target_end: The end (exclusive) of the target region for gap fill.
        :param original_transcript_sequence: The original transcript sequence if the target region is modified. Only applicable for mutation calling.
        :param original_target_start: The original target start if the target region is modified. Only applicable for mutation calling.
        :param original_target_end: The original target end if the target region is modified. Only applicable for mutation calling.
        :param n_probes: The number of probes to output.
        :param max_iterations: The maximum number of iterations to perform within the simulated annealing.
        :param initial_temp: The initial temperature of the simulated annealing.
        :param fast: If True, we will use a fast approximate search for only "optimal" probes, rather than a full search.
        :return: The list of probes.
        """
        if fast and n_probes > 1:
            print("WARNING: Fast mode is only compatible with a single probe."
                  "Falling back to slow mode.")
            fast = False
        is_gapfill = target_start is not None
        target_length = None

        if self.config.flex_0bp_genotyping:
            if (target_end is None) or (target_start is None):
                is_gapfill = False
            else:
                target_length = target_end - target_start
                is_gapfill = True

        if is_gapfill and not self.config.flex_0bp_genotyping:
            if target_end is None:  # Implicit end
                target_end = target_start
            # Convert from 1-indexed, inclusive to 0-indexed, exclusive range
            target_start -= 1
            target_length = target_end - target_start
            min_bridge_length = max(self.config.min_bridge_length, target_length)
            max_bridge_length = self.config.max_bridge_length
            if min_bridge_length > max_bridge_length:
                return []
            if max_bridge_length < target_length:
                return []
            if target_start - self.config.lhs_probe_length < 0:
                return []
            if len(transcript_sequence) < target_end:
                return []
            if len(transcript_sequence) < self.config.lhs_probe_length + self.config.rhs_probe_length + min_bridge_length - self.config.flexible_probe_length_range:
                return []
        else:
            if len(transcript_sequence) < self.config.lhs_probe_length + self.config.rhs_probe_length - self.config.flexible_probe_length_range:
                return []  # The sequence is too short to create probes

        if fast:  # Approximate search
            probe = self._fast_probe_search(
                transcript_name=transcript_name,
                transcript_sequence=transcript_sequence,
                target_start=target_start,
                target_end=target_end,
                original_transcript_sequence=original_transcript_sequence,
                original_target_start=original_target_start - 1 if original_target_start else None,  # Convert to 0-indexed
                original_target_end=original_target_end,
            )
            if probe is None:
                return []
            return [probe]

        # We will perform a simulated annealing search to find the best probes
        # Probe bounds based on the LHS start index
        visited = dict()
        has_flexible_probe_length = self.config.flexible_probe_length_range > 0
        bounds = []
        if is_gapfill and not self.config.flex_0bp_genotyping:
            # Start position bounds
            bounds.append((
                max(0, target_start - self.config.lhs_probe_length - self.config.max_bridge_length - self.config.flexible_probe_length_range + target_length),
                min(target_start - self.config.lhs_probe_length, len(transcript_sequence)-1)
            ))
            if has_flexible_probe_length:
                # LHS probe length adjustment bounds
                bounds.append((-self.config.flexible_probe_length_range, self.config.flexible_probe_length_range))
                # RHS probe length adjustment bounds
                bounds.append((-self.config.flexible_probe_length_range, self.config.flexible_probe_length_range))
            bounds.append((
                # Bounds for the gap length
                min(max(self.config.min_bridge_length, target_length), self.config.max_bridge_length),
                self.config.max_bridge_length
            ))
        elif is_gapfill and self.config.flex_0bp_genotyping:
            length_diff = np.abs(len(transcript_sequence) - len(original_transcript_sequence)) + 1
            ## fix bounds for 0bp genotyping when variant side is set because otherwise desired side may not be visited
            # Start position bounds
            bounds.append((
                max(0, target_start - self.config.lhs_probe_length - self.config.max_bridge_length - self.config.flexible_probe_length_range - target_length - length_diff),
                min(target_start - self.config.lhs_probe_length + target_length + length_diff, len(transcript_sequence)-1)
            ))
            if has_flexible_probe_length:
                # LHS probe length adjustment bounds
                bounds.append((-self.config.flexible_probe_length_range, self.config.flexible_probe_length_range))
                # RHS probe length adjustment bounds
                bounds.append((-self.config.flexible_probe_length_range, self.config.flexible_probe_length_range))
            bounds.append((
                # Bounds for the probe shift
                -target_length - 1,
                target_length + 1
            ))
        else:
            bounds.append((
                0,
                len(transcript_sequence) - self.config.lhs_probe_length - self.config.rhs_probe_length + self.config.flexible_probe_length_range
            ))
            if has_flexible_probe_length:
                # LHS probe length adjustment bounds
                bounds.append((-self.config.flexible_probe_length_range, self.config.flexible_probe_length_range))
                # RHS probe length adjustment bounds
                bounds.append((-self.config.flexible_probe_length_range, self.config.flexible_probe_length_range))
        # Compute the number of possibilities if we were to search the entire space
        possible_variations = np.prod([b[1] - b[0] + 1 for b in bounds])
        # If the number of variations is low, just do a brute force search
        if possible_variations <= max_iterations:
            brute_force = True
        else:
            brute_force = False

        probes = []

        for _ in range(n_probes):
            try:
                if brute_force:
                    # Brute force will be more efficient with slice objects specifying the exact bounds
                    bounds = [slice(b[0], b[1]+1, 1) for b in bounds]
                    result = brute(
                        func=self._test_flex_probe,
                        args=(transcript_name, transcript_sequence, original_transcript_sequence, expect_hits, is_gapfill, target_start if is_gapfill else None, target_length if is_gapfill else None, visited, existing_probes, flex_probes, original_target_start if is_gapfill else None, original_target_end if is_gapfill else None, predefined_lhs),
                        ranges=bounds,
                        full_output=True
                    )
                    result_args = [int(x) for x in np.round(result[0]).astype(int)]
                    score = result[1]
                else:
                    x0 = [np.random.randint(b[0], b[1]+1) for b in bounds]
                    print(max_iterations)
                    result = dual_annealing(
                        func=self._test_flex_probe,
                        bounds=bounds,
                        x0=np.array(x0),
                        args=(transcript_name, transcript_sequence, original_transcript_sequence, expect_hits, is_gapfill, target_start if is_gapfill else None, target_length if is_gapfill else None, visited, existing_probes, flex_probes, original_target_start if is_gapfill else None, original_target_end if is_gapfill else None, predefined_lhs),
                        maxiter=max_iterations,
                        initial_temp=initial_temp
                    )
                    print(result)
                    result_args = [int(x) for x in np.round(result.x).astype(int)]
                    score = result.fun
            except ValueError as e:
                if "negative dimensions" in str(e):
                    continue
                else:
                    raise

            # Note that we flip the RHS and LHS to follow 10x nomenclature here
            if is_gapfill:
                start_idx = result_args[0]
                if has_flexible_probe_length:
                    rhs_probe_length_adjustment = result_args[1]
                    lhs_probe_length_adjustment = result_args[2]
                else:
                    rhs_probe_length_adjustment = 0
                    lhs_probe_length_adjustment = 0
                gap_length = result_args[-1]
            else:
                start_idx = result_args[0]
                if has_flexible_probe_length:
                    rhs_probe_length_adjustment = result_args[1]
                    lhs_probe_length_adjustment = result_args[2]
                else:
                    rhs_probe_length_adjustment = 0
                    lhs_probe_length_adjustment = 0
                gap_length = 0
            if self.config.flex_0bp_genotyping and is_gapfill:
                shift = result_args[1]
                gap_length = 0
            # Again we are careful to flip the RHS and LHS to follow 10x nomenclature
            rhs_length = self.config.rhs_probe_length + rhs_probe_length_adjustment
            lhs_length = self.config.lhs_probe_length + lhs_probe_length_adjustment
            # We will extract the probes in 10X orientation (i.e. LHS probe is on the right of the gene and vice versa)
            rhs_gene_seq = transcript_sequence[start_idx:start_idx+rhs_length]
            lhs_gene_seq = transcript_sequence[start_idx+rhs_length+gap_length:start_idx+rhs_length+gap_length+lhs_length]
            gap_start = start_idx + rhs_length
            gap_end = start_idx+rhs_length+gap_length
            gapfill_gene_seq = transcript_sequence[start_idx+rhs_length:start_idx+rhs_length+gap_length]

            rhs_probe_seq = reverse_complement(rhs_gene_seq)
            lhs_probe_seq = reverse_complement(lhs_gene_seq)
            gapfill_probe_seq = reverse_complement(gapfill_gene_seq)

            if self.config.flex_0bp_genotyping and is_gapfill:
                nonmutated_rhs_gene_seq = original_transcript_sequence[start_idx+shift:start_idx+rhs_length+shift]
                nonmutated_lhs_gene_seq = original_transcript_sequence[start_idx+rhs_length+shift:start_idx+lhs_length+rhs_length+shift]
                nonmutated_rhs_probe_seq = reverse_complement(nonmutated_rhs_gene_seq)
                nonmutated_lhs_probe_seq = reverse_complement(nonmutated_lhs_gene_seq)

            if score < self.config.invalid_score:
                # rhs_offtarget_hits, rhs_offtarget_names = self.blast_hits(rhs_gene_seq)
                # lhs_offtarget_hits, lhs_offtarget_names = self.blast_hits(lhs_gene_seq)
                original_target_start=original_target_start - 1 if original_target_start else None  # Convert to 0-indexed
                probes.append(Probe(
                    transcript_name=transcript_name,
                    rhs_probe=rhs_probe_seq,
                    lhs_probe=lhs_probe_seq,
                    rhs_gene_start=start_idx,
                    lhs_gene_start=start_idx+rhs_length+gap_length,
                    rhs_gene_sequence=rhs_gene_seq,
                    lhs_gene_sequence=lhs_gene_seq,
                    rhs_GC=(rhs_probe_seq.count("G") + rhs_probe_seq.count("C")) / len(rhs_probe_seq) if len(rhs_probe_seq) > 0 else 0.5,
                    lhs_GC=(lhs_probe_seq.count("G") + lhs_probe_seq.count("C")) / len(lhs_probe_seq) if len(lhs_probe_seq) > 0 else 0.5,
                    expected_blast_hits=expect_hits,
                    total_rhs_blast_hits=[],
                    total_lhs_blast_hits=[],
                    rhs_blast_hits=[],
                    lhs_blast_hits=[],
                    score=score,
                    transcript_sequence=transcript_sequence,
                    gap_length=gap_length if is_gapfill else None,
                    gap_probe_sequence=gapfill_probe_seq if is_gapfill else None,
                    gap_gene_sequence=gapfill_gene_seq if is_gapfill else None,
                    target_start_gap=gap_start,
                    target_end_gap=gap_end,
                    original_transcript_sequence=original_transcript_sequence,
                    # Add the offset we have to the original target start and end
                    original_target_start_gap=original_target_start + (gap_start - target_start) if original_target_start is not None else None,
                    original_target_end_gap=original_target_end + (gap_end - target_end) if original_target_end is not None else None,
                    nonmutated_rhs_probe = nonmutated_rhs_probe_seq if self.config.flex_0bp_genotyping and is_gapfill else None,
                    nonmutated_lhs_probe = nonmutated_lhs_probe_seq if self.config.flex_0bp_genotyping and is_gapfill else None
                ))
        return probes

    def _fast_probe_test(self, lhs_seq: str, rhs_seq: str, gap_seq: str = None) -> bool:
        """
        Basic filter to determine if a probe has potential.
        :param lhs_seq: The left-hand side sequence (gene space).
        :param rhs_seq: The right-hand side sequence (gene space).
        :param gap_seq: The gap sequence (gene space).
        :return: True if passes, False otherwise.
        """
        if len(lhs_seq) != self.config.lhs_probe_length or len(rhs_seq) != self.config.rhs_probe_length:
            return False

        lhs_probe_seq = reverse_complement(lhs_seq)
        rhs_probe_seq = reverse_complement(rhs_seq)
        # Test junction
        if gap_seq is not None:
            if len(gap_seq) < self.config.min_bridge_length or len(gap_seq) > self.config.max_bridge_length:
                return False

            gap_probe_seq = reverse_complement(gap_seq)
            junction = gap_probe_seq[-1].upper() + rhs_probe_seq[0].upper()
        else:
            junction = lhs_probe_seq[-1].upper() + rhs_probe_seq[0].upper()

        # Empirically "good" junctions
        if junction not in ["TA", "TC", "TG", "TT", "CT", "CA", "AT"]:
            return False

        rhs_gc = (rhs_probe_seq.count("G") + rhs_probe_seq.count("C")) / len(rhs_probe_seq)
        lhs_gc = (lhs_probe_seq.count("G") + lhs_probe_seq.count("C")) / len(lhs_probe_seq)

        if lhs_gc < self.config.tx_min_gc or lhs_gc > self.config.tx_max_gc:
            return False
        if rhs_gc < self.config.tx_min_gc or rhs_gc > self.config.tx_max_gc:
            return False

        return True

    def _fast_probe_search(self,
                           transcript_name: str,
                           transcript_sequence: str,
                           target_start: int = None,
                           target_end: int = None,
                           original_transcript_sequence: str = None,
                           original_target_start: int = None,
                           original_target_end: int = None) -> Optional[Probe]:
        """
        Only creates a single that are optimal for the target sequence without a BLAST search.
        For gapfilling probes, we search from the target site, outwards. With no gapfilling, we search from the start
            of the sequence until something is found.
        :param transcript_name: The name of the sequence.
        :param transcript_sequence: The sequence to target.
        :param target_start: The start of the target region for gap fill.
        :param target_end: The end (exclusive) of the target region for gap fill. If None, assume a SNV.
        :param original_transcript_sequence: The original sequence being targetted. Only applicable for mutation calling.
        :param original_target_start: The original start of the target region for gap fill. Only applicable for mutation calling.
        :param original_target_end: The original end of the target region for gap fill. Only applicable for mutation calling.
        :return: A probe if found, else None.
        """
        if target_start is not None:  # Gapfilling
            target_length = target_end - target_start
            min_bridge = max(self.config.min_bridge_length, target_length)
            max_bridge = self.config.max_bridge_length

            # Search from the target site outwards
            for gap_length in range(min_bridge, max_bridge+1):
                # Compute all the windows to test
                # Do not allow the target to be at the start of the gap sequence since it can effect junctions.
                # Case 1: Gap length == target length, only one probe to test  <- Won't be allowed to happen since the start of the target is not allowed to be at the start of the gap sequence
                # Case 2: Gap length > target length, test all windows containing the target
                # Case 3: Gap length < target length, test all windows containing the target  <- Shouldn't be possible

                free_bases = gap_length - target_length
                for i in range(-free_bases, 0):  # Note that we exclude the last offset to prevent the target from starting in the junction
                    start = target_start + i
                    end = start + gap_length
                    if start < 0 or end > len(transcript_sequence):
                        continue
                    lhs_probe = transcript_sequence[start - self.config.lhs_probe_length:start]
                    rhs_probe = transcript_sequence[end:end + self.config.rhs_probe_length]
                    gap_probe = transcript_sequence[start:end]
                    if self._fast_probe_test(lhs_probe, rhs_probe, gap_probe):
                        return Probe(
                            transcript_name=transcript_name,
                            rhs_probe=reverse_complement(lhs_probe),
                            lhs_probe=reverse_complement(rhs_probe),
                            rhs_gene_start=start - self.config.lhs_probe_length,
                            lhs_gene_start=end,
                            rhs_gene_sequence=lhs_probe,
                            lhs_gene_sequence=rhs_probe,
                            rhs_GC=(lhs_probe.count("G") + lhs_probe.count("C")) / len(lhs_probe),
                            lhs_GC=(rhs_probe.count("G") + rhs_probe.count("C")) / len(rhs_probe),
                            expected_blast_hits=None,  # BLAST info will be filled in later with a bulk BLAST search
                            total_rhs_blast_hits=None,
                            total_lhs_blast_hits=None,
                            rhs_blast_hits=None,
                            lhs_blast_hits=None,
                            score=1.0,
                            transcript_sequence=transcript_sequence,
                            gap_length=gap_length,
                            gap_probe_sequence=reverse_complement(gap_probe),
                            gap_gene_sequence=gap_probe,
                            target_start_gap=start,
                            target_end_gap=end,
                            original_transcript_sequence=original_transcript_sequence,
                            # Add the offset we are using for our probe
                            original_target_start_gap=original_target_start + i if original_target_start else None,
                            # Add the distance between the target region and the end of the gap to the original target end
                            original_target_end_gap=original_target_end + i + (end - target_start) if original_target_end else None,
                        )
                    else:
                        continue
        else:  # Non-gapfilling
            for start in range(len(transcript_sequence) - self.config.lhs_probe_length - self.config.rhs_probe_length + 1):
                lhs_probe = transcript_sequence[start:start+self.config.lhs_probe_length]
                rhs_probe = transcript_sequence[start+self.config.lhs_probe_length:start+self.config.lhs_probe_length+self.config.rhs_probe_length]
                if self._fast_probe_test(lhs_probe, rhs_probe):
                    return Probe(
                        transcript_name=transcript_name,
                        rhs_probe=reverse_complement(lhs_probe),
                        lhs_probe=reverse_complement(rhs_probe),
                        rhs_gene_start=start,
                        lhs_gene_start=start+self.config.lhs_probe_length,
                        rhs_gene_sequence=lhs_probe,
                        lhs_gene_sequence=rhs_probe,
                        rhs_GC=(lhs_probe.count("G") + lhs_probe.count("C")) / len(lhs_probe),
                        lhs_GC=(rhs_probe.count("G") + rhs_probe.count("C")) / len(rhs_probe),
                        expected_blast_hits=None,  # BLAST info will be filled in later with a bulk BLAST search
                        total_rhs_blast_hits=None,
                        total_lhs_blast_hits=None,
                        score=1.0,
                    )
                else:
                    continue
        return None

    def bulk_blast_filter_probes(self,
                                 lhs_seqs: list[str],
                                 rhs_seqs: list[str],
                                 expect_hits: list[int] = None) -> tuple[list[int], list[tuple[int, int]], list[tuple[list[str], list[str]]]]:
        """
        Perform a bulk BLAST search to filter out probes.
        :param lhs_seqs: The left-hand side probe sequences.
        :param rhs_seqs: The right-hand side probe sequences.
        :param expect_hits: The expected number of hits to the transcriptome.
        :return: The indices of probes to remove, and the list of all LHS and RHS hits as well as a list of the gene names.
        """
        to_remove = []
        all_hits = []
        all_names = []

        assert len(lhs_seqs) == len(rhs_seqs)

        print(f"Performing a bulk BLAST search on {len(lhs_seqs)} probes.")

        lhs_results = blast_search(self.blast_db, lhs_seqs, evalue=self.config.evalue_cutoff)
        rhs_results = blast_search(self.blast_db, rhs_seqs, evalue=self.config.evalue_cutoff)

        assert len(lhs_results) == len(lhs_seqs) and len(rhs_results) == len(rhs_seqs)

        for i, (lhs_result, rhs_result, expect_hits) in enumerate(zip(lhs_results, rhs_results, expect_hits or [1]*len(lhs_seqs))):
            lhs_hits = defaultdict(list)
            rhs_hits = defaultdict(list)
            for aln in lhs_result.alignments:
                name = aln.hit_def.split(" ")[0]
                # 10X recommends at least 5 mismatches for a probe to not hybridize efficiently, so select alignments with < 5 mismatches to be "hits"
                matches = sum([hsp.identities for hsp in aln.hsps])
                if matches >= len(lhs_seqs[i]) - 5:
                    lhs_hits[name].append(aln)
            for aln in rhs_result.alignments:
                name = aln.hit_def.split(" ")[0]
                matches = sum([hsp.identities for hsp in aln.hsps])
                if matches >= len(rhs_seqs[i]) - 5:
                    rhs_hits[name].append(aln)

            all_hits.append((len(lhs_hits), len(rhs_hits)))
            all_names.append((list(lhs_hits.keys()), list(rhs_hits.keys())))
            lhs_hits = len(lhs_hits) - expect_hits
            rhs_hits = len(rhs_hits) - expect_hits

            if lhs_hits > 0 or rhs_hits > 0:
                to_remove.append(i)

        return to_remove, all_hits, all_names

    def create_flex_probes(self,
                           transcript_name: str,
                           transcript_sequence: str,
                           expect_hits: int = 1,
                           n_probes: int = 1,
                           max_iterations: int = 1_000,
                           initial_temp: float = 500.,
                           fast: bool = False,
                           existing_probes: dict[str, list[Probe]] = None,
                           flex_probes: dict[str, list[Probe]] = None
                           ) -> list[Probe]:
        """
        Create vanilla non-gap-filling probes for a transcript. Because of the large search space, we will perform
            dual simulated annealing to find the best probes.
        :param transcript_name: The name of the sequence.
        :param transcript_sequence: The sequence to target.
        :param expect_hits: The expected number of hits to the transcriptome.
        :param n_probes: The number of probes to output.
        :param max_iterations: The maximum number of iterations to perform within the simulated annealing.
        :param initial_temp: The initial temperature of the simulated annealing.
        :param fast: If True, we will use a fast approximate search for only "optimal" probes, rather than a full search. Also skips BLAST searches.
        :return: The list of probes.
        """
        return self._simulated_annealing_probe_search(
            transcript_name=transcript_name,
            transcript_sequence=transcript_sequence,
            expect_hits=expect_hits,
            n_probes=n_probes,
            max_iterations=max_iterations,
            initial_temp=initial_temp,
            fast=fast,
            existing_probes=existing_probes,
            flex_probes=flex_probes
        )

    def create_gapfilling_flex_probes(self,
                                      transcript_name: str,
                                      transcript_sequence: str,
                                      target_start: int,
                                      target_end: int = None,
                                      original_transcript_sequence: Optional[str] = None,
                                      original_target_start: Optional[int] = None,
                                      original_target_end: Optional[int] = None,
                                      expect_hits: int = 1,
                                      n_probes: int = 1,
                                      max_iterations: int = 1_000,
                                      initial_temp: float = 500.,
                                      fast: bool = False,
                                      existing_probes: dict[str, list[Probe]] = None,
                                      flex_probes: dict[str, list[Probe]] = None,
                                      predefined_lhs: str = None
                                      ) -> list[Probe]:
        """
        Create gap-filling probes for a transcript. Because of the large search space, we will perform
            dual simulated annealing to find the best probes.
        :param transcript_name: The name of the sequence.
        :param transcript_sequence: The sequence to target.
        :param target_start: The start of the target region for gap fill. 1-indexed.
        :param target_end: The end (exclusive) of the target region for gap fill. If None, assume a SNV.
        :param original_transcript_sequence: The original sequence being targetted. Only applicable for mutation calling.
        :param original_target_start: The original start of the target region for gap fill. Only applicable for mutation calling.
        :param original_target_end: The original end of the target region for gap fill. Only applicable for mutation calling.
        :param expect_hits: The expected number of hits to the transcriptome.
        :param n_probes: The number of probes to output.
        :param max_iterations: The maximum number of iterations to perform within the simulated annealing.
        :param initial_temp: The initial temperature of the simulated annealing.
        :param fast: If True, we will use a fast approximate search for only "optimal" probes, rather than a full search. Also skips BLAST searches.
        :return: The list of probes.
        """
        return self._simulated_annealing_probe_search(
            transcript_name=transcript_name,
            transcript_sequence=transcript_sequence,
            expect_hits=expect_hits,
            target_start=target_start,
            target_end=target_end,
            original_transcript_sequence=original_transcript_sequence,
            original_target_start=original_target_start,
            original_target_end=original_target_end,
            n_probes=n_probes,
            max_iterations=max_iterations,
            initial_temp=initial_temp,
            fast=fast,
            existing_probes=existing_probes,
            flex_probes=flex_probes,
            predefined_lhs=predefined_lhs
        )

    def convert_probe_set_to_df(self, probes: dict[str, list[Probe]], visium: bool = False, truseq: bool = False, barcode: int | list[int] = 0) -> pd.DataFrame:
        """
        Convert a set of probes to a dataframe.
        :param probes: The dictionary of probe identifier -> list of designed probes.
        :param visium: Whether to design probes for Visium.
        :param truseq: Whether to replace the standard handle with the TruSeq handle.
        :param barcode: The barcode(s) to use for the probe.
        :return: The dataframe.
        """
        df = dict(
            name=[],  # The name of the transcript

            lhs_probe=[],  # The left-hand side probe sequence
            rhs_probe=[],  # The right-hand side probe sequence
            lhs_gene_sequence=[],  # The left-hand side gene sequence
            rhs_gene_sequence=[],  # The right-hand side gene sequence

            lhs_IDT_order=[],  # The IDT order for the left-hand side probe
            rhs_IDT_order=[],  # The IDT order for the right-hand side probe

            lhs_GC=[],  # The GC content of the left-hand side probe
            rhs_GC=[],  # The GC content of the right-hand side probe
            GC_difference=[],  # The difference in GC content between the left-hand and right-hand side probes

            gap_length=[],  # The length of the gap
            gap_probe_sequence=[],  # The sequence of the gap
            gap_gene_sequence=[],  # The gene sequence of the gap
            target_gap_start=[],  # The start of the gap in the target sequence
            target_gap_end=[],  # The end of the gap in the target sequence

            expected_blast_hits=[],  # The expected number of hits to the transcriptome for the LHS or RHS
            total_lhs_blast_hits=[],  # The total number of hits to the transcriptome for the LHS
            total_rhs_blast_hits=[],  # The total number of hits to the transcriptome for the RHS
            rhs_blast_hits=[],  # The names of the off-target hits
            lhs_blast_hits=[],  # The names of the off-target hits

            score=[],  # The score of the probe

            barcode=[],  # The barcode used for the probe
            transcript=[],  # The full transcript sequence we targetted again
            # For gapfilling SNV probes:
            original_transcript_sequence=[],  # The original transcript sequence if the target region is modified. Only applicable for mutation calling.
            original_target_gap_start=[],  # The original start of the target region for gap fill. Only applicable for mutation calling.
            original_target_gap_end=[],  # The original end of the target region for gap fill. Only applicable for mutation calling.
            original_gap_probe_sequence=[],  # The original gap probe sequence if the target region is modified. Only applicable for mutation calling.
            original_gap_gene_sequence=[],  # The original gap gene sequence if the target region is modified. Only applicable for mutation calling.
            nonmutated_lhs_probe=[],  # The non-mutated left-hand side probe sequence. Only applicable for 0bp mutation calling.
            nonmutated_rhs_probe=[]  # The non-mutated right-hand side probe sequence. Only applicable for 0bp mutation calling.
        )
        if isinstance(barcode, int):
            barcodes = [barcode]
        else:
            barcodes = barcode
        for barcode in barcodes:
            has_gapfill = False
            for transcript_name, probe_list in probes.items():
                for probe in probe_list:
                    has_gapfill = has_gapfill or probe.is_gapfill

                    df["name"].append(transcript_name)
                    df["lhs_probe"].append(probe.lhs_probe)
                    df["lhs_gene_sequence"].append(probe.lhs_gene_sequence)
                    df["rhs_probe"].append(probe.rhs_probe)
                    df["rhs_gene_sequence"].append(probe.rhs_gene_sequence)
                    df["lhs_IDT_order"].append(probe.lhs_probe_order(truseq))
                    df["rhs_IDT_order"].append(probe.rhs_probe_order(barcode, visium))
                    df["lhs_GC"].append(probe.lhs_GC)
                    df["rhs_GC"].append(probe.rhs_GC)
                    df["GC_difference"].append(abs(probe.lhs_GC - probe.rhs_GC))
                    df["expected_blast_hits"].append(probe.expected_blast_hits)
                    df["total_lhs_blast_hits"].append(probe.total_lhs_blast_hits)
                    df["total_rhs_blast_hits"].append(probe.total_rhs_blast_hits)
                    df["rhs_blast_hits"].append(";".join(probe.rhs_blast_hits))
                    df["lhs_blast_hits"].append(";".join(probe.lhs_blast_hits))
                    df["score"].append(probe.score)
                    df["barcode"].append(barcode+1)
                    df["transcript"].append(probe.transcript_sequence)

                    if not probe.is_gapfill:
                        df["gap_length"].append(0)
                        df["gap_probe_sequence"].append("N/A")
                        df["gap_gene_sequence"].append("N/A")
                        df["target_gap_start"].append("N/A")
                        df["target_gap_end"].append("N/A")
                        df["original_transcript_sequence"].append("N/A")
                        df["original_target_gap_start"].append("N/A")
                        df["original_target_gap_end"].append("N/A")
                        df["original_gap_probe_sequence"].append("N/A")
                        df["original_gap_gene_sequence"].append("N/A")
                    else:
                        df["gap_length"].append(probe.gap_length)
                        df["gap_probe_sequence"].append(probe.gap_probe_sequence)
                        df["gap_gene_sequence"].append(probe.gap_gene_sequence)
                        df["target_gap_start"].append(probe.target_start_gap)
                        df["target_gap_end"].append(probe.target_end_gap)
                        has_original = probe.original_transcript_sequence is not None
                        df["original_transcript_sequence"].append(probe.original_transcript_sequence if has_original else "N/A")
                        df["original_target_gap_start"].append(probe.original_target_start_gap if has_original else "N/A")
                        df["original_target_gap_end"].append(probe.original_target_end_gap if has_original else "N/A")
                        df["original_gap_gene_sequence"].append(probe.original_transcript_sequence[probe.original_target_start_gap:probe.original_target_end_gap] if has_original else "N/A")
                        df["original_gap_probe_sequence"].append(reverse_complement(probe.original_transcript_sequence[probe.original_target_start_gap:probe.original_target_end_gap]) if has_original else "N/A")
                    if self.config.flex_0bp_genotyping and probe.is_gapfill:
                        df["nonmutated_lhs_probe"].append(probe.nonmutated_lhs_probe)
                        df["nonmutated_rhs_probe"].append(probe.nonmutated_rhs_probe)
                    else:
                        df["nonmutated_lhs_probe"].append("N/A")
                        df["nonmutated_rhs_probe"].append("N/A")


        if not has_gapfill:
            del df["gap_length"]
            del df["gap_probe_sequence"]
            del df["gap_gene_sequence"]
            del df['target_gap_start']
            del df['target_gap_end']
            del df["original_transcript_sequence"]
            del df["original_target_gap_start"]
            del df["original_target_gap_end"]
            del df["original_gap_probe_sequence"]
            del df["original_gap_gene_sequence"]

        return pd.DataFrame(df)

    def generate_flex_probe_set_df(self,
                                   transcript_sequences: dict[str, str],
                                   expect_hits: list[int] = None,
                                   n_probes: int = 1,
                                   max_iterations_per_probe: int = 1_000,
                                   initial_temp: float = 500.,
                                   visium: bool = False,
                                   truseq: bool = False,
                                   fast: bool = False,
                                   barcode: int | list[int] = 0) -> pd.DataFrame:
        """
        Create a set of probes for a given set of sequences.
        :param transcript_sequences: The sequences to target.
        :param expect_hits: The number of expected hits for each target. 1 if not specified.
        :param n_probes: The number of probes per target.
        :param max_iterations_per_probe: The max iterations per probe for simulated annealing.
        :param initial_temp: The initial temperature of each simulated annealing run.
        :param visium: Whether to design probes for Visium.
        :param truseq: Whether to replace the standard handle with the TruSeq handle.
        :param fast: If True, we will use a fast approximate search for only "optimal" probes, rather than a full search.
        :param barcode: The barcode(s) to use for the probe.
        :return: The dataframe result. Or if visium is True, a tuple of the flex probe dataframe and the visium probe dataframe.
        """
        all_probes = dict()

        if expect_hits is None:
            expect_hits = self.generate_expect_hits(list(transcript_sequences.values()))
        elif len(expect_hits) == 0:
            expect_hits = [1] * len(transcript_sequences)

        flex_probes = dict()
        flex_list = self.reference_probes
        checked_sequences = []
        i = 0
        for (name, sequence) in transcript_sequences.items():
            sequence = sequence[0]
            name = name.split(' ')[0]
            if sequence not in checked_sequences:
                checked_sequences.append(sequence)
                for probe_seq in flex_list['probe_seq']:
                    probe_gene_seq = reverse_complement(probe_seq)
                    if probe_gene_seq in sequence:
                        i += 1
                        flex_probes[i] = Probe(
                            transcript_name=name,
                            rhs_probe = probe_seq[25:],
                            lhs_probe = probe_seq[:25],
                            rhs_gene_start = sequence.index(probe_gene_seq),
                            lhs_gene_start = sequence.index(probe_gene_seq) + len(probe_seq[25:]),
                            rhs_gene_sequence = reverse_complement(probe_seq[25:]),
                            lhs_gene_sequence = reverse_complement(probe_seq[:25]),
                            rhs_GC = (probe_seq[25:].count("G") + probe_seq[25:].count("C")) / len(probe_seq[25:]) if len(probe_seq[25:]) > 0 else 0.5,
                            lhs_GC = (probe_seq[:25].count("G") + probe_seq[:25].count("C")) / len(probe_seq[:25]) if len(probe_seq[:25]) > 0 else 0.5,
                            expected_blast_hits = 1,
                            total_rhs_blast_hits = [],
                            total_lhs_blast_hits = [],
                            rhs_blast_hits = [],
                            lhs_blast_hits = [],
                            score = 0,
                            transcript_sequence = sequence,
                            gap_length = None,
                            gap_probe_sequence = None,
                            gap_gene_sequence = None,
                            target_start_gap = None,
                            target_end_gap = None,
                            original_transcript_sequence = sequence,
                            original_target_start_gap = None,
                            original_target_end_gap = None,
                        )

        for i, (transcript_name, transcript_sequence) in tqdm(enumerate(transcript_sequences.items()),
                                                              desc='Generating probes',
                                                              total=len(transcript_sequences),
                                                              unit='target',
                                                              dynamic_ncols=True):
            expected = expect_hits[i]
            print(transcript_name)
            print(transcript_sequence)
            print(len(transcript_sequence))
            probes = self.create_flex_probes(
                transcript_name,
                transcript_sequence,
                expected,
                n_probes,
                max_iterations_per_probe,
                initial_temp,
                existing_probes=all_probes,
                flex_probes=flex_probes,
                fast=fast
            )
            if probes is None or len(probes) == 0:
                print(f"Warning: No probe for {transcript_name}")
                all_probes[transcript_name] = []
            else:
                all_probes[transcript_name] = probes
        # If we are running with approximate search, we should filter out probes that have mismatches
        if fast:
            lhs_seqs = [probe.lhs_gene_sequence for probes in all_probes.values() for probe in probes]
            rhs_seqs = [probe.rhs_gene_sequence for probes in all_probes.values() for probe in probes]
            expect_hits = None if not expect_hits else [hits for hits, probes in zip(expect_hits, all_probes.values()) for _ in probes]
            to_remove, all_hits, all_names = self.bulk_blast_filter_probes(lhs_seqs, rhs_seqs, expect_hits)
            print(f"Removing {len(to_remove)} probes due to off-target activity.")
            curr_idx = 0
            for transcript_name, probes in list(all_probes.items()):
                removed = []
                for probe in probes:
                    if curr_idx in to_remove:
                        removed.append(probe)
                    else:
                        probe.total_lhs_blast_hits = all_hits[curr_idx][0]
                        probe.total_rhs_blast_hits = all_hits[curr_idx][1]
                        probe.expected_blast_hits = expect_hits[curr_idx]
                        probe.lhs_blast_hits = all_names[curr_idx][0]
                        probe.rhs_blast_hits = all_names[curr_idx][1]
                    curr_idx += 1
                all_probes[transcript_name] = [probe for probe in probes if probe not in removed]

        if visium:
            return self.convert_probe_set_to_df(all_probes, truseq=truseq, barcode=barcode), self.convert_probe_set_to_df(all_probes, truseq=truseq, barcode=barcode, visium=True)
        return self.convert_probe_set_to_df(all_probes, truseq=truseq, barcode=barcode)

    def generate_gapfilling_flex_probe_set_df(self,
                                                transcript_sequences: dict[str, str | tuple[str, str]],
                                                target_starts: list[int | tuple[int, int]],
                                                target_ends: list[int | tuple[int, int]],
                                                expect_hits: list[int] = None,
                                                n_probes: int = 1,
                                                max_iterations_per_probe: int = 100_000,
                                                initial_temp: float = 500.,
                                                visium: bool = False,
                                                truseq: bool = True,
                                                fast: bool = False,
                                                barcode: int | list[int] = 0,
                                                lhs_probes: list[str] = None) -> pd.DataFrame:
        """
        Create a set of gap-filling probes for a given set of sequences.
        :param transcript_sequences: The sequences to target. If a tuple, first element is the original transcript, the second is the mutated transcript.
        :param target_starts: The start of the target regions for gap fill. If a tuple, the first element is the start of the original target, the second is the start of the mutated target. 1-indexed.
        :param target_ends: The end (exclusive) of the target regions for gap fill. If a tuple, the first element is the end of the original target, the second is the end of the mutated target.
        :param expect_hits: The expected hits for each target to the transcriptome. 1 if not specified.
        :param n_probes: The number of probes per target.
        :param max_iterations_per_probe: The max iterations per probe for simulated annealing.
        :param initial_temp: The initial temperature of each simulated annealing run.
        :param visium: Include probes for Visium.
        :param truseq: Replace the standard handle with the TruSeq handle.
        :param fast: If True, we will use a fast approximate search for only "optimal" probes, rather than a full search.
        :param barcode: The barcode(s) to use for the probe.
        :param lhs_probes: used to pre-define desired LHS probes.
        :return: The dataframe result. Or if visium is True, a tuple of the flex probe dataframe and the visium probe dataframe.
        """
        all_probes = dict()

        missing = []
        if expect_hits is None:
            expect_hits = self.generate_expect_hits(list(transcript_sequences.values()))
        elif len(expect_hits) == 0:
            expect_hits = [1] * len(transcript_sequences)

        flex_probes = dict()
        flex_list = self.reference_probes
        checked_sequences = []
        i = 0
        for (name, sequence) in transcript_sequences.items():
            sequence = sequence[0]
            name = name.split(' ')[0]
            if sequence not in checked_sequences:
                checked_sequences.append(sequence)
                for probe_seq in flex_list['probe_seq']:
                    probe_gene_seq = reverse_complement(probe_seq)
                    if probe_gene_seq in sequence:
                        i += 1
                        flex_probes[i] = Probe(
                            transcript_name=name,
                            rhs_probe = probe_seq[25:],
                            lhs_probe = probe_seq[:25],
                            rhs_gene_start = sequence.index(probe_gene_seq),
                            lhs_gene_start = sequence.index(probe_gene_seq) + len(probe_seq[25:]),
                            rhs_gene_sequence = reverse_complement(probe_seq[25:]),
                            lhs_gene_sequence = reverse_complement(probe_seq[:25]),
                            rhs_GC = (probe_seq[25:].count("G") + probe_seq[25:].count("C")) / len(probe_seq[25:]) if len(probe_seq[25:]) > 0 else 0.5,
                            lhs_GC = (probe_seq[:25].count("G") + probe_seq[:25].count("C")) / len(probe_seq[:25]) if len(probe_seq[:25]) > 0 else 0.5,
                            expected_blast_hits = 1,
                            total_rhs_blast_hits = [],
                            total_lhs_blast_hits = [],
                            rhs_blast_hits = [],
                            lhs_blast_hits = [],
                            score = 0,
                            transcript_sequence = sequence,
                            gap_length = None,
                            gap_probe_sequence = None,
                            gap_gene_sequence = None,
                            target_start_gap = None,
                            target_end_gap = None,
                            original_transcript_sequence = sequence,
                            original_target_start_gap = None,
                            original_target_end_gap = None,
                        )
        for i, (transcript_name, transcript_sequence) in tqdm(enumerate(transcript_sequences.items()),
                                                              desc='Generating probes',
                                                              total=len(transcript_sequences),
                                                              unit='target',
                                                              dynamic_ncols=True):
            if isinstance(transcript_sequence, tuple):
                original_transcript_sequence, transcript_sequence = transcript_sequence
            else:
                original_transcript_sequence = None

            target_start = target_starts[i]
            if isinstance(target_start, tuple):
                original_target_start, target_start = target_start
            else:
                original_target_start = None

            target_end = target_ends[i]
            if isinstance(target_end, tuple):
                original_target_end, target_end = target_end
            else:
                original_target_end = None

            expected = expect_hits[i]

            probes = self.create_gapfilling_flex_probes(
                transcript_name,
                transcript_sequence,
                target_start,
                target_end,
                original_transcript_sequence,
                original_target_start,
                original_target_end,
                expected,
                n_probes,
                max_iterations_per_probe,
                initial_temp,
                existing_probes=all_probes,
                flex_probes=flex_probes,
                fast=fast,
                predefined_lhs=lhs_probes[i] if lhs_probes else None
            )
            if probes is None or len(probes) == 0:
                missing.append(transcript_name)
                all_probes[transcript_name] = []
            else:
                all_probes[transcript_name] = probes

        if len(missing) > 0:
            print(f"Warning: No probes for")
            for name in missing:
                print(f"-{name}")

        # If we are running with approximate search, we should filter out probes that have mismatches
        if fast:
            lhs_seqs = [probe.lhs_gene_sequence for probes in all_probes.values() for probe in probes]
            rhs_seqs = [probe.rhs_gene_sequence for probes in all_probes.values() for probe in probes]
            expect_hits = None if not expect_hits else [hits for hits, probes in zip(expect_hits, all_probes.values()) for _ in probes]
            to_remove, all_hits, all_names = self.bulk_blast_filter_probes(lhs_seqs, rhs_seqs, expect_hits)
            print(f"Removing {len(to_remove)} probes due to off-target activity.")
            curr_idx = 0
            for transcript_name, probes in list(all_probes.items()):
                removed = []
                for probe in probes:
                    if curr_idx in to_remove:
                        removed.append(probe)
                    else:
                        probe.total_lhs_blast_hits = all_hits[curr_idx][0]
                        probe.total_rhs_blast_hits = all_hits[curr_idx][1]
                        probe.expected_blast_hits = expect_hits[curr_idx]
                        probe.lhs_blast_hits = all_names[curr_idx][0]
                        probe.rhs_blast_hits = all_names[curr_idx][1]
                    curr_idx += 1
                all_probes[transcript_name] = [probe for probe in probes if probe not in removed]

        if visium:
            return self.convert_probe_set_to_df(all_probes, truseq=truseq, barcode=barcode, visium=True)
        return self.convert_probe_set_to_df(all_probes, truseq=truseq, barcode=barcode)

    def generate_expect_hits(self, sequences: list[str]) -> list[int]:
        print("Searching for baseline hits to the transcriptome.")
        expect_hits = []
        for transcript_sequence in tqdm(sequences, desc='Baseline hits', unit='target', dynamic_ncols=True):
            n_hits, _ = self.blast_hits(transcript_sequence)
            expect_hits.append(n_hits)
        return expect_hits
