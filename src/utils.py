import itertools
import os
import subprocess
from collections import namedtuple, Counter
from io import StringIO
from typing import Union
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Blast import NCBIXML
from pyensembl import Genome, Transcript, Exon
from Bio import SeqIO
import niquests as requests

Probe = namedtuple("Probe", [
    "name",
    'strand',
    "lhs",
    "rhs",
    "lhs_start",
    "rhs_start",
    "bridge",
    "lhs_gene",
    "rhs_gene",
    "bridge_gene",
    "bridge_length",
    "score",
    "lhs_GC_content",
    "rhs_GC_content",
    "bridge_snv_start",
    "bridge_snv_end"
])


# From table 2: https://kb.10xgenomics.com/hc/en-us/articles/17623693026445-How-do-I-design-custom-probes-for-a-Single-Cell-Gene-Expression-Flex-i-e-Fixed-RNA-Profiling-for-multiplexed-samples-experiment
barcode_seqs = [
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


def fetch_csv_from_url(filename: str, url: str, force_download: bool = False) -> pd.DataFrame:
    """
    Fetch a CSV from `url` and cache it to `filename`.
    """
    filename = Path(filename)
    if force_download or not filename.exists():
        filename.parent.mkdir(parents=True, exist_ok=True)
        response = requests.get(url)
        with open(filename, "wb") as f:
            f.write(response.content)
    return pd.read_csv(filename, comment='#')


def fetch_human_flex_v2_probeset(filename: str = "./human_v2_probes.csv") -> pd.DataFrame:
    """
    Fetch the Flex v2 probe set from 10X Genomics.
    """
    url = "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv"
    return fetch_csv_from_url(filename, url)


def fetch_mouse_flex_v2_probeset(filename: str = "./mouse_v2_probes.csv") -> pd.DataFrame:
    """
    Fetch the Flex v2 probe set from 10X Genomics.
    """
    url = "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv"
    return fetch_csv_from_url(filename, url)


def fetch_human_flex_v1_probeset(filename: str = "./human_v1_probes.csv") -> pd.DataFrame:
    """
    Fetch the Flex v1 probe set from 10X Genomics.
    """
    url = "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
    return fetch_csv_from_url(filename, url)


def fetch_mouse_flex_v1_probeset(filename: str = "./mouse_v1_probes.csv") -> pd.DataFrame:
    """
    Fetch the Flex v1 probe set from 10X Genomics.
    """
    url = "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.1_GRCm39-2024-A.csv"
    return fetch_csv_from_url(filename, url)


def fetch_human_visiumhd_probeset(filename: str = "./human_visiumhd_probes.csv") -> pd.DataFrame:
    """
    Fetch the Visium HD probe set from 10X Genomics.
    """
    url = "https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv"
    return fetch_csv_from_url(filename, url)


def fetch_mouse_visiumhd_probeset(filename: str = "./mouse_visiumhd_probes.csv") -> pd.DataFrame:
    """
    Fetch the Visium HD probe set from 10X Genomics.
    """
    url = "https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv"
    return fetch_csv_from_url(filename, url)


def parse_snv_info(dna_snv_string: str) -> tuple[int, int, str, str, int]:
    """
    Parse the SNV information from the input strings.
    :param dna_snv_string: The DNA SNV string (ex. c.1849G>T).
    :param hg38_location_string: The hg38 location string (ex. 9:133256215-133256215).
    :return: The SNV start, SNV end, action, and data.
    """
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
            start = end = int(dna_range[0][:-1])
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
            start = end = int(dna_range[0])
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
            start = end = int(dna_range[0])
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
            start = end = int(dna_range[0])
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
            start = end = int(dna_range[0])
        else:
            start = int(dna_range[0])
            end = int(dna_range[1])
        data = f"inv{split[1]}"
        snv_length_change = 0
    else:
        print(f"WARNING: Unrecognized SNV action: {dna_snv_string}, assuming genomic range.")
        split = dna_snv_string.split("_")
        if len(split) == 1:
            start = end = int(split[0])
        else:
            start = int(split[0])
            end = int(split[1])
        data = "undefined"
        action = "undefined"
        snv_length_change = 0

    return start, end, action, data, snv_length_change


def parse_chromosome_location(hg38_location_string: str) -> tuple[str, int, int]:
    """
    Parse the chromosome location from the input string.
    :param hg38_location_string: The hg38 location string (ex. 9:133256215-133256215).
    :return: The chromosome, hg38 start, and hg38 end.
    """
    hg38_location_string = hg38_location_string.replace(" ", "").replace("..", "-")
    split_string = hg38_location_string.split(":")
    chromosome = split_string[0]
    location_range = split_string[1].split("-")
    hg38_start = int(location_range[0])
    hg38_end = int(location_range[1]) if len(location_range) == 2 else hg38_start

    return chromosome, hg38_start, hg38_end


def get_intron(transcript: Transcript, previous_exon: Exon = None, next_exon: Exon = None) -> str:
    intron_start = previous_exon.end + 1 if previous_exon is not None else 0
    intron_end = next_exon.start - 1 if next_exon is not None else transcript.length
    genome: Genome = transcript.genome
    genome.db.query_locus(transcript.contig, intron_start, intron_end)
    return genome.db.get_sequence()


def max_homopolymer_length(sequence: str) -> int:
    if len(sequence) == 0:
        return 0
    groups = itertools.groupby(sequence)
    return max(len(list(group)) for key, group in groups)

def get_melting_temp(sequence: str):
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    melting_temp = 64.9 + 41*(g_count+c_count-16.4)/(a_count+t_count+g_count+c_count)
    return melting_temp

def setup_blast(blast_db: str, sequences: dict[str, str]):
    """
    Set up the BLAST database for a genome. Note that the name of each sequence should begin with the gene name.
    :param sequences: The sequences to add to the BLAST database. Keys are the names of the sequences. Values are the sequences.
    """
    with open("temp.fasta", "w") as fasta_file:
        for name, sequence in sequences.items():
            fasta_file.write(f">{name}\n{sequence}\n")

    subprocess.run(["makeblastdb", "-in", "temp.fasta", "-dbtype", "nucl", "-out", str(blast_db)])
    os.remove("temp.fasta")


def blast_search(blast_db: str,
                 sequence: str | list[str],
                 evalue: float = 1.0) -> Union['NCBIXML', 'list[NCBIXML]']:
    """
    Perform a BLAST search to identify the number of genes this matches to.
    :param blast_db: The BLAST database path to search.
    :param sequence: The sequence(s) to search.
    :param evalue: The maximum E-value to consider.
    :return: The BLAST results.
    """
    multi_sequences = not isinstance(sequence, str)
    # BLAST drops a '*' wildcard from a query, which shifts every later position.
    # 'N' keeps the position and aligns as one mismatch.
    if multi_sequences:
        fasta = "\n".join([f">{i}\n{seq.replace('*', 'N')}" for i, seq in enumerate(sequence)])
    else:
        fasta = sequence.replace('*', 'N')
    blastn_command = ["blastn", "-db", str(blast_db), "-evalue", str(evalue), "-outfmt", "5", "-num_threads", "4", "-task", "blastn-short"]
    blast_out = subprocess.run(blastn_command, input=fasta, text=True, stdout=subprocess.PIPE).stdout
    if len(blast_out) == 0:
        return None

    # Parse the BLAST results
    if multi_sequences:
        blast_hits = list(sorted(NCBIXML.parse(StringIO(blast_out)), key=lambda rec: int(rec.query)))
    else:
        blast_hits = NCBIXML.read(StringIO(blast_out))

    return blast_hits


def probe_blast_search(blast_db: str,
                       name: str,
                       lhs_sequence: str,
                       rhs_sequence: str,
                       evalue: float = 1.0
                       ) -> tuple[int, int]:
    """
    Perform a BLAST search to identify the number of genes this matches to.
    :param blast_db: The BLAST database path to search.
    :param name: The name of the gene the probe is targeting.
    :param lhs_sequence: The left-hand side probe sequence.
    :param rhs_sequence: The right-hand side probe sequence.
    :param evalue: The maximum E-value to consider.
    :return: The number of hits for the LHS and RHS probes.
    """
    lhs_hits = blast_search(blast_db, lhs_sequence, evalue)
    rhs_hits = blast_search(blast_db, rhs_sequence, evalue)

    if lhs_hits is None or rhs_hits is None:
        return 0, 0

    # 10X recommends at least 5 mismatches for a probe to be considered unique, so select alignments with < 5 mismatches to be "hits"
    # Also compute hits wrt gene name only
    n_lhs_hits = len(set([aln.hit_def.split(" ")[0] for aln in lhs_hits.alignments if max(hsp.identities for hsp in aln.hsps) >= len(lhs_sequence) - 5]))
    n_rhs_hits = len(set([aln.hit_def.split(" ")[0] for aln in rhs_hits.alignments if max(hsp.identities for hsp in aln.hsps) >= len(rhs_sequence) - 5]))

    return n_lhs_hits, n_rhs_hits


class MismatchScanner:
    """
    Find every place in a background where a short query matches with at most a given number of
    mismatches, on either strand. Unlike BLAST, the search is exact: it misses no such match.

    An 'N' (in either case) is ambiguous, so it always counts as a match, in the query and in the
    background. Any other non-ACGT character, including the other IUPAC codes, counts as a mismatch
    unless it faces an 'N'.

    It uses the pigeonhole principle. A query with at most m mismatches has at most m // 2 mismatches
    in its first k bases or in its last k bases (2k <= query length). Every variant of these two seeds
    with at most m // 2 substitutions is looked up in one pass over the background, and each candidate
    is then checked base by base over the whole query. Insertions and deletions are not counted as matches.
    """
    _N = 4  # An ambiguous base. It matches every base.
    _OTHER = 5  # Any other character. It matches no base.
    _SEPARATOR = 6  # Between two background sequences
    _PAD = 7  # Fills a query up to the length of the longest query
    _ENCODE = np.full(256, _OTHER, dtype=np.uint8)
    for _i, _base in enumerate("ACGTN"):
        _ENCODE[ord(_base)] = _i
        _ENCODE[ord(_base.lower())] = _i
    _COMPLEMENT = np.array([3, 2, 1, 0, _N, _OTHER], dtype=np.uint8)
    # A background seed window with more non-ACGT bases than this is checked against every query,
    # rather than looked up once for each base its non-ACGT bases could be
    _MAX_EXPANDED = 4

    def __init__(self, sequences: dict[str, str], seed_length: int = 12, chunk_size: int = 1 << 24):
        """
        :param sequences: The background. Key is the gene name, a space, and the transcript id; value is the sequence.
        :param seed_length: The maximal seed length. 12 keeps one lookup table at 16 MB.
        :param chunk_size: The number of background positions handled at once.
        """
        self.names = list(sequences.keys())
        self.genes = [name.split(" ")[0] for name in self.names]
        text = "#".join(sequences[name] for name in self.names).encode("ascii", errors="replace")
        self.background = self._ENCODE[np.frombuffer(text, dtype=np.uint8)]
        lengths = np.array([len(sequences[name]) for name in self.names], dtype=np.int64)
        self.starts = np.concatenate([[0], np.cumsum(lengths + 1)[:-1]]).astype(np.int64)
        self.ends = self.starts + lengths
        # One separator (never a match) between sequences, so no window crosses two of them
        self.background[self.ends[:-1]] = self._SEPARATOR
        self.seed_length = seed_length
        self.chunk_size = chunk_size

    def _encode_query(self, query: str) -> np.ndarray:
        return self._ENCODE[np.frombuffer(query.encode("ascii", errors="replace"), dtype=np.uint8)]

    @staticmethod
    def _seed_variants(seed: np.ndarray, max_mismatches: int) -> np.ndarray:
        """
        All codes of the seed with at most max_mismatches substitutions. A non-ACGT base takes every base
        at no cost. This is exact for an 'N' and a superset for any other base, which the full check removes.
        """
        k = len(seed)
        wild = [i for i in range(k) if seed[i] > 3]
        known = [i for i in range(k) if seed[i] <= 3]
        weights = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)
        base_code = int(np.sum(np.where(seed > 3, 0, seed).astype(np.int64) * weights))
        deltas = []
        for n_extra in range(max_mismatches + 1):
            for extra in itertools.combinations(known, n_extra):
                choices = [[b for b in range(4) if b != seed[p]] for p in extra]
                for bases in itertools.product(*choices):
                    deltas.append(sum((b - int(seed[p])) * int(weights[p]) for p, b in zip(extra, bases)))
        if len(deltas) * 4 ** len(wild) >= 4 ** k // 4:
            return np.arange(4 ** k, dtype=np.int64)  # Close to every code, so take them all
        codes = base_code + np.array(deltas, dtype=np.int64)
        for p in wild:
            codes = (codes[:, None] + np.arange(4, dtype=np.int64)[None, :] * weights[p]).ravel()
        return np.unique(codes)

    def search(self, queries: list[str], max_mismatches: int = 4) -> list[dict[str, int]]:
        """
        :param queries: The query sequences, in any case. Each must be at least 2 bases long.
        :param max_mismatches: The highest number of mismatches that still counts as a match.
        :return: For each query, the genes it matches, each with its lowest number of mismatches.
        """
        if len(queries) == 0:
            return []
        # Both strands: the reverse complement of a query is searched as one more query
        encoded = [self._encode_query(q) for q in queries]
        encoded += [self._COMPLEMENT[q[::-1]] for q in encoded]
        lengths = np.array([len(q) for q in encoded], dtype=np.int64)
        k = int(min(self.seed_length, lengths.min() // 2))
        seed_mismatches = max_mismatches // 2
        weights = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)

        # Seed table: code -> (query, offset of the seed in the query)
        key_codes, key_query, key_offset = [], [], []
        for qi, q in enumerate(encoded):
            for offset in (0, len(q) - k):
                variants = self._seed_variants(q[offset:offset + k], seed_mismatches)
                key_codes.append(variants)
                key_query.append(np.full(len(variants), qi, dtype=np.int64))
                key_offset.append(np.full(len(variants), offset, dtype=np.int64))
        key_codes = np.concatenate(key_codes)
        order = np.argsort(key_codes, kind="stable")
        key_codes, key_query, key_offset = key_codes[order], np.concatenate(key_query)[order], np.concatenate(key_offset)[order]
        unique_codes, first_key, key_counts = np.unique(key_codes, return_index=True, return_counts=True)
        in_table = np.zeros(4 ** k, dtype=bool)
        in_table[unique_codes] = True
        # Every (query, seed offset) pair, for windows that are checked against every query
        pair_query = np.repeat(np.arange(len(encoded), dtype=np.int64), 2)
        pair_offset = np.stack([np.zeros(len(encoded), dtype=np.int64), lengths - k], axis=1).ravel()

        max_len = int(lengths.max())
        padded = np.full((len(encoded), max_len), self._PAD, dtype=np.uint8)  # Masked below
        for qi, q in enumerate(encoded):
            padded[qi, :len(q)] = q

        def candidates(seed_positions: np.ndarray, seed_codes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            """Expand each seed hit to every (query, offset) that owns its code. Return queries and query starts."""
            slot = np.searchsorted(unique_codes, seed_codes)
            counts = key_counts[slot]
            positions = np.repeat(seed_positions, counts)
            first = np.repeat(first_key[slot], counts)
            within = np.arange(len(positions)) - np.repeat(np.cumsum(counts) - counts, counts)
            keys = first + within
            return key_query[keys], positions - key_offset[keys]

        best = [dict() for _ in queries]
        n = len(self.background)
        for chunk_start in range(0, n - k + 1, self.chunk_size):
            chunk_end = min(chunk_start + self.chunk_size, n - k + 1)
            window = self.background[chunk_start:chunk_end + k - 1]
            m = chunk_end - chunk_start
            codes = np.zeros(m, dtype=np.int64)
            n_ambiguous = np.zeros(m, dtype=np.int64)  # 'N' and other non-ACGT bases
            crosses = np.zeros(m, dtype=bool)
            for j in range(k):
                column = window[j:j + m]
                codes = codes * 4 + np.where(column <= 3, column, 0)
                n_ambiguous += (column == self._N) | (column == self._OTHER)
                crosses |= column == self._SEPARATOR
            found_query, found_start = [], []

            # Seed windows of ACGT only: one lookup each
            hits = np.flatnonzero((n_ambiguous == 0) & ~crosses & in_table[codes])
            q, s = candidates(hits + chunk_start, codes[hits])
            found_query.append(q)
            found_start.append(s)

            # Seed windows with a few non-ACGT bases: one lookup for each base they could be. For an 'N'
            # this keeps it a match. Any other base is a mismatch whatever it is replaced with, so the
            # replacement that agrees with the query seed counts it once, as the full check below does.
            ambiguous = np.flatnonzero((n_ambiguous > 0) & ~crosses)
            for c in range(1, self._MAX_EXPANDED + 1):
                rows = ambiguous[n_ambiguous[ambiguous] == c]
                if len(rows) == 0:
                    continue
                block = window[rows[:, None] + np.arange(k)[None, :]]
                variable = np.nonzero(block > 3)[1].reshape(len(rows), c)
                combos = np.array(list(itertools.product(range(4), repeat=c)), dtype=np.int64)
                expanded = codes[rows][:, None] + weights[variable] @ combos.T
                row_idx, combo_idx = np.nonzero(in_table[expanded])
                q, s = candidates(rows[row_idx] + chunk_start, expanded[row_idx, combo_idx])
                found_query.append(q)
                found_start.append(s)

            # Seed windows with many non-ACGT bases (rare): check against every query
            many = ambiguous[n_ambiguous[ambiguous] > self._MAX_EXPANDED]
            found_query.append(np.tile(pair_query, len(many)))
            found_start.append(np.repeat(many + chunk_start, len(pair_query)) - np.tile(pair_offset, len(many)))

            query_idx, starts = np.concatenate(found_query), np.concatenate(found_start)
            if len(starts) == 0:
                continue
            # The whole query must fall inside one background sequence
            seq_idx = np.searchsorted(self.starts, starts, side="right") - 1
            inside = (starts >= 0) & (seq_idx >= 0)
            inside[inside] &= starts[inside] + lengths[query_idx[inside]] <= self.ends[seq_idx[inside]]
            query_idx, starts, seq_idx = query_idx[inside], starts[inside], seq_idx[inside]
            for batch in range(0, len(starts), 1 << 20):
                qb, sb, tb = query_idx[batch:batch + (1 << 20)], starts[batch:batch + (1 << 20)], seq_idx[batch:batch + (1 << 20)]
                target = self.background[sb[:, None] + np.minimum(np.arange(max_len)[None, :], lengths[qb][:, None] - 1)]
                query = padded[qb]
                match = ((target == query) & (query <= 3)) | (query == self._N) | (target == self._N)
                mismatches = (~match & (query != self._PAD)).sum(axis=1)
                keep = mismatches <= max_mismatches
                for qi, ti, mm in zip(qb[keep], tb[keep], mismatches[keep]):
                    query_number = int(qi) % len(queries)
                    gene = self.genes[ti]
                    if mm < best[query_number].get(gene, max_mismatches + 1):
                        best[query_number][gene] = int(mm)
        return best


def rank_and_filter_transcripts(transcripts: list[Transcript], use_support=True, filter_transcripts=True) -> list[Transcript]:
    """
    Rank the transcripts by support level and length. Additionally only retain complete protein coding transcripts.
    :param transcripts: The transcripts to rank.
    :param use_support: Whether to use the support level in ranking.
    :param filter_transcripts: Whether to filter out incomplete transcripts.
    :return: The ranked transcripts.
    """
    filtered = list(filter(lambda x: x.complete or not filter_transcripts, sorted(transcripts, key=lambda x: (x.length, -(x.support_level or 0) if use_support else 1), reverse=True)))
    return filtered


def transcriptome(genome: Genome, filter_incomplete: bool = True, filter_biotypes: list[str] = None) -> dict[str, str]:
    """
    Select the largest transcript for each gene to create a non-redundant transcriptome.
    :param genome: The genome to work with.
    :param filter_incomplete: Whether to filter out incomplete transcripts.
    :param filter_biotypes: The biotypes to filter out.
    :return: The transcriptome. Key is gene name and transcript ID, value is the full transcript sequence.
    """
    transcriptome = dict()
    unnamed_count = 0
    for gene in genome.genes():
        try:
            transcripts = rank_and_filter_transcripts(gene.transcripts, filter_transcripts=filter_incomplete)
        except Exception:
            continue  # Skip genes with annotation issues
        if len(transcripts) == 0:
            continue
        for transcript in transcripts:
            gene_name = gene.gene_name.strip() or gene.name.strip() or gene.gene_id.strip()
            if gene_name == "":
                gene_name = f"Unnamed{unnamed_count}"
                unnamed_count += 1
            if filter_biotypes is not None and transcript.biotype in filter_biotypes:
                continue
            try:
                # Use the full transcript (including UTRs). Probes bind the mature mRNA, so an
                # off-target match in an untranslated region is as real as one in a coding region.
                sequence = transcript.sequence
            except Exception:
                continue  # Skip transcripts with annotation issues (e.g., duplicate exon IDs)
            if sequence is not None:
                transcriptome[gene_name + " " + transcript.transcript_id] = sequence
    return transcriptome


def export_probes_to_xlsx(probes: list[Probe], filename: str, probe_order_format: bool = True, n_barcodes: int = 1) -> pd.DataFrame:
    results = dict(
        name=[],
        strand=[],
        lhs_probe=[],
        rhs_probe=[],
        lhs_seq=[],
        rhs_seq=[],
        lhs_gene=[],
        rhs_gene=[],
        lhs_GC_content=[],
        rhs_GC_content=[],
        GC_difference=[],
        gap_filled=[],
        gap_filled_gene=[],
        gap_length=[],
        gap_snv_start=[],
        gap_snv_end=[],
        score=[],
    )
    # Note that the LHS
    rhs_prefix = "/5Phos/"
    rhs_probe_barcode_bridge = "ACGCGGTTAGCACGTANN"
    rhs_suffix = "CGGTCCTAGCAA"

    lhs_prefix = "CCTTGGCACCCGAGAATTCCA"

    for probe in probes:
        if probe is None:
            continue
        for barcode in barcode_seqs[:n_barcodes]:
            results["name"].append(probe.name)
            results["strand"].append(probe.strand)
            results["lhs_probe"].append(lhs_prefix + probe.lhs)
            results["rhs_probe"].append(rhs_prefix + probe.rhs + rhs_probe_barcode_bridge + barcode + rhs_suffix)
            results["lhs_seq"].append(probe.lhs)
            results["rhs_seq"].append(probe.rhs)
            results["lhs_gene"].append(probe.lhs_gene)
            results["rhs_gene"].append(probe.rhs_gene)
            results["lhs_GC_content"].append(probe.lhs_GC_content)
            results["rhs_GC_content"].append(probe.rhs_GC_content)
            results["GC_difference"].append(abs(probe.lhs_GC_content - probe.rhs_GC_content))
            results["gap_filled"].append(probe.bridge if probe.bridge_length > 0 else "N/A")
            results["gap_filled_gene"].append(probe.bridge_gene if probe.bridge_length > 0 else "N/A")
            results["gap_length"].append(probe.bridge_length)
            results["gap_snv_start"].append(probe.bridge_snv_start)
            results["gap_snv_end"].append(probe.bridge_snv_end)
            results["score"].append(probe.score)

    df = pd.DataFrame(results)
    if not probe_order_format:
        df.to_excel(filename, index=False)
    else:
        exported_df = dict(
            name=[],
            value=[]
        )
        for i, row in df.iterrows():
            for j, barcode in enumerate(barcode_seqs[:n_barcodes]):
                exported_df['name'].append(row['name'] + " " + f"#{barcode}")
                # exported_df['value'].append("3'-5'" if row['strand'] == '-' else "5'-3'")
                exported_df['value'].append("5'-3'")

                exported_df["name"].append("RHSprobe")
                exported_df["value"].append(rhs_prefix + row['rhs_seq'] + rhs_probe_barcode_bridge + barcode + rhs_suffix)

                exported_df["name"].append("LHSprobe")
                exported_df["value"].append(lhs_prefix + row['lhs_seq'])
        pd.DataFrame(exported_df).to_excel(filename, index=False, header=False)
    return df

def generate_seq_dict(reference):
    seq_dict = {}
    for long_sequence_record in SeqIO.parse(open(reference), 'fasta'):
        id = long_sequence_record.description
        long_sequence_record = long_sequence_record.reverse_complement()
        long_sequence_record.id = id
        seq_dict[long_sequence_record.id] = long_sequence_record
    return seq_dict

def check_overlap_by_alignment(seq1,seq2,long_sequence=None,seq_dict=None):
    if seq_dict is not None:
        for seq_id, long_sequence_record in seq_dict.items():
            long_sequence = str(long_sequence_record.seq)
            if seq1 in long_sequence:
                if seq2 in long_sequence:
                    ## check for overlap: is end of seq1 before start of seq2 or vice versa
                    if (long_sequence.index(seq1) + len(seq1) <= long_sequence.index(seq2)) or (long_sequence.index(seq2) + len(seq2) <= long_sequence.index(seq1)):
                        continue
                    else:
                        return True
    elif long_sequence is not None:
        if seq1 in long_sequence:
            if seq2 in long_sequence:
                ## check for overlap: is end of seq1 before start of seq2 or vice versa
                if (long_sequence.index(seq1) + len(seq1) <= long_sequence.index(seq2)) or (long_sequence.index(seq2) + len(seq2) <= long_sequence.index(seq1)):
                    pass
                else:
                    return True
    else:
        raise ValueError("Either long_sequence or seq_dict must be provided")
    return False

def find_overlaps_with_flex(lhs_probe, rhs_probe, long_sequence=None, seq_dict=None, flex_list=None,return_list=False):
    if flex_list is None:
        raise ValueError("flex_list must be provided")
    if seq_dict is not None:
        seq_dict = {seq_id: long_sequence_record for seq_id, long_sequence_record in seq_dict.items() if (lhs_probe in str(long_sequence_record.seq)) or (rhs_probe in str(long_sequence_record.seq))}
        lhs_overlaps = flex_list['probe_seq'].apply(lambda seq2: check_overlap_by_alignment(lhs_probe, seq2, seq_dict = seq_dict))
        rhs_overlaps = flex_list['probe_seq'].apply(lambda seq2: check_overlap_by_alignment(rhs_probe, seq2, seq_dict = seq_dict))
    elif long_sequence is not None:
        lhs_overlaps = flex_list['probe_seq'].apply(lambda seq2: check_overlap_by_alignment(lhs_probe, seq2, long_sequence = long_sequence))
        rhs_overlaps = flex_list['probe_seq'].apply(lambda seq2: check_overlap_by_alignment(rhs_probe, seq2, long_sequence = long_sequence))
    else:
        raise ValueError("Either long_sequence or seq_dict must be provided")
    overlap_list = lhs_overlaps | rhs_overlaps
    if return_list:
        return overlap_list
    else:
        return overlap_list.sum()


def has_tandem_repeat(sequence):
    n_repeat = 0
    # Find all substrings of length 3 or more
    for start in [0,1,2]:
        substrings = [sequence[i:i + length] for length in range(3, len(sequence) + 1) for i in range(start, len(sequence) - length + 1, length)]
        # Count occurrences of each substring
        substring_counts = Counter(substrings)
        n_repeat = max(n_repeat,max(substring_counts.values()))
    if n_repeat >= 4:
        return n_repeat
    return 0
