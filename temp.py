import pandas as pd
import os

os.environ["PYENSEMBL_CACHE_DIR"] = "./ensembl_cache"

from pyensembl import EnsemblRelease
from pyensembl.species import Species, human
from difflib import SequenceMatcher


genome = None

def find_best_ensembl_transcript(
        gene_name: str,
        transcript_sequence: str,
        ensembl_release: int = 111,
        species: Species = human,
        min_similarity: float = 0.0
) -> tuple[str | None, float]:
    """
    Find the Ensembl transcript ID that best matches a given sequence for a gene.

    Args:
        gene_name: The gene symbol (e.g., "TP53", "BRCA1")
        transcript_sequence: The transcript sequence to match against
        ensembl_release: Ensembl release version (default: 111)
        species: Species for genome lookup (default: human)
        min_similarity: Minimum similarity ratio to consider a match (0.0-1.0)

    Returns:
        Tuple of (best_transcript_id, similarity_ratio)
        Returns (None, 0.0) if no transcripts found or none meet min_similarity
    """
    # Build/load genome
    global genome
    if genome is None:
        genome = EnsemblRelease(release=ensembl_release, species=species)
        genome.download(overwrite=False)
        genome.index(overwrite=False)

    # Get all transcript IDs for the gene
    try:
        transcript_ids = genome.transcript_ids_of_gene_name(gene_name)
    except ValueError:
        # Gene not found
        return None, 0.0

    if not transcript_ids:
        return None, 0.0

    # Normalize input sequence
    query_seq = transcript_sequence.upper().replace(" ", "").replace("\n", "")

    best_transcript_id = None
    best_similarity = 0.0

    for transcript_id in transcript_ids:
        try:
            transcript = genome.transcript_by_id(transcript_id)
            try:
                ensembl_seq = transcript.coding_sequence
            except:
                ensembl_seq = transcript.sequence
        except Exception:
            # Some transcripts may not have sequence data
            continue

        if ensembl_seq is None:
            continue

        # Check for exact match first (fastest)
        if query_seq == ensembl_seq:
            return transcript_id, 1.0

        # Calculate similarity ratio
        similarity = SequenceMatcher(None, query_seq, ensembl_seq).ratio()

        if similarity > best_similarity:
            best_similarity = similarity
            best_transcript_id = transcript_id

    if best_similarity < min_similarity:
        return None, best_similarity

    return best_transcript_id, best_similarity


def write_new_inputs(inputs):
    inputs = pd.read_csv(inputs, header=None)
    inputs.columns = ['name']
    with open("new_inputs.csv", 'w') as f:
        for i, row in inputs.iterrows():
            try:
                transcript = merged[row['name']]
                f.write(transcript)
                f.write(" ")
                f.write(" ".join(row['name'].split(" ")[1:]))
                f.write("\n")
            except:
                f.write(row['name'])
                f.write("\n")


def compare_transcripts(new, old, release=111):
    new = pd.read_table(new)
    if old is None:
        old = new.copy()
        old = old[['name', 'transcript']]
        old['name'] = old['name'] + "_"
    else:
        old = pd.read_csv(old, header=None)
    new['gene'] = new.name.str.split(" ").str[0]
    old.columns = ["gene", "transcript"]
    merged = pd.merge(new, old, on=['gene'], how='left', suffixes=("_new", "_old"))
    merged['matches'] = merged['transcript_new'] == merged['transcript_old']
    probe2enst = dict()
    missing = 0
    imperfect = 0
    for i, row in merged.iterrows():
        if row.matches:
            probe2enst[row['name']] = row['transcript_id']
        else:
            best_transcript, similarity = find_best_ensembl_transcript(row['gene'],
                                                                       row['transcript_old'] if not pd.isna(
                                                                           row['transcript_old']) else row[
                                                                           'transcript_new'],
                                                                       ensembl_release=release)
            if best_transcript is None:
                print("WARNING CANT FIND FOR", row['name'])
                missing += 1
            elif similarity < 1.:
                imperfect += 1
            else:
                probe2enst[row['name']] = best_transcript
    print("MISSING", missing, "IMPERFECT", imperfect)
    return probe2enst


merged = compare_transcripts("test_probes.tsv", None, #"../gapped_flex_probe_design/v2_MPN_Patient_Probes/gene2transcript.csv",
                             111)

input()

write_new_inputs("probe_sets/impact/inputs.csv")