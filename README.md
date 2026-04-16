# GIFT Probe Design Pipeline

This repository contains the pipeline for designing custom gap-filling probes compatible with 10x Genomics Flex (single-cell gene expression) and Visium HD (spatial transcriptomics) platforms. Probes are designed to target specific single nucleotide variants (SNVs) in transcripts using HGVSc notation.

## Requirements

- Python 3.12+
- NCBI BLAST+ command-line tools (`makeblastdb`, `blastn`) — only required when using `--blast` for off-target filtering

Note: The repo also supports installing dependencies with [pixi](https://pixi.prefix.dev/latest/) which will automatically
install python and blast into a virtual environment.

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd probeset_design_pipeline

# Install dependencies using uv or pip
uv sync
# or
pip install -e .

# Install dependencies with pixi (includes python and blast for --blast support)
pixi install
```

On first run, the pipeline will download Ensembl genome annotations to `./ensembl_cache/`. This may take several minutes.

## Usage

```bash
python src/app.py [OPTIONS] <targets> <output_name>
```

### Arguments

| Argument | Description |
|----------|-------------|
| `targets` | Path to CSV file containing target variants |
| `output_name` | Base name for output files (without extension) |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--config_file` | None | Path to JSON configuration file |
| `--organism` | human | Target organism: `human` or `mouse` |
| `--technology` | FlexV1 | Platform: `FlexV1`, `FlexV2`, or `VisiumHD` |
| `--barcodes` | 1 | Number of barcodes for multiplexing (Flex v1 only, max 16) |
| `--output_format` | tsv | Output format: `csv`, `tsv`, or `xlsx` |
| `--skip_errors` | false | Skip targets that fail (gene not found, invalid variant) instead of aborting |
| `--mane` | false | Use MANE Select canonical transcripts (recommended for reproducibility) |
| `--msk` | false | Use MSK IMPACT canonical transcript overrides |
| `--release` | 111 | Ensembl release version to use |
| `--blast` | false | Enable BLAST-based off-target filtering (requires BLAST+ in PATH) |

### Example

```bash
python src/app.py --mane --skip_errors --config_file probe_sets/GBM/config.json probe_sets/GBM/inputs.csv GBM_probes
```

## Input Format

The targets file is a CSV without a header row. Three formats are supported:

**Single column** - HGVSc notation only (gene name extracted from first word):
```
EGFR c.2582T>G
TP53 c.524G>A
PIK3CA c.1633G>A
```

**Two columns** - HGVSc notation and custom transcript sequence:
```
EGFR c.2582T>G,ATGCGA...
```

**Three columns** - Gene name, HGVSc notation, and transcript sequence:
```
EGFR,c.2582T>G,ATGCGA...
```

### Supported Variant Types

- Point mutations: `c.2582T>G`
- Deletions: `c.227_228del`
- Insertions: `c.1791dup`, `c.100_101insACG`
- Deletion-insertions: `c.1889_1894delinsTAGGAT`
- Inversions: `c.100_105inv`

Use `0bp` for wildtype probes targeting a gene without a specific variant (e.g., `PIK3R1 0bp`).

Ensembl transcript IDs are also accepted in place of gene symbols (e.g., `ENST00000275493 c.2582T>G`).

## Transcript Selection

The pipeline fetches transcript sequences from Ensembl. Transcript isoform selection significantly affects probe sequences since HGVSc positions are relative to the coding sequence.

| Mode | Flag | Description |
|------|------|-------------|
| Default | (none) | Selects longest complete transcript for each gene |
| MANE Select | `--mane` | Uses NCBI/Ensembl agreed-upon canonical transcripts (recommended) |
| MSK IMPACT | `--msk` | Uses MSK-IMPACT canonical transcript overrides |
| Explicit | (none) | Use Ensembl transcript ID instead of gene symbol |

**Recommendation**: Use `--mane` for reproducibility. MANE Select provides stable, standardized transcript selection agreed upon by NCBI and Ensembl.

## Configuration

Configuration files control probe design parameters. See `probe_sets/GBM/config.json` for an example.

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_bridge_length` | 10 | Maximum gap/bridge length in bp |
| `min_bridge_length` | 1 | Minimum gap/bridge length in bp |
| `tx_min_gc` | 0.2 | Minimum GC content for probes |
| `tx_max_gc` | 0.8 | Maximum GC content for probes |
| `strict_gc_content` | true | Reject probes outside GC bounds (vs. penalize) |
| `add_probes_to_blast` | true | Add accepted probes to BLAST database for off-target checking |
| `max_probe_overlap` | 0 | Max overlap (bp) allowed between designed probes. `0` disallows any overlap. `-1` disables the check entirely (probes may overlap freely). |
| `exclude_probes` | [] | List of [lhs, rhs] sequence pairs to exclude |

### Scoring Penalties

The pipeline uses a heuristic scoring system where lower scores indicate better probes. Penalties can be tuned:

- `suboptimal_gap_penalty` - Non-ideal junction chemistry
- `invalid_gap_penalty` - Poor junction chemistry
- `unbalanced_gc_penalty` - GC content difference between LHS/RHS
- `homopolymer_penalty` - Long homopolymer runs
- `lhs_tm_penalty` / `rhs_tm_penalty` - Melting temperature deviation
- `flex_overlap_penalty` - Overlap with existing 10x Flex probes
- `tandem_repeat_penalty` - Tandem repeat sequences

## Output

The output file contains designed probes with IDT-compatible ordering sequences. Each probe consists of:

- **LHS probe**: Left-hand side probe with 5' adapter sequence
- **RHS probe**: Right-hand side probe with 5' phosphate, barcode bridge, and 3' adapter

- There is also an expected gap sequence extracted and included, though this is not ordered from IDT as GIFT-seq will detect endogenous RNA sequence for the gap.

Probes are ordered 5' to 3' and include all necessary adapter sequences for the selected technology.
