# MPN 16-plex Probes

This directory contains inputs and outputs for MPN (Myeloproliferative Neoplasm) variant probes designed for 16-plex multiplexing across 8 barcodes.

## Reference Outputs

The `BC*.tsv` files contain the original probe designs:
- `BC1.tsv` through `BC8.tsv` - Designed probes for each barcode

Each file includes probes from multiple sources:
- `make_MPN_probes` - Pipeline-generated probes
- `competition` - Competition probes (JAK2 c.1849G>T, ASXL1 c.1934dupG)
- `splice_manual` - Manually designed splice site probes
- `manual` - Other manually designed probes

## Input Files

The `inputs_BC*.csv` files contain the HGVSc identifiers extracted from the pipeline-generated probes. Each barcode has slightly different variant coverage (~84-87 variants per barcode).

## Generating Probes

To regenerate probes for each barcode:

```bash
# Barcode 1
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 1 probe_sets/MPN_16plex/inputs_BC1.csv MPN_BC1_probes

# Barcode 2
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 2 probe_sets/MPN_16plex/inputs_BC2.csv MPN_BC2_probes

# Barcode 3
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 3 probe_sets/MPN_16plex/inputs_BC3.csv MPN_BC3_probes

# Barcode 4
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 4 probe_sets/MPN_16plex/inputs_BC4.csv MPN_BC4_probes

# Barcode 5
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 5 probe_sets/MPN_16plex/inputs_BC5.csv MPN_BC5_probes

# Barcode 6
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 6 probe_sets/MPN_16plex/inputs_BC6.csv MPN_BC6_probes

# Barcode 7
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 7 probe_sets/MPN_16plex/inputs_BC7.csv MPN_BC7_probes

# Barcode 8
python src/app.py --release 109 --skip_errors --config_file probe_sets/MPN_16plex/config.json --barcodes 8 probe_sets/MPN_16plex/inputs_BC8.csv MPN_BC8_probes
```

**Note**: The original probes were generated using custom transcript FASTA sequences (via `gene2fasta.pkl`) rather than Ensembl transcripts. To reproduce exactly, provide the transcript sequences in the input file (column 2 or 3).

### Command Options

| Option | Description |
|--------|-------------|
| `--skip_errors` | Skip variants that fail (recommended) |
| `--barcodes N` | Specify which barcode(s) to use for multiplexing |
| `--config_file` | Use the MPN-specific configuration |

## Configuration

The `config.json` contains scoring parameters optimized for MPN probe design:
- Bridge length: 2-10 bp
- GC content: 0.2-0.8 (lenient, with penalty)
- Excluded probes: Known problematic sequences

## Variant Coverage

Each barcode targets ~84-87 variants across 27 genes commonly mutated in MPN:

**Core Genes**: ASXL1, CALR, CBL, DNMT3A, EZH2, GATA2, HRAS, IDH1, IDH2, JAK2, KRAS, MPL, NFE2, NRAS, PPM1D, PTPN11, RRAS, RUNX1, SETBP1, SETD2, SF3B1, SH2B3, SRSF2, TET2, TP53, U2AF1, ZRSR2

Each gene has a wildtype (0bp) probe plus mutation-specific probes.

## Special Cases

- **ASXL1 c.2081_2099dup_0**: Has both a 0bp wildtype and a mutant probe (ASXL1 c.2081_2099dup_1)
- **Competition probes**: JAK2 c.1849G>T and ASXL1 c.1934dupG were designed externally
- **Splice site probes**: NFE2, TET2, and TP53 splice variants were manually designed