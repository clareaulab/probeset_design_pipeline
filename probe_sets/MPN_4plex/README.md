# MPN 4-plex Probes

This directory contains inputs and outputs for MPN (Myeloproliferative Neoplasm) variant probes designed for 4-plex experiments (single barcode).

## Reference Output

The `MPN_4plex_probes.tsv` file contains the original probe designs with 95 probes across multiple sources:
- `MPN_v2_probes` - Pipeline-generated probes
- `Patient_specific_panel_4plex` - Patient-specific panel probes
- `competition` - Competition probes (JAK2 c.1849G>T)
- `manual` - Manually designed probes (CALR c.1099_1150del)
- `MPN_additional_probe_23bp` - Additional probe with modified length

## Input File

The `inputs.csv` file contains all 95 HGVSc identifiers from the reference output.

## Generating Probes

To regenerate probes:

```bash
python src/app.py --mane --skip_errors --config_file probe_sets/MPN_4plex/config.json probe_sets/MPN_4plex/inputs.csv MPN_4plex_probes
```

**Note**: The original probes were generated using custom transcript FASTA sequences (via `gene2fasta.pkl`) rather than Ensembl transcripts. To reproduce exactly, provide the transcript sequences in the input file (column 2 or 3).

### Command Options

| Option | Description |
|--------|-------------|
| `--mane` | Use MANE Select canonical transcripts (recommended) |
| `--skip_errors` | Skip variants that fail (recommended) |
| `--config_file` | Use the MPN 4-plex specific configuration |

## Configuration

The `config.json` contains scoring parameters optimized for MPN 4-plex probe design. Key differences from the 16-plex configuration:

| Parameter | 4-plex | 16-plex |
|-----------|--------|---------|
| `suboptimal_gap_penalty` | 10 | 11 |
| `invalid_gap_penalty` | 1e8 | 16 |
| `distance_from_max_bridge_length_penalty` | 8 | 5 |
| `unbalanced_gc_penalty` | 1 | 2 |
| `homopolymer_penalty` | 2 | 0 |
| `lhs_tm_penalty` | 1 | 2 |
| `rhs_tm_penalty` | -1 | -2 |

## Variant Coverage

The panel targets 95 variants across 26 genes commonly mutated in MPN:

**Genes**: ASXL1, CALR, CBL, CHEK2, CUX1, DNMT3A, ETV6, EZH2, IDH1, IDH2, JAK2, KIT, KRAS, MPL, NFE2, NRAS, PPM1D, PTPN11, RUNX1, SETBP1, SF3B1, SH2B3, SRSF2, TET2, TP53, U2AF1

Each gene has a wildtype (0bp) probe plus mutation-specific probes.

## Special Cases

- **JAK2 c.1849G>T**: Competition probe (externally designed)
- **CALR c.1099_1150del**: Manually designed for large deletion
- **EZH2 c.2214_2218dup**: Additional probe with 23bp LHS (non-standard length)
