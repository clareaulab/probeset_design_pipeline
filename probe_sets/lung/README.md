# Lung Probes

This directory contains inputs and outputs for lung cancer variant probes designed for single-barcode experiments.

## Reference Output

The `lung_probes.tsv` file contains 384 probes from multiple sources:
- `make_lung_met_probes` (378) - Pipeline-generated standard probes (25bp LHS)
- `CL_impact_panel` (4) - CL IMPACT panel probes
- `make_short_lung_met_probe` (2) - Short probes (23bp LHS) for XRCC2

## Input Files

| File | Probes | Description |
|------|--------|-------------|
| `inputs.csv` | 378 | Standard probe variants (25bp LHS) |
| `inputs_short.csv` | 2 | Short probe variants (23bp LHS) for XRCC2 |

## Configuration Files

| File | Description |
|------|-------------|
| `config.json` | Standard probes (25bp LHS, 25bp RHS) |
| `config_short.json` | Short probes (23bp LHS, 25bp RHS) with extra excluded probe |

## Generating Probes

### Standard Probes (25bp)

```bash
python src/app.py --mane --skip_errors --config_file probe_sets/lung/config.json probe_sets/lung/inputs.csv lung_probes
```

### Short Probes (23bp LHS)

```bash
python src/app.py --mane --skip_errors --config_file probe_sets/lung/config_short.json probe_sets/lung/inputs_short.csv lung_probes_short
```

**Note**: The original probes were generated using custom transcript FASTA sequences rather than Ensembl transcripts. To reproduce exactly, provide the transcript sequences in the input file.

### Command Options

| Option | Description |
|--------|-------------|
| `--mane` | Use MANE Select canonical transcripts (recommended) |
| `--skip_errors` | Skip variants that fail (recommended) |
| `--config_file` | Use the lung-specific configuration |

## Configuration Details

Both configurations share most parameters (same as MPN 4-plex):
- Bridge length: 2-10 bp
- GC content: 0.2-0.8 (lenient, with penalty)
- `invalid_gap_penalty`: 1e8 (strict gap chemistry)
- `homopolymer_penalty`: 2

**Short probe differences** (`config_short.json`):
- `lhs_probe_length`: 23 (vs 25)
- Extra excluded probe pair

## Variant Coverage

The panel targets 378+ variants across many genes commonly mutated in lung cancer, including:

**Key Genes**: AGO2, AMER1, ARID2, ATM, ATR, BAP1, BCL2, BRAF, CCNE1, CHEK1, DICER1, DROSHA, EGFR (implied), EPHA3, EPHA5, ERBB4, ETV1, FANCC, FGF3, FGF4, FGFR4, FH, HGF, INPP4A, INHBA, JAK3, KEAP1, KLF4, KRAS, MDC1, MDM4, MYCN, MYD88, NCOA3, NCOR1, NF1, NFE2L2, NOTCH1, NOTCH2, NSD1, NTHL1, NTRK1, NTRK2, NTRK3, NUP93, PAK5, PARP1, PHOX2B, PIK3C2G, PIK3CG, POLE, PREX2, PRKD1, PRKN, PTPRD, PTPRS, RAC1, RB1, RPTOR, RRAS, SCG5, SHOC2, SLX4, SMARCA4, SOX9, SOS1, SPOP, STAG2, STK11, TGFBR2, TET1, TP53, XRCC2, YAP1

Each gene typically has a wildtype (0bp) probe plus mutation-specific probes.

## Special Cases

- **XRCC2 probes**: Use short configuration (23bp LHS) - `inputs_short.csv`
- **CL IMPACT panel probes**: 4 probes from external panel (PTEN variants)
