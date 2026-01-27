# Dual Probe IMPACT SNV Probes

This directory contains inputs and outputs for 0bp genotyping probes designed for dual-barcode experiments where the variant is positioned on either the LHS (left-hand side) or RHS (right-hand side) of the probe.

## Reference Outputs

| File | Description |
|------|-------------|
| `cl_dual_lhs_impact_snv_probes.tsv` | LHS genotyping probes (variant on left side) |
| `cl_dual_rhs_impact_snv_probes.tsv` | RHS genotyping probes (variant on right side) |

### Probe Sources

| Source | LHS Count | RHS Count | Description |
|--------|-----------|-----------|-------------|
| `make_0bp_genotyping_probes_lhs` | 189 | - | Pipeline-generated LHS probes |
| `make_0bp_genotyping_probes_rhs` | - | 189 | Pipeline-generated RHS probes |
| `cl_impact_snv` | 68 | 68 | CL IMPACT panel wildtype probes |
| `manual` | 4 | 4 | Manually designed probes |

## Input Files

| File | Variants | Description |
|------|----------|-------------|
| `inputs.csv` | 96 | Target variants for genotyping probes |

## Configuration Files

| File | Description |
|------|-------------|
| `config_lhs.json` | LHS genotyping (variant on left probe side) |
| `config_rhs.json` | RHS genotyping (variant on right probe side) |

## Generating Probes

### LHS Genotyping Probes (Barcode 1)

```bash
python src/app.py --mane --skip_errors --barcodes 1 --config_file probe_sets/dual_probe_impact/config_lhs.json probe_sets/dual_probe_impact/inputs.csv dual_lhs_probes
```

### RHS Genotyping Probes (Barcode 2)

```bash
python src/app.py --mane --skip_errors --barcodes 2 --config_file probe_sets/dual_probe_impact/config_rhs.json probe_sets/dual_probe_impact/inputs.csv dual_rhs_probes
```

### Command Options

| Option | Description |
|--------|-------------|
| `--mane` | Use MANE Select canonical transcripts (recommended) |
| `--skip_errors` | Skip variants that fail (recommended) |
| `--barcodes N` | Barcode assignment (1 for LHS, 2 for RHS) |
| `--config_file` | Use the appropriate genotyping configuration |

## Configuration Details

Both configurations share the same parameters except for `variant_side`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `flex_0bp_genotyping` | true | Enable 0bp genotyping mode |
| `variant_side` | "lhs" or "rhs" | Which probe side contains the variant |
| `max_bridge_length` | 0 | No gap between probes |
| `min_bridge_length` | 0 | No gap between probes |
| `last_base_penalty` | 0 | No penalty for last base |
| `first_base_penalty` | 0 | No penalty for first base |
| `distance_from_max_bridge_length_penalty` | 0 | No bridge distance penalty |

## 0bp Genotyping Mode

In 0bp genotyping mode:
- The variant is placed directly on either the LHS or RHS probe (no gap)
- LHS genotyping: Mutation is at the 3' end of the left probe
- RHS genotyping: Mutation is at the 5' end of the right probe
- Dual-barcode design allows both probe orientations to be used together

## Variant Coverage

The panel targets 96 variants across genes commonly found in IMPACT panels, including:

**Key Genes**: APC, AR, BRAF, CDKN2A, CTNNB1, EGFR, ERBB2, FBXW7, FGFR3, GNAS, IDH2, KRAS, MYD88, NRAS, PIK3CA, PTEN, SF3B1, SMAD4, TP53, U2AF1, XPO1

Special variants include synonymous SNPs (`c.XXX>*`) for wildtype detection.

## Special Cases

- **Wildtype probes**: KRAS 0bp, MYD88 0bp, FBXW7 0bp
- **CL IMPACT panel probes**: 68 wildtype probes from external panel
- **Manual probes**: JAK2 c.1849G>T, BCR-ABL c.fusion
