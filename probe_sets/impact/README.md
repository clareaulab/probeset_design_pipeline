# MSK IMPACT Probe Set

This probe set contains 163 probes from multiple sources for detecting cancer-associated mutations.

## Probe Sources

| Source | Count | Description |
|--------|-------|-------------|
| impact | 49 | MSK IMPACT panel SNVs |
| snv | 45 | Cell line-specific SNVs (K562, MM1S, SET2, SUPT1) |
| 0bp | 68 | Wildtype (0bp) probes for gene detection |
| competition | 1 | JAK2 c.1849G>T with non-standard probe lengths |

## Recreating the Probes

### 1. MSK IMPACT SNVs (49 probes)

```bash
python src/app.py --msk --skip_errors --config_file probe_sets/impact/config.json probe_sets/impact/inputs_impact.csv impact_probes
```

**Note:** Uses `--msk` flag for MSK IMPACT canonical transcript overrides.

### 2. Cell Line SNVs (45 probes)

```bash
python src/app.py --mane --skip_errors --config_file probe_sets/impact/config.json probe_sets/impact/inputs_snv.csv snv_probes
```

### 3. Wildtype 0bp Probes (68 probes)

```bash
python src/app.py --mane --skip_errors --config_file probe_sets/impact/config_0bp.json probe_sets/impact/inputs_0bp.csv zerobp_probes
```

### 4. BCR-ABL Fusion Probe (1 probe)

```bash
python src/app.py --skip_errors --config_file probe_sets/impact/config.json probe_sets/impact/inputs_bcr_abl.csv bcr_abl_probes
```

**Note:** BCR-ABL uses a custom fusion sequence provided in the input file.

### 5. JAK2 Competition Probe (1 probe)

```bash
python src/app.py --msk --skip_errors --config_file probe_sets/impact/config_jak2.json probe_sets/impact/inputs_jak2.csv jak2_probes
```

**Note:** JAK2 c.1849G>T uses non-standard probe lengths (28bp LHS, 44bp RHS).

## Configuration Files

- `config.json` - Standard config (max_bridge_length=10, strict_gc_content=true)
- `config_0bp.json` - 0bp probes config (max/min_bridge_length=0)
- `config_jak2.json` - JAK2 special config (lhs=28bp, rhs=44bp)

## Input Files

- `inputs_impact.csv` - MSK IMPACT panel variants
- `inputs_snv.csv` - Cell line SNVs from WIP_cell_mixing_probes
- `inputs_0bp.csv` - Gene names for wildtype probes
- `inputs_bcr_abl.csv` - BCR-ABL fusion with custom sequence
- `inputs_jak2.csv` - JAK2 competition variant

## Known Issues (from Probe_Notes.txt)

- Outdated HGNC names were replaced with current names
- RPS6KA4 c.2308del - Missing due to not enough room at the RHS
- STK19 c.265G>A - Renamed to c.264G>A due to transcript update
- SMAD2 c.1391C>G - Missing due to not enough room at the RHS
- CASP8 c.1596delA - Missing due to not enough room at the RHS
- MYD88 c.794T>C - Can't find the T in the reference transcript

## Cell Line SNV Sources

The SNV probes originated from these cell lines (from WIP_cell_mixing_probes/):
- K562_heterozygous (7 variants)
- K562_homozygous (5 variants)
- K562_RNAedits (1 variant)
- MM1S_heterozygous (8 variants)
- MM1S_homozygous (6 variants)
- SET2_heterozygous (3 variants)
- SET2_homozygous (5 variants)
- SUPT1_heterozygous (3 variants)
- SUPT1_homozygous (7 variants)
