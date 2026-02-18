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

### 1. MSK IMPACT SNVs

```bash
python src/app.py --release 109 --skip_errors --config_file probe_sets/impact/config.json probe_sets/impact/inputs.csv impact_probes
```

**Note:** Uses `--msk` flag for MSK IMPACT canonical transcript overrides.


### 2. JAK2 and BCR-ABL Fusions (2 probes)

These were generated manually

## Configuration Files

- `config.json` - Standard config (max_bridge_length=10, strict_gc_content=true)

## Input Files

- `inputs.csv` - MSK IMPACT panel variants
