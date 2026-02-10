# GBM Probes

This directory contains inputs to generate probes for GBM variants.

## Generating Probes

```bash
python src/app.py --release 75 --skip_errors --config_file probe_sets/GBM/config.json probe_sets/GBM/inputs.csv GBM_probe
```

Note that EGFRvIII probes were created manually and the PTEN c.37_41del where generated, but not used in any experiments.