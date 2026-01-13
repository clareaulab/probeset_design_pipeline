# GBM Probes

This directly allows you to generate probes for GBM except for one manually designed splice variant probe.

To generate the probes:

```bash
python src/app.py --config_file probe_sets/GBM/config.json probe_sets/GBM/inputs.csv GBM_probes
```