# NK interferon gamma context analysis

Reproducible secondary analyses of published tumor–natural killer cell experiments. The code compares net IFNγ responses with matched functional perturbations, checks a fixed ligand-prediction benchmark, and describes genetic-screen and host contexts.

## Quick start

This public repository contains the analysis code and methods documentation. The numerical inputs are supplied separately in `Analysis_Code_and_Data.zip`, the anonymized archive accompanying the manuscript submission. That archive includes the scripts shown here, verified numerical inputs, frozen predictions, source coordinates, genetic-screen results and figure input tables. Extract it into a new directory and run the following commands from that directory:

```text
python -m pip install -r requirements.txt
python scripts/verify_frozen_core.py
python scripts/verify_inference.py
python scripts/verify_extension_tables.py
python scripts/make_figures.py
python scripts/make_extension_figures.py
```

The first three commands independently check the numerical claims. The last two regenerate all seven figures and require Arial to reproduce the original typography. The archive already contains the reference figure PNGs. Python 3.13.5, NumPy 2.4.3, pandas 2.3.3, SciPy 1.16.3, matplotlib 3.10.0, openpyxl 3.1.5, Pillow 11.1.0 and lxml 5.3.0 were used for the scientific validation.

## Contents of the accompanying analysis archive

- `source_data/`: extracted functional observations, published donor source workbooks, matched molecular inputs, model specifications, all frozen prediction comparisons and sensitivity outputs used by the retained analysis.
- `analysis/`: calibrated source coordinates, reconstructed functional comparisons, 222,609 genetic-screen rows, the prespecified 21-gene display and 28 mouse observations.
- `provenance/`: source-study identities, hashes, screen designs and ambiguous identifiers.
- `scripts/`: independent recalculation, source extraction and figure generation.
- `methods_history/`: original algorithm and acquisition scripts with their historical relative-layout assumptions; all frozen sensitivity outputs remain in `provenance/phase1_analysis_archive/`.
- `tables/` and `figures/`: corresponding numerical display tables and reference figures.

The ordinary workflow starts with the included verified extracted inputs. `reconstruct_new_functional.py`, `extract_host_context.py` and `process_screens.py` additionally support reconstruction from the original publisher figures/workbooks. Those files must be retrieved from the sources and placed at the paths in `provenance/ORIGINAL_SOURCE_FILES.csv` before these three extraction scripts are run. Full articles and the 42 MB original genetic-screen workbook are not bundled. The derived full genetic-screen table is included. The large original DepMap expression release is not duplicated; exact matched model inputs are retained.

## Interpretation

These are secondary analyses of previously published experiments. Image-reading bounds are not biological confidence intervals. Shared controls, repeated donors, related tumor models and unpaired source observations retain their source-specific dependence. Cytokine addition, genetic perturbation, baseline cytotoxicity and multiday selection are distinct outcomes. Negative prediction comparisons and all prespecified alternatives are retained.

The archive reproduces the retained manuscript analysis. Later exploratory strengthening analyses were not adopted into the manuscript and are not presented as independent validation here. Original publications remain the source of their experiments and biological mechanisms; see `provenance/STUDY_SOURCES.md`.
