# WMH analysis code structure

This repository keeps the original notebook intact and also exposes each executable notebook cell as an ordered Python script under `scripts/notebook_steps/`. The exported scripts preserve the notebook order and do not remove analysis functionality.

## Directory layout

```text
.
├── WMH_Analysis.ipynb              # Original end-to-end analysis notebook
├── README.md                       # Project overview and run instructions
├── requirements.txt                # Python package list used by the analysis
├── run_apoe.sh                     # APOE extraction script for UKB RAP
├── docs/
│   └── analysis_pipeline.md        # This structure guide
└── scripts/
    ├── export_notebook_steps.py    # Regenerate scripts from WMH_Analysis.ipynb
    ├── run_notebook_steps.sh       # Sequential runner for exported steps
    └── notebook_steps/             # Ordered Python exports of notebook cells
```


## Beginner FAQ: where are the changes?

The changes are regular repository files, not separate attachments. Once this branch or pull request is checked out, you will see the added files directly in the project tree: `scripts/export_notebook_steps.py`, `scripts/run_notebook_steps.sh`, `scripts/notebook_steps/`, `requirements.txt`, and this documentation file.

You do not need to download code manually from the chat. If you are using a local clone, update it with your normal Git workflow, for example by pulling or merging the branch that contains these commits.

## Exported analysis steps

| Step | Script | Purpose |
|---:|---|---|
| 01 | `scripts/notebook_steps/01_flowchart.py` | Build the study flowchart outputs. |
| 02 | `scripts/notebook_steps/02_remove_withdrawn_consent.py` | Remove participants who withdrew consent. |
| 03 | `scripts/notebook_steps/03_clean_final_dataset.py` | Clean variables, derive analysis fields, and write final processed data. |
| 04 | `scripts/notebook_steps/04_matching_ato_weighting.py` | Create PS-matched and ATO-weighted cohorts with balance outputs. |
| 05 | `scripts/notebook_steps/05_wmh_transformation_diagnostics.py` | Compare raw vs. log1p head-size-normalized WMH models. |
| 06 | `scripts/notebook_steps/06_main_sensitivity_psm_wmh.py` | Run primary and sensitivity matched-cohort WMH models and tables. |
| 07 | `scripts/notebook_steps/07_cognitive_secondary_outcomes.py` | Run cognitive secondary outcome analyses. |
| 08 | `scripts/notebook_steps/08_wmh_cognitive_mediation.py` | Run SA → WMH → cognition mediation analyses. |
| 09 | `scripts/notebook_steps/09_mediation_plot_addon.py` | Generate mediation summary plots. |
| 10 | `scripts/notebook_steps/10_diagnostics_tables_5_cohorts.py` | Generate diagnostic tables across five cohorts, including ATO support. |
| 11 | `scripts/notebook_steps/11_additional_sensitivity_analysis.py` | Run additional sensitivity analyses. |
| 12 | `scripts/notebook_steps/12_stratify_cohorts.py` | Build stratified primary and sensitivity cohort files. |
| 13 | `scripts/notebook_steps/13_stratified_ols_fdr.py` | Run stratified OLS analyses with in-stratum BH-FDR. |
| 14 | `scripts/notebook_steps/14_interaction_ols_fdr.py` | Run interaction OLS analyses with BH-FDR on interaction tests. |
| 15 | `scripts/notebook_steps/15_annual_wmh_change.py` | Run annual WMH change analyses. |
| 16 | `scripts/notebook_steps/16_duration_vs_wmh_loess.py` | Run diagnosis-duration vs. WMH LOESS and spline analyses. |
| 17 | `scripts/notebook_steps/17_print_package_versions.py` | Print package versions for reproducibility. |

## Running the structured scripts

Regenerate the exported scripts after editing `WMH_Analysis.ipynb`:

```bash
python scripts/export_notebook_steps.py
```

Run all exported steps from the repository root:

```bash
bash scripts/run_notebook_steps.sh
```

Run a selected contiguous range:

```bash
START_STEP=3 END_STEP=6 bash scripts/run_notebook_steps.sh
```

The scripts expect the same input files and working directory assumptions as the original notebook. Edit the CONFIG sections in the relevant script or notebook cell when your local input filenames differ.


## What changed versus the original notebook-only layout

The original notebook remains the authoritative interactive workflow so no analysis functionality is removed. The practical code-structure change is that the same code is now available as ordered Python files, with a regeneration command for keeping those files synchronized with the notebook. If you edit the notebook, run `python scripts/export_notebook_steps.py`; if you edit an exported script directly, make the same change in the notebook before regenerating so the two views do not diverge.
