# HeLa Pre/Post LC-MS QC Workflow

This repository contains an R-based quality-control workflow for comparing two theoretically identical HeLa samples acquired at different timepoints:
- preRun: acquired before the broader LC-MS/MS sequence
- postRun: acquired after the broader LC-MS/MS sequence

The goal is to quantify whether data quality declined over run order and provide outputs that are both machine-readable (CSV) and visually interpretable (plots + HTML report).

## What This Workflow Does

The script ingests paired pre/post text exports and computes:
- run-level summary deltas (counts, abundance summaries, MS/MS-level metrics)
- shared-feature intensity drift (peptide-level and protein-group level)
- run-order bias indicators (direction and magnitude of change)
- QC visualization panels for fast manual review

Outputs include:
- consolidated metric summary CSV
- detailed shared-feature tables
- QC plots
- self-contained HTML summary report

## Input Requirements

The input directory must contain paired pre/post files with these file types:
- MSMSSpectrumInfo
- PSMs
- PeptideGroups
- ProteinGroups
- ResultStatistics
- StatisticsInsights

Expected naming pattern includes:
- *_preRun_*
- *_postRun_*

The script validates that exactly one preRun and one postRun file exist for each required file type.

## Quick Start

From the repo root:

```bash
Rscript qc_HeLa.R \
  --input-dir='/absolute/path/to/data'
```

Optional arguments:

```bash
Rscript qc_HeLa.R \
  --input-dir='/absolute/path/to/data' \
  --output-dir='/absolute/path/to/output' \
  --render-report=true \
  --write-preprocess-script=true
```

## Output Location Behavior

If --output-dir is not provided, output is written to:
- dirname(input_dir)/output

This keeps output beside the data folder by default.

## Generated Files

Main output files:
- qc_metrics_summary.csv
- qc_shared_peptides.csv
- qc_shared_protein_groups.csv
- qc_overlap_summary.csv
- qc_peptide_bias_summary.csv
- qc_psm_rt_bias.csv
- qc_report.html
- preprocess_tables.sh

Plot files (output/plots):
- 01_high_confidence_counts.png
- 02_shared_peptide_scatter.png
- 03_ma_plot.png
- 04_log2fc_histogram.png
- 05_rt_binned_log2fc_trend.png
- 06_psm_precursor_abundance_boxplot.png
- 07_isolation_interference_violin.png
- 08_psm_rt_bias_bar.png

## How To Interpret Key Outputs

### 1) qc_metrics_summary.csv
Columns:
- metric: normalized metric identifier
- preRun, postRun: run values
- delta_abs: postRun - preRun
- delta_pct: percent change relative to preRun
- decline_flag:
  - major_decline: <= -20%
  - minor_decline: <= -5%
  - stable: between -5% and +5%
  - increase: > +5%

### 2) qc_peptide_bias_summary.csv
Summarizes shared peptide abundance drift:
- median_log2_fc: median log2(post/pre)
- frac_negative_fc: fraction of shared peptides with lower postRun abundance
- frac_abs_fc_gt1: fraction with |log2 FC| > 1
- wilcox_p: Wilcoxon signed-rank p-value vs 0 shift

### 3) qc_report.html
Provides:
- run-order bias summary
- Top metric changes table sorted by event severity first:
  1. major_decline
  2. minor_decline
  3. stable
  4. increase
- embedded QC plots

## Understanding The R Script

This section is intended for users maintaining or extending qc_HeLa.R.

### Script Layout

The script is organized in functional blocks:
1. CLI/config parsing
2. Input discovery and validation
3. Data loading helpers
4. Metric extraction
5. Shared-feature analyses
6. Bias quantification
7. Plot generation
8. HTML report generation
9. Main orchestration

### Function Groups

Configuration and utilities:
- parse_cli_args: reads command-line options and applies defaults
- stop_if_missing_columns: defensive schema checks
- sanitize_metric_name / normalize_sample_metric: harmonizes metric naming

Input and loading:
- discover_input_files: validates required pre/post pairs
- load_tables / get_table: loads and retrieves typed tables

Quantification:
- extract_statistics_insights
- extract_result_statistics
- build_summary_metrics
- analyze_shared_peptides
- analyze_protein_groups
- build_psm_quality_table
- build_msms_quality_table
- compute_bias_metrics

Output:
- save_plots
- write_html_report
- write_preprocess_bash_script
- main

### Syntax and Style Notes

The script uses modern tidyverse conventions:
- native R pipe: |>
- column-safe data access: .data[[col_name]]
- grouped summaries: group_by + summarise
- reshaping: pivot_wider
- pattern matching: stringr::str_detect

Error handling patterns:
- explicit stop(...) for missing paths, files, or required schema
- tryCatch(...) for optional statistical tests (avoids hard failure)

### How Shared Features Are Defined

Peptide-level matching uses:
- Annotated Sequence
- Modifications

Key format:
- Annotated Sequence || Modifications

This is intentionally conservative to reduce ambiguous mapping.

### Why Isolation Interference Uses MSMSSpectrumInfo

In many exports, PSM tables may not include complete isolation interference columns.
The workflow now reads isolation interference from MSMSSpectrumInfo so the violin plot is populated reliably.

## Performance and Preprocessing

The generated preprocess helper script:
- output/preprocess_tables.sh

It optionally gzips large text files (>10M) while keeping originals.
This can reduce storage pressure for archival workflows.

## Common Troubleshooting

1) No files detected
- Confirm input path points to the data folder containing .txt files.
- Confirm filenames include preRun/postRun markers.

2) Missing required file pair(s)
- Ensure each of the 6 required file types exists for both preRun and postRun.

3) Empty or sparse plots
- Check source tables for corresponding columns.
- Verify export completeness from the upstream proteomics software.

4) Report not generated
- Ensure --render-report is true.
- Check write permissions for the output directory.

## Minimal Reproducible Command

```bash
Rscript qc_HeLa.R --input-dir='/home/jkg/OneDrive/Documents/All Files - personal/School_UT/Work/SG-01-008_BlackSoldierFlyGutMicrobiome/qc_hela/data'
```

## Repository Contents

- qc_HeLa.R: main workflow script
- README.md: documentation and usage guide
- output/: generated results (if run in-repo)
