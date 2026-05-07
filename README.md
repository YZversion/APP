# HD RNA-seq Explorer

**An integrated R Shiny application for interactive exploration of Huntington's Disease transcriptomic data.**

---

## Overview

HD RNA-seq Explorer is a single, self-contained R Shiny application built for BF591 (Bioinformatics with R) at Boston University. The app provides four fully integrated interactive modules for exploring a publicly available RNA-seq dataset comparing post-mortem brain tissue from Huntington's Disease (HD) patients and neurologically normal controls.

The application allows users to:

- Explore sample-level clinical and technical metadata
- Filter and visualize a normalized gene expression counts matrix
- Examine pre-computed differential expression results through an interactive table and volcano plot
- Investigate the expression of any individual gene across sample groups

All analysis is performed reactively in-browser — no command-line R or pre-generated plots are required to use the app.

---

## Dataset

| Property | Value |
|---|---|
| GEO Accession | [GSE64810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810) |
| SRA Accession | SRP051844 |
| Tissue | Post-mortem human prefrontal cortex (Brodmann Area 9) |
| Sequencing Platform | Illumina HiSeq 2000 |
| Total Samples | 69 (20 HD, 49 neurologically normal controls) |
| Comparison | HD vs. Control (HD is the numerator for log₂FC) |
| Normalization | DESeq2 size-factor normalization (pre-computed by study authors) |

### Sample Metadata

Key clinical and technical variables per sample:

| Column | Type | Description |
|---|---|---|
| `sample_id` | character | Study identifier (C_#### or H_####) |
| `geo_accession` | character | GEO sample accession (GSM...) |
| `diagnosis` | categorical | `Control` or `HD` |
| `age_at_death` | numeric | Age at death in years |
| `pmi_hours` | numeric | Post-mortem interval in hours |
| `rin` | numeric | RNA integrity number (0–10) |
| `mrna_seq_reads` | numeric | mRNA-seq read count |
| `vonsattel_grade` | numeric | HD neuropathology grade (HD samples only) |
| `cag_repeat` | numeric | HTT CAG repeat length (HD samples only) |
| `hv_striatal_score` | numeric | Hadzi–Vonsattel striatal score |
| `hv_cortical_score` | numeric | Hadzi–Vonsattel cortical score |
| `sex` | — | Not available in the GEO series matrix extract (all NA) |

> **Note:** `vonsattel_grade`, `cag_repeat`, `age_of_onset`, and `disease_duration_years` are NA for all Control samples — these are HD-specific covariates.

### Processed Files

The preprocessing pipeline produces five app-ready CSV files in `data/processed/`:

| File | Approx. Size | Description |
|---|---|---|
| `sample_info.csv` | < 1 MB | 69 rows × 19 columns of sample metadata |
| `counts_wide.csv` | ~30 MB | Genes × samples normalized counts matrix |
| `differential_expression.csv` | ~6 MB | Pre-computed DESeq2 DE results (HD vs. Control) |
| `gene_expression_long.csv` | ~147 MB | Long-format pivot of counts × metadata (pipeline output only) |
| `app_data_dictionary.csv` | < 1 MB | Column-level documentation for all files |

> **Important:** `gene_expression_long.csv` is **not loaded** by the Shiny app. The app pivots reactively from `counts_wide.csv` + `sample_info.csv` to avoid the 147 MB overhead.

> **TODO:** Confirm full dataset citation. Based on GEO accession, likely: Labadorf A et al. (2015) *RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression.* PLOS ONE 10(12): e0143563. https://doi.org/10.1371/journal.pone.0143563 — verify before submission.

---

## App Modules

### Tab 1 — Sample Information Exploration

**Purpose:** Understand the clinical and technical composition of the 69-sample cohort before interpreting any expression results.

**Input:** `sample_info.csv`

**Sub-tabs and controls:**

| Sub-tab | Output | Description |
|---|---|---|
| Summary | Text + table | Row/column counts; per-column type, mean ± SD (numeric) or up to 5 distinct values (categorical) |
| Data Table | Interactive `DT` table | Full sample metadata, sortable and filterable by any column |
| Distributions | Violin or histogram | Choose any numeric column; optionally group by any categorical column — renders a grouped violin + jitter when grouped, a histogram when ungrouped |

**Interpretation note:** Use this tab to check batch variables (PMI, RIN) for potential confounding before drawing biological conclusions from the expression tabs.

---

### Tab 2 — Counts Matrix Exploration

**Purpose:** Assess data quality and identify informative genes before differential expression interpretation.

**Input:** `counts_wide.csv`

**Controls:**

| Control | Default | Description |
|---|---|---|
| Variance percentile slider | 50% | Retains genes at or above the Nth variance percentile |
| Minimum non-zero samples slider | 10 | Retains genes expressed in at least N samples |
| Log₂-transform toggle | On | Applies log₂(x + 1) to the heatmap matrix |
| PCA display type | Scatter | Switch between PC-vs-PC scatter and beeswarm view of top N PCs |

**Sub-tabs:**

| Sub-tab | Output | Description |
|---|---|---|
| Filter Summary | Summary table | Total genes, passing genes (count + %), failing genes (count + %), sample count |
| Diagnostic Plots | Two scatter plots | (1) log₁₀ median vs. log₁₀ variance; (2) log₁₀ median vs. zero-count samples — passing genes colored darker, filtered genes lighter |
| Heatmap | `pheatmap` | Hierarchically clustered heatmap of filtered genes (capped at 500 highest-variance genes); RdBu color scale |
| PCA | Scatter or beeswarm | `prcomp` on the transposed filtered matrix (samples as rows); axis labels include % variance explained |

**Interpretation note:** The diagnostic scatter plots reveal the relationship between expression level and variability. Low-expression, high-zero genes filtered out here are typically noise rather than biology.

---

### Tab 3 — Differential Expression

**Purpose:** Identify genes significantly dysregulated in HD relative to controls, and explore DE results interactively.

**Input:** `differential_expression.csv` (pre-computed DESeq2 results, HD vs. Control)

**Controls:**

| Control | Description |
|---|---|
| X-axis column | Choose from: `log2FoldChange`, `baseMean`, `lfcSE`, `stat`, `pvalue`, `padj` |
| Y-axis column | Choose from: `neg_log10_padj`, `pvalue`, `padj`, `stat`, `baseMean`, `log2FoldChange` |
| Base point color | Color picker for non-significant points (default: grey `#999999`) |
| Highlight color | Color picker for significant points (default: red `#d7191c`) |
| Highlight threshold slider | Points where −log₁₀(padj) exceeds this value are highlighted (range 0–300) |
| Gene symbol search | Case-insensitive partial match filter on the table |
| Update Volcano Plot button | Renders the volcano only on click, preventing layout thrash on slider interaction |

**Sub-tabs:**

| Sub-tab | Output | Description |
|---|---|---|
| Table | Interactive `DT` table | Filtered DE results, sortable by any column |
| Volcano Plot | ggplot2 scatter | Customizable volcano; highlighted points exceed the padj threshold |

**Interpretation note:** The comparison is HD vs. Control — positive log₂FC means higher expression in HD. The `neg_log10_padj` column is precomputed (−log₁₀(padj)) so it can be used directly as the y-axis.

---

### Tab 4 — Individual Gene Expression Visualization

**Purpose:** Examine the normalized expression of any single gene across all 69 samples, stratified by any clinical grouping variable.

**Inputs:** `counts_wide.csv` + `sample_info.csv` (both uploaded independently)

**Controls:**

| Control | Description |
|---|---|
| Gene selector | Searchable selectize populated from row names of the counts matrix |
| Group by | Dropdown of categorical columns from sample metadata (defaults to `diagnosis`) |
| Plot type | Boxplot, Violin, Bar (mean ± SE), or Beeswarm |
| Plot button | Renders the plot on click (`bindEvent`) |

**Outputs:**

| Output | Description |
|---|---|
| Expression plot | Selected plot type, colored by group |
| Gene table | Reactive `DT` table of normalized counts joined to sample metadata for the selected gene |

**Interpretation note:** Bar plots show group mean ± SE. Violin and beeswarm plots show the full distribution. For HD-specific columns such as `vonsattel_grade`, grouping by diagnosis first is recommended since control samples have NA for those columns.

---

## Repository Structure

```
.
├── app/
│   └── app.R                      # Complete Shiny application (~723 lines)
├── data/
│   ├── raw/                       # .gitignore'd — downloaded GEO source files
│   ├── metadata/
│   │   └── samples_clean.csv      # Committed — 69-row sample annotation
│   └── processed/                 # .gitignore'd — app-ready CSVs (generated by pipeline)
├── docs/
│   └── data_preparation_plan.md   # Workflow notes and CSV schemas
├── results/                       # .gitignore'd — analysis outputs
├── scripts/
│   ├── 01_download_data.R         # Downloads 3 GEO files into data/raw/
│   ├── 02_prepare_processed_csvs.R # Parses GEO files → 5 app-ready CSVs (~2 min)
│   ├── 03_validate_processed_csvs.R # Sanity-checks all processed files
│   ├── run_all_data_prep.R        # Runs all three steps in order
│   └── README.md                  # Scripts usage notes
├── tests/
│   ├── test_tab1.R                # 13 assertions — CSV loading, column summary
│   ├── test_tab2.R                # 30 assertions — matrix filtering, PCA, heatmap
│   ├── test_tab3.R                # 17 assertions — DE table filtering, volcano logic
│   └── test_tab4.R                # 21 assertions — gene extraction, pivoting, plotting
├── .gitignore
├── CLAUDE.md                      # Developer guidance (not for end users)
├── LICENSE
└── README.md
```

> `data/raw/` and `data/processed/` are excluded from version control. Only `data/metadata/samples_clean.csv` is committed. The processed files must be regenerated locally before the app will run.

---

## Installation and Requirements

### R Version

R ≥ 4.2.0 is recommended.

### CRAN Packages

Install all required packages with:

```r
install.packages(c(
  # Shiny app
  "shiny",
  "bslib",
  "ggplot2",
  "DT",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "ggbeeswarm",
  "colourpicker",
  # Preprocessing pipeline
  "dplyr",
  "readr",
  # Tests
  "testthat",
  "here"
))
```

### Bioconductor Packages

This app uses **pre-computed** DESeq2 results provided by the original study authors and downloaded directly from NCBI GEO. The Shiny application itself does not call any Bioconductor packages. No BiocManager installation is required to run the app.

> If you wish to re-run the differential expression analysis independently, you would need `BiocManager::install("DESeq2")`, but this is not part of the current workflow.

---

## How to Run the App

### Step 1 — Clone the repository

```bash
git clone https://github.com/YZversion/APP.git
cd APP
```

### Step 2 — Generate the processed data files

Run from the project root (RStudio terminal or PowerShell). This downloads three GEO files (~50 MB total) and generates the five processed CSVs:

```r
Rscript scripts/run_all_data_prep.R
```

Expected runtime: approximately 2 minutes. Individual steps:

```r
Rscript scripts/01_download_data.R         # Download GEO files
Rscript scripts/02_prepare_processed_csvs.R # Build app-ready CSVs
Rscript scripts/03_validate_processed_csvs.R # Verify outputs
```

### Step 3 — Install dependencies

```r
install.packages(c(
  "shiny", "bslib", "ggplot2", "DT", "tidyr",
  "pheatmap", "RColorBrewer", "ggbeeswarm", "colourpicker"
))
```

### Step 4 — Launch the app

```r
shiny::runApp("app/app.R")
```

Or from within the `app/` directory:

```r
setwd("app")
shiny::runApp()
```

The app will open in your default browser. Upload the files from `data/processed/` using the file inputs in each tab.

---

## Input File Formats

All inputs are CSV files. The app validates that each upload is a `.csv` with more than one row and more than one column, showing an inline error message otherwise.

### `sample_info.csv`

| Requirement | Detail |
|---|---|
| Required columns | `sample_id`, `diagnosis` |
| `sample_id` format | Character: `C_####` (control) or `H_####` (HD) |
| `diagnosis` values | Exactly `Control` or `HD` |
| Rows | One row per sample (69 total) |

### `counts_wide.csv`

| Requirement | Detail |
|---|---|
| Column 1 | `gene_id` — Ensembl ID (e.g., `ENSG00000000003.15`) |
| Column 2 | `gene_symbol` — HGNC gene symbol |
| Column 3 | `gene_name` — Gene description (may be empty) |
| Columns 4–72 | One numeric column per sample; column names must match `sample_id` values in `sample_info.csv` |
| Values | DESeq2 size-factor normalized counts (non-negative real numbers) |

The app strips columns 1–3 automatically before any analysis and sets `gene_symbol` as row names for gene lookup.

### `differential_expression.csv`

| Requirement | Detail |
|---|---|
| Required columns | `gene_id`, `gene_symbol`, `log2FoldChange`, `pvalue`, `padj`, `neg_log10_padj` |
| `neg_log10_padj` | Precomputed as −log₁₀(padj); used as the default volcano y-axis |
| `direction` | `Up in HD`, `Down in HD`, or `No significant change` |
| `significance` | `FDR<0.05` or `Not significant` |

### `gene_expression_long.csv`

This file is produced by the pipeline and used for reference only. The app does **not** accept it as an upload — it reconstructs the long-format data reactively from `counts_wide.csv` and `sample_info.csv`.

---

## Data Processing Workflow

The preprocessing pipeline (`scripts/02_prepare_processed_csvs.R`) performs the following steps:

1. **Metadata extraction** — Parses the GEO series matrix (`GSE64810_series_matrix.txt.gz`) to extract sample-level characteristics: diagnosis, age at death, PMI, RIN, CAG repeat, Vonsattel grade, and Hadzi–Vonsattel scores. Diagnosis labels are standardized to `Control` / `HD`.

2. **Counts matrix construction** — Reads the pre-computed DESeq2 normalized counts file (`GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz`). Sample columns are coerced to numeric and subset to match the 69 samples in the metadata. Gene symbols are joined from the DE table.

3. **Differential expression table** — Reads the pre-computed DESeq2 DE results (`GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz`). Adds derived columns: `neg_log10_padj`, `significance` (FDR < 0.05), `direction` (Up/Down in HD), and `rank` (sorted by padj ascending, then |log₂FC| descending).

4. **Long-format pivot** — Pivots the wide counts matrix to long format and joins sample metadata. Output: `gene_expression_long.csv` (~147 MB). Not used by the app at runtime.

5. **Validation** — `03_validate_processed_csvs.R` checks that all required columns are present, sample IDs match between counts and metadata, diagnosis labels are well-formed, and gene IDs are unique.

> **Normalization note:** Counts were normalized by the original study authors using DESeq2 size-factor normalization. This app does not re-run normalization. Do not apply additional normalization to the `counts_wide.csv` values.

> **DE method note:** Differential expression was computed by the original study authors using DESeq2 with outlier trimming (Cook's distance filtering). The comparison is HD vs. Control; positive log₂FC indicates higher expression in HD. This app displays those results — it does not re-run DESeq2.

---

## Screenshots / Demo

> Add screenshots after recording the submission video. Save images to `www/screenshots/`.

![Sample Information Tab](www/screenshots/tab1_sample_info.png)

![Counts Exploration — Diagnostic Plots](www/screenshots/tab2_diagnostic_plots.png)

![Counts Exploration — Heatmap](www/screenshots/tab2_heatmap.png)

![Counts Exploration — PCA](www/screenshots/tab2_pca.png)

![Differential Expression — Volcano Plot](www/screenshots/tab3_volcano.png)

![Gene Expression — Individual Gene](www/screenshots/tab4_gene_expression.png)

---

## Video Presentation

The submission requires a screen-recorded presentation of ≤ 5 minutes. Suggested outline:

| Segment | Duration | Content |
|---|---|---|
| Introduction | ~30 sec | State the biological question (HD vs. Control), dataset (GSE64810, 69 samples, BA9), and app structure (4 tabs) |
| Tab 1 — Sample Info | ~60 sec | Load `sample_info.csv`; show Summary sub-tab (row/column counts, column types); switch to Data Table; show a grouped violin plot (e.g., `age_at_death` grouped by `diagnosis`) |
| Tab 2 — Counts | ~90 sec | Load `counts_wide.csv`; explain the two filter sliders; show Filter Summary; walk through Diagnostic Plots and explain pass/fail coloring; show the heatmap and toggle log-transform; show PCA scatter with % variance labels |
| Tab 3 — DE | ~60 sec | Load `differential_expression.csv`; search for a gene (e.g., "HTT") in the table; switch to Volcano Plot; adjust colors and threshold slider; click "Update Volcano Plot" |
| Tab 4 — Gene Expression | ~60 sec | Load both files; search for "HTT"; select "diagnosis" as group; try all four plot types; explain the reactive table below the plot |

---

## Testing

The project includes a formal test suite of 81 assertions across four files, covering all pure helper functions in `app/app.R`. Tests do not launch the Shiny application.

### Run all tests

```r
library(testthat)
test_file("tests/test_tab1.R")   # 13 assertions
test_file("tests/test_tab2.R")   # 30 assertions
test_file("tests/test_tab3.R")   # 17 assertions
test_file("tests/test_tab4.R")   # 21 assertions
```

### What is tested

| File | Functions tested |
|---|---|
| `test_tab1.R` | `load_csv_from_path` (extension validation, empty-file errors), `build_col_summary` (numeric/categorical formatting), `get_numeric_cols`, `get_categorical_cols` |
| `test_tab2.R` | `strip_count_matrix`, `pass_var_filter`, `pass_nonzero_filter`, `compute_gene_stats`, `cap_heatmap_genes`; PCA sanity check (PC count, variance sums to 1) |
| `test_tab3.R` | `filter_de_table` (empty search, whitespace, case-insensitive partial match, no-match), `make_volcano` (returns gg object, correct axis labels, highlighted-flag threshold) |
| `test_tab4.R` | `get_gene_counts` (named numeric vector, missing-gene error), `pivot_gene_to_long` (columns present, values join correctly, order-independent), `make_gene_plot` (all four plot types return gg; labels match arguments) |

All 81 assertions pass with `FAIL 0 | WARN 0 | SKIP 0`.

---

## Known Limitations

- **Sex metadata is unavailable.** The GEO series matrix for GSE64810 does not include biological sex. The `sex` column exists in `sample_info.csv` but is all `NA`. Grouping by sex in Tab 1 or Tab 4 will produce an empty/NA-only plot.

- **HD-specific covariates are NA for controls.** `vonsattel_grade`, `cag_repeat`, `age_of_onset`, and `disease_duration_years` are meaningful only for HD samples. Grouping by these columns will show all controls as NA.

- **Large file load times.** `counts_wide.csv` is ~30 MB. Initial upload and parsing may take several seconds depending on hardware.

- **Heatmap is capped at 500 genes.** After filtering, only the top 500 highest-variance genes are passed to `pheatmap`. This is a deliberate performance trade-off; the full filtered set may contain thousands of genes.

- **Column name assumptions.** The app expects `counts_wide.csv` to have gene metadata in exactly columns 1–3 (`gene_id`, `gene_symbol`, `gene_name`) and numeric sample counts in columns 4 onward. Files not matching this structure will produce incorrect results or errors.

- **Tabs do not share loaded data.** Each tab has its own independent file input. Uploading `counts_wide.csv` in Tab 2 does not populate Tab 4; both files must be uploaded separately in Tab 4.

- **No download buttons.** Filtered tables, plots, and PCA results cannot currently be exported from within the app.

- **No deployment.** The app is designed to run locally. It has not been tested on shinyapps.io or a Shiny Server.

---

## Future Improvements

- **Cross-tab data sharing.** Allow a single file upload to propagate to all tabs that require the same file.
- **Downloadable results.** Add `downloadButton` for filtered DE tables, PCA scores, and plots.
- **Robust file validation.** Detect and report mismatched sample IDs between counts and metadata at upload time.
- **Gene set enrichment analysis.** Extend Tab 3 with an fgsea module using the DE rank column already present in `differential_expression.csv`.
- **Linked gene selection.** Clicking a point on the volcano plot could pre-populate the gene selector in Tab 4.
- **shinyapps.io deployment.** Publish a live demo with pre-loaded data.
- **Sex metadata recovery.** Cross-reference GEO sample accessions against SRA metadata to recover sex information.

---

## Citation / References

**Dataset:**
> TODO: Verify — Labadorf A, Hoss AG, Lagomarsino V, Latourelle JC, Hadzi TC, et al. (2015) RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. *PLOS ONE* 10(12): e0143563. https://doi.org/10.1371/journal.pone.0143563

**GEO entry:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810

**Key R packages:**

| Package | Purpose | Citation |
|---|---|---|
| [shiny](https://shiny.posit.co/) | Web application framework | Chang W et al., Posit Software |
| [bslib](https://rstudio.github.io/bslib/) | Bootstrap themes (Flatly) | Sievert C et al., Posit Software |
| [ggplot2](https://ggplot2.tidyverse.org/) | Data visualization | Wickham H (2016), Springer |
| [DT](https://rstudio.github.io/DT/) | Interactive data tables | Xie Y et al. |
| [pheatmap](https://cran.r-project.org/package=pheatmap) | Hierarchical heatmap | Kolde R |
| [ggbeeswarm](https://github.com/eclarke/ggbeeswarm) | Beeswarm plots | Clarke E, Sherrill-Mix S |
| [colourpicker](https://daattali.com/shiny/colourpicker/) | Color picker widget | Attali D |
| [RColorBrewer](https://colorbrewer2.org/) | Color palettes | Neuwirth E |

---

## Author

**Name:** TODO: add name  
**Email:** TODO: add email  
**Course:** BF591 — Bioinformatics with R, Boston University  
**Project type:** Final project — Integrated R Shiny Bioinformatics Application  
**Repository:** https://github.com/YZversion/APP
