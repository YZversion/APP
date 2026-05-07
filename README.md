# HD RNA-seq Explorer

**An integrated R Shiny application for interactive exploration of Huntington's Disease transcriptomic data.**

BF591 — Bioinformatics with R | Boston University | Final Project

---

## Overview

HD RNA-seq Explorer is a single, self-contained R Shiny application that provides four integrated interactive modules for exploring a publicly available RNA-seq dataset comparing post-mortem brain tissue from Huntington's Disease (HD) patients and neurologically normal controls.

| Module | Purpose |
|---|---|
| Sample Information | Explore clinical and technical metadata for all 69 samples |
| Counts Matrix | Filter, cluster, and decompose the normalized gene expression matrix |
| Differential Expression | Browse and visualize pre-computed DESeq2 results interactively |
| Gene Expression *(optional module)* | Examine any single gene's expression across user-defined sample groups |

All analysis is performed reactively in-browser — no command-line R or pre-generated plots are required to use the app.

---

## Quick Start

```r
# 1. Clone the repository
#    (bash / PowerShell terminal)
git clone https://github.com/YZversion/APP.git
cd APP

# 2. Install R dependencies  (run once in R or RStudio)
install.packages(c(
  "shiny", "bslib", "ggplot2", "DT", "tidyr",
  "pheatmap", "RColorBrewer", "ggbeeswarm", "colourpicker"
))

# 3. Generate processed data files  (~2 min; downloads ~50 MB from NCBI GEO)
Rscript scripts/run_all_data_prep.R

# 4. Launch the app
shiny::runApp("app/app.R")
```

Upload the CSV files from `data/processed/` to each tab when prompted.

> `data/raw/` and `data/processed/` are not committed to this repository. They are generated locally by the preprocessing pipeline.

---

## Project Requirements Coverage

This table maps the app implementation to the BF591 Final Project rubric.

| Requirement | Implementation | Status |
|---|---|---|
| Single integrated app (not four separate apps) | One `app/app.R`; `navbarPage` with four `tabPanel`s | ✅ |
| Sample Information Exploration | Tab 1: Summary, Data Table, Distributions | ✅ |
| Counts Matrix Exploration | Tab 2: Filter Summary, Diagnostic Plots, Heatmap, PCA | ✅ |
| Differential Expression | Tab 3: Sortable table with gene search, volcano plot | ✅ |
| Choose-your-own module | **Tab 4: Individual Gene Expression Visualization** | ✅ |
| File validation | `.csv` extension check; row/column count check; inline `validate()` messages | ✅ |
| Labeled buttons and input descriptions | All controls use `helpText()` and labeled `actionButton`s | ✅ |
| Unit test suite | 81 assertions across 4 `testthat` files; all passing | ✅ |

### Optional Module Choice

Tab 4 implements **Individual Gene Expression Visualization**, the "choose-your-own" module for this project. It accepts `counts_wide.csv` and `sample_info.csv`, lets the user search for any gene by symbol, choose a grouping variable from the sample metadata, and render the expression distribution as a boxplot, violin, bar chart, or beeswarm plot. A reactive data table of normalized counts joined to sample metadata is shown below the plot.

---

## Dataset

| Property | Value |
|---|---|
| GEO Accession | [GSE64810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810) |
| SRA Study | SRP051844 |
| Tissue | Post-mortem human prefrontal cortex (Brodmann Area 9) |
| Sequencing platform | Illumina HiSeq 2000 |
| Total samples | 69 (20 HD, 49 neurologically normal controls) |
| Comparison direction | HD vs. Control — positive log₂FC = higher in HD |
| Normalization | DESeq2 size-factor normalization (pre-computed by the original study authors; this app does not re-run DESeq2) |

> **Citation (TODO: verify):** Labadorf A, Hoss AG, Lagomarsino V, Latourelle JC, Hadzi TC, et al. (2015) RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. *PLOS ONE* 10(12): e0143563. <https://doi.org/10.1371/journal.pone.0143563>

### Sample Metadata (`sample_info.csv`)

| Column | Type | Notes |
|---|---|---|
| `sample_id` | character | `C_####` (control) or `H_####` (HD) |
| `geo_accession` | character | GEO sample accession (GSM…) |
| `diagnosis` | categorical | `Control` or `HD` |
| `age_at_death` | numeric | Years |
| `pmi_hours` | numeric | Post-mortem interval in hours |
| `rin` | numeric | RNA integrity number (0–10) |
| `mrna_seq_reads` | numeric | mRNA-seq read count |
| `vonsattel_grade` | numeric | HD neuropathology grade — **NA for all controls** |
| `cag_repeat` | numeric | HTT CAG repeat length — **NA for all controls** |
| `hv_striatal_score` | numeric | Hadzi–Vonsattel striatal score |
| `hv_cortical_score` | numeric | Hadzi–Vonsattel cortical score |
| `sex` | — | Not available in the GEO series matrix for this accession (all NA; automatically excluded from grouping dropdowns) |

### Processed Files

| File | Size | Description |
|---|---|---|
| `sample_info.csv` | < 1 MB | 69 rows × 19 columns of sample metadata |
| `counts_wide.csv` | ~30 MB | Genes × samples normalized counts (cols 1–3 gene metadata; cols 4–72 samples) |
| `differential_expression.csv` | ~6 MB | Pre-computed DESeq2 results with derived columns |
| `gene_expression_long.csv` | ~147 MB | Long-format pivot — **pipeline output only; not loaded by the app** |
| `app_data_dictionary.csv` | < 1 MB | Column-level documentation for all files |

---

## App Modules

### Tab 1 — Sample Information Exploration

**Input:** `sample_info.csv`

**Sub-tabs:**

| Sub-tab | Output |
|---|---|
| Summary | Row/column count; column type; mean ± SD (numeric) or up to 5 distinct values (categorical) |
| Data Table | Full sample metadata — sortable, filterable `DT` table |
| Distributions | Grouped violin + jitter when a categorical column is selected; histogram when no group is chosen |

**Controls (sidebar):** Numeric column selector; optional group-by selector. Both are populated dynamically from the uploaded file; the group-by selector automatically excludes all-NA columns (e.g., `sex`).

*Interpretation note: Use PMI and RIN distributions to assess potential technical confounders before interpreting expression results.*

---

### Tab 2 — Counts Matrix Exploration

**Input:** `counts_wide.csv`

**Controls:**

| Control | Default | Effect |
|---|---|---|
| Variance percentile | 50 % | Retains genes at or above the Nth percentile of per-gene variance |
| Minimum non-zero samples | 10 | Retains genes with non-zero counts in at least N samples |
| Log₂-transform toggle | On | Applies log₂(x + 1) to the heatmap matrix |
| PCA display type | Scatter | Scatter (PC vs. PC) or beeswarm (top N PCs) |

**Sub-tabs:**

| Sub-tab | Output |
|---|---|
| Filter Summary | Passing/failing gene counts and percentages; total sample count |
| Diagnostic Plots | (1) log₁₀ median vs. log₁₀ variance; (2) log₁₀ median vs. zero-sample count — pass/fail colored |
| Heatmap | `pheatmap` hierarchical heatmap, RdBu palette, capped at 500 highest-variance passing genes |
| PCA | `prcomp` on transposed filtered matrix (samples as rows); axis labels show % variance explained |

*Interpretation note: Filtered genes in the diagnostic plots are rendered lighter and drawn first so passing genes are visible on top. Low-expression, high-zero genes removed here are typically technical noise.*

---

### Tab 3 — Differential Expression

**Input:** `differential_expression.csv` (pre-computed DESeq2 results — the app does not re-run DESeq2)

**Controls:**

| Control | Description |
|---|---|
| X-axis / Y-axis selectors | Choose any DE column for each axis (default: `log2FoldChange` vs. `neg_log10_padj`) |
| Base / Highlight color pickers | `colourpicker` widgets for non-significant and significant points |
| Threshold slider | Points where −log₁₀(padj) exceeds this value are drawn in the highlight color |
| Gene symbol search | Case-insensitive partial-match filter applied to the table |
| Update Volcano Plot button | Plot renders only on click — prevents re-rendering on every slider adjustment |

**Sub-tabs:**

| Sub-tab | Output |
|---|---|
| Table | Filterable, sortable `DT` table of DE results |
| Volcano Plot | ggplot2 scatter with user-defined axes, colors, and significance threshold |

---

### Tab 4 — Individual Gene Expression Visualization *(optional module)*

**Inputs:** `counts_wide.csv` + `sample_info.csv` (uploaded independently from Tab 2)

**Controls:**

| Control | Description |
|---|---|
| Gene selector | Searchable `selectizeInput` populated from all gene symbols in the counts matrix |
| Group by | Dropdown of categorical columns from `sample_info.csv` (defaults to `diagnosis`) |
| Plot type | Boxplot, Violin, Bar (mean ± SE), or Beeswarm |
| Plot button | Plot renders on click via `bindEvent()` |

**Outputs:**
- Expression plot colored by group
- Reactive `DT` table showing normalized counts joined to sample metadata for the selected gene

*Interpretation note: The bar plot shows group mean ± SE. HD-specific columns (`vonsattel_grade`, `cag_repeat`) are NA for all controls; grouping by `diagnosis` is recommended before exploring HD-specific covariates.*

---

## Repository Structure

```
.
├── app/
│   └── app.R                        # Complete Shiny app (~720 lines, all four tabs)
├── data/
│   ├── raw/                         # gitignored — downloaded GEO source files
│   ├── metadata/
│   │   └── samples_clean.csv        # committed — 69-row sample annotation
│   └── processed/                   # gitignored — app-ready CSVs (generated by pipeline)
├── docs/
│   └── data_preparation_plan.md     # Planning document (CSV schemas, workflow notes)
├── results/                         # gitignored — analysis outputs
├── scripts/
│   ├── 01_download_data.R           # Downloads 3 GEO files into data/raw/
│   ├── 02_prepare_processed_csvs.R  # Parses GEO files → 5 app-ready CSVs
│   ├── 03_validate_processed_csvs.R # Validates columns, sample IDs, and row counts
│   └── run_all_data_prep.R          # Runs all three steps in sequence
├── tests/
│   ├── test_tab1.R                  # 13 assertions
│   ├── test_tab2.R                  # 30 assertions
│   ├── test_tab3.R                  # 17 assertions
│   └── test_tab4.R                  # 21 assertions
├── .gitignore
├── CLAUDE.md                        # Developer notes (not for end users)
├── LICENSE
└── README.md
```

Only `data/metadata/samples_clean.csv` and the scripts/tests/app are committed. Raw GEO files and processed CSVs must be regenerated locally.

---

## Installation and Requirements

### R version

R ≥ 4.2.0 recommended. Tested on R 4.x with RStudio on macOS and Windows.

### Required packages

```r
# Shiny app
install.packages(c(
  "shiny",
  "bslib",
  "ggplot2",
  "DT",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "ggbeeswarm",
  "colourpicker"
))

# Preprocessing pipeline (run_all_data_prep.R only)
install.packages(c("dplyr", "readr"))

# Test suite
install.packages(c("testthat", "here"))
```

### Bioconductor

No Bioconductor packages are required to run the Shiny app or the preprocessing pipeline. The normalized counts and differential expression results are downloaded directly from NCBI GEO as pre-computed files produced by the original study authors. If you wish to independently reproduce the DESeq2 analysis, you would need `BiocManager::install("DESeq2")`, but this is outside the scope of the current workflow.

---

## How to Run the App

### Step 1 — Clone the repository

```bash
git clone https://github.com/YZversion/APP.git
cd APP
```

### Step 2 — Generate processed data files

Run from the **project root** (RStudio terminal, bash, or PowerShell). Downloads ~50 MB from NCBI GEO and writes five CSVs to `data/processed/`.

```r
Rscript scripts/run_all_data_prep.R
```

To run steps individually:

```r
Rscript scripts/01_download_data.R          # Download GEO files
Rscript scripts/02_prepare_processed_csvs.R  # Build app-ready CSVs (~2 min)
Rscript scripts/03_validate_processed_csvs.R # Verify all outputs
```

### Step 3 — Launch the app

From the project root in R or RStudio:

```r
shiny::runApp("app/app.R")
```

The app opens in your default browser. Upload the appropriate CSV from `data/processed/` using the file input in each tab.

---

## Input File Formats

All file inputs accept `.csv` only. The app validates the extension, minimum row count (> 1), and minimum column count (> 1), and shows an inline error message on any failure.

### `sample_info.csv`

| Requirement | Detail |
|---|---|
| `sample_id` | Required; unique; format `C_####` or `H_####` |
| `diagnosis` | Required; values must be exactly `Control` or `HD` |
| Rows | One per sample |

### `counts_wide.csv`

| Requirement | Detail |
|---|---|
| Column 1 | `gene_id` — Ensembl ID |
| Column 2 | `gene_symbol` — HGNC gene symbol (used as row names for gene lookup) |
| Column 3 | `gene_name` — description (may be empty) |
| Columns 4+ | One numeric column per sample; names must match `sample_id` in `sample_info.csv` |
| Values | DESeq2 size-factor normalized counts (non-negative) |

Columns 1–3 are stripped automatically before any computation.

### `differential_expression.csv`

| Requirement | Detail |
|---|---|
| Required columns | `gene_id`, `gene_symbol`, `log2FoldChange`, `pvalue`, `padj`, `neg_log10_padj` |
| `neg_log10_padj` | Pre-computed as −log₁₀(padj); default volcano y-axis |
| `significance` | `FDR<0.05` or `Not significant` |
| `direction` | `Up in HD`, `Down in HD`, or `No significant change` |

---

## Data Processing Workflow

The preprocessing pipeline (`scripts/02_prepare_processed_csvs.R`) performs these steps:

1. **Metadata extraction** — Parses the GEO series matrix to extract diagnosis, age at death, PMI, RIN, CAG repeat, Vonsattel grade, and Hadzi–Vonsattel scores. Diagnosis labels are standardized to `Control` / `HD`.

2. **Counts matrix construction** — Reads the pre-computed normalized counts file; coerces sample columns to numeric; subsets to the 69 samples in the metadata; joins gene symbols from the DE table.

3. **Differential expression table** — Reads the pre-computed DESeq2 results file; adds derived columns: `neg_log10_padj`, `significance` (FDR < 0.05), `direction` (Up/Down in HD), and `rank` (sorted by padj ascending, then |log₂FC| descending).

4. **Long-format pivot** — Pivots the wide matrix to long format and joins sample metadata. Writes `gene_expression_long.csv` (~147 MB) for reference. This file is **not loaded by the app at runtime**.

5. **Validation** — Checks required columns, unique gene IDs, matching sample IDs between counts and metadata, and well-formed diagnosis labels.

> **Normalization:** Counts were normalized by the original study authors using DESeq2 size-factor normalization. Do not apply additional normalization to `counts_wide.csv`.

> **Differential expression:** Results were computed by the original study authors using DESeq2 (Wald test, Benjamini–Hochberg correction, Cook's distance outlier filtering). The comparison is HD vs. Control. This app displays those results; it does not re-run DESeq2.

---

## Screenshots

Screenshots have not been captured yet. After recording the submission video, save images to `www/screenshots/` and replace the placeholders below.

| Tab | Placeholder path |
|---|---|
| Sample Information | `www/screenshots/tab1_sample_info.png` |
| Counts — Diagnostic Plots | `www/screenshots/tab2_diagnostic_plots.png` |
| Counts — Heatmap | `www/screenshots/tab2_heatmap.png` |
| Counts — PCA | `www/screenshots/tab2_pca.png` |
| Differential Expression — Volcano | `www/screenshots/tab3_volcano.png` |
| Gene Expression | `www/screenshots/tab4_gene_expression.png` |

---

## Video Presentation

Submission requires a screen-recorded presentation of ≤ 5 minutes.

| Segment | Time | Suggested content |
|---|---|---|
| Introduction | ~30 s | Biological question (HD vs. Control), dataset (GSE64810, 69 samples, BA9 cortex), app overview |
| Tab 1 — Sample Info | ~60 s | Load `sample_info.csv`; show Summary; switch to Data Table; show a violin of `age_at_death` by `diagnosis` |
| Tab 2 — Counts | ~90 s | Load `counts_wide.csv`; adjust sliders; walk through Filter Summary and Diagnostic Plots; show heatmap with log-transform on/off; show PCA with % variance labels |
| Tab 3 — DE | ~60 s | Load `differential_expression.csv`; search for "HTT" in the table; switch to Volcano Plot; change colors and threshold; click Update |
| Tab 4 — Gene Expression | ~60 s | Load both files; search for "HTT"; cycle through all four plot types; point out the reactive table |

---

## Testing

The project includes a formal test suite covering all pure helper functions extracted from `app/app.R`. Tests use `testthat` and source `app/app.R` directly — they do not launch the Shiny server.

### Run all tests

```r
library(testthat)
test_file("tests/test_tab1.R")   # 13 assertions
test_file("tests/test_tab2.R")   # 30 assertions
test_file("tests/test_tab3.R")   # 17 assertions
test_file("tests/test_tab4.R")   # 21 assertions
```

All 81 assertions pass: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 81 ]`

### Functions covered

| Test file | Functions tested |
|---|---|
| `test_tab1.R` | `load_csv_from_path` (extension check, malformed-file errors), `build_col_summary`, `get_numeric_cols`, `get_categorical_cols` |
| `test_tab2.R` | `strip_count_matrix`, `pass_var_filter`, `pass_nonzero_filter`, combined filter intersection, `compute_gene_stats`, `cap_heatmap_genes`; PCA variance check |
| `test_tab3.R` | `filter_de_table` (empty, whitespace, case-insensitive, no-match), `make_volcano` (gg object, axis labels, highlighted-flag logic) |
| `test_tab4.R` | `get_gene_counts` (named vector, missing-gene error), `pivot_gene_to_long` (columns, join correctness, order independence), `make_gene_plot` (all four plot types, axis labels) |

---

## Known Limitations

- **Sex metadata is unavailable.** `sex` is all-NA in the GEO series matrix for this accession. The app's `get_categorical_cols()` helper automatically excludes all-NA columns, so `sex` does not appear in any grouping dropdown.

- **HD-specific covariates are NA for all controls.** `vonsattel_grade`, `cag_repeat`, `age_of_onset`, and `disease_duration_years` are meaningful only for HD samples. They do appear in the grouping dropdown; selecting one will show controls as `NA`.

- **Heatmap capped at 500 genes.** After filtering, only the top 500 highest-variance genes are passed to `pheatmap`. This is an intentional performance trade-off.

- **File structure assumptions.** The app expects `counts_wide.csv` to have gene metadata in exactly columns 1–3 and numeric sample counts in columns 4+. Files not matching this layout will produce incorrect results.

- **Tabs do not share loaded data.** Each tab has its own independent file input. Uploading `counts_wide.csv` in Tab 2 does not populate Tab 4.

- **No export buttons.** Filtered tables and plots cannot currently be downloaded from within the app.

- **Local deployment only.** The app has not been tested on shinyapps.io or Shiny Server. The `counts_wide.csv` file (~30 MB) may be slow to upload in a remote session.

---

## Future Improvements

- Cross-tab data sharing so a file uploaded once is available in all tabs that need it
- Download buttons for filtered tables, PCA scores, and plots
- Upstream file validation: detect mismatched sample IDs between counts and metadata at upload time
- GSEA extension using the `rank` column already present in `differential_expression.csv`
- Volcano-to-Tab-4 linking: clicking a point pre-populates the gene selector
- shinyapps.io deployment with pre-loaded demo data

---

## Citations

**Dataset:**
> TODO: Verify — Labadorf A, Hoss AG, Lagomarsino V, Latourelle JC, Hadzi TC, et al. (2015) RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. *PLOS ONE* 10(12): e0143563. <https://doi.org/10.1371/journal.pone.0143563>

**GEO entry:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810>

**Key R packages:**

| Package | Use in this app |
|---|---|
| [shiny](https://shiny.posit.co/) | Web application framework |
| [bslib](https://rstudio.github.io/bslib/) | Bootstrap 5 / Flatly theme |
| [ggplot2](https://ggplot2.tidyverse.org/) | All plots |
| [DT](https://rstudio.github.io/DT/) | Interactive data tables |
| [pheatmap](https://cran.r-project.org/package=pheatmap) | Hierarchical heatmap |
| [ggbeeswarm](https://github.com/eclarke/ggbeeswarm) | Beeswarm plots (Tab 2 PCA and Tab 4) |
| [colourpicker](https://daattali.com/shiny/colourpicker/) | Color pickers in Tab 3 |
| [RColorBrewer](https://colorbrewer2.org/) | RdBu heatmap palette |

---

## Author

**Name:** TODO: add your name  
**Course:** BF591 — Bioinformatics with R, Boston University  
**Project:** Final project — Integrated R Shiny Bioinformatics Application  
**Repository:** <https://github.com/YZversion/APP>
