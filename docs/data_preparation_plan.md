# Data Preparation Plan

## Dataset Choice

Use GEO series `GSE64810`: mRNA-seq expression profiling of human post-mortem BA9 brain tissue for Huntington's disease and neurologically normal individuals.

Reasons:
- Directly matches the project theme: Huntington's disease, human post-mortem prefrontal cortex / BA9, RNA-seq.
- Has a manageable final-project sample size: 69 total samples, 20 HD and 49 control.
- GEO provides processed normalized counts and a DESeq2 result table, so the project can focus on Shiny exploration and interpretable RNA-seq analysis without requiring full FASTQ alignment.

Key links:
- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810
- SRA study: https://www.ncbi.nlm.nih.gov/sra?term=SRP051844

## Files To Download Or Prepare

Place original downloaded files in `data/raw/`.

Required:
- `GSE64810_series_matrix.txt.gz`
  - Source: GEO Series Matrix File.
  - Purpose: sample metadata, GSM accessions, sample titles, characteristics.
- `GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz`
  - Source: GEO supplementary file.
  - Purpose: normalized gene counts matrix.
- `GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz`
  - Source: GEO supplementary file.
  - Purpose: published DESeq2 differential expression table.

Recommended helper files:
- Gene annotation table from Bioconductor, Ensembl BioMart, or org.Hs.eg.db.
  - Purpose: map Ensembl IDs to gene symbols and gene descriptions if the count table does not already include clear symbols.
- `metadata_manual_notes.csv`
  - Purpose: manually record any metadata decisions, renamed columns, ambiguous sample labels, or excluded samples.

Optional advanced raw-data route:
- SRA run table from SRP051844.
- FASTQ files from SRA.
- Reference genome and gene annotation, for example GRCh38 FASTA and GTF.
- This is not necessary for the first Shiny version unless the instructor requires raw-read preprocessing.

## Recommended Project Folder Structure

```text
data/raw/
  GSE64810_series_matrix.txt.gz
  GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz
  GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz

data/metadata/
  samples_clean.csv
  gene_annotation.csv
  metadata_manual_notes.csv

data/processed/
  sample_info.csv
  counts_wide.csv
  counts_long.csv
  counts_summary_by_sample.csv
  counts_summary_by_gene.csv
  differential_expression.csv
  gene_expression_long.csv
  app_data_dictionary.csv

scripts/
  01_download_data.R
  02_clean_metadata.R
  03_prepare_counts.R
  04_differential_expression.R
  05_export_shiny_csvs.R
```

## Preprocessing Workflow

### 1. Download Data

Download GEO processed files into `data/raw/`. Keep these files unchanged for reproducibility.

Expected sources:
- GEO series matrix for metadata.
- GEO supplementary normalized counts.
- GEO supplementary DESeq2 results.

### 2. Clean Sample Metadata

Parse `GSE64810_series_matrix.txt.gz` and extract one row per sample.

Minimum required metadata:
- GEO sample accession, for example `GSM1580869`.
- Original sample ID, for example `C_0002` or `H_0001`.
- Diagnosis group: `Control` or `HD`.
- Tissue/region: `BA9 prefrontal cortex`.

Recommended metadata if present in GEO characteristics:
- Sex.
- Age at death.
- RNA integrity number / RIN.
- Post-mortem interval / PMI.
- Vonsattel grade.
- CAG repeat size.
- Batch/library variables if available.

Quality checks:
- Confirm exactly 69 samples.
- Confirm 49 controls and 20 HD.
- Confirm sample IDs in metadata match count matrix columns.
- Standardize missing values as blank or `NA`, not mixed strings.

### 3. Prepare Counts Matrix

Input: `GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz`.

Steps:
- Read the normalized count matrix.
- Standardize the gene ID column name.
- Add gene symbols if needed.
- Ensure all sample columns match `sample_info.csv`.
- Convert numeric columns to numeric.
- Remove duplicate gene rows only if they are exact duplicates; otherwise keep stable gene IDs as the unique key.
- Create both wide and long count tables.

Quality checks:
- No duplicated `gene_id` values in final app tables.
- No sample columns missing from metadata.
- No all-zero or all-missing sample columns.
- Record whether values are raw counts, normalized counts, or transformed counts.

### 4. Differential Expression

Two acceptable routes:

Default route for the first app:
- Import GEO's published DESeq2 result table.
- Standardize column names.
- Add significance labels using clear thresholds, for example `padj < 0.05` and direction based on `log2FoldChange`.

More reproducible route:
- Use the normalized count table for exploration only.
- If raw integer counts are available or generated, run DESeq2 locally with design `~ diagnosis`.
- Export the local DE table with the same schema as below.

For BF591, the default route is acceptable if clearly documented as using GEO's processed DESeq2 output.

### 5. Export App-Ready CSVs

All final Shiny inputs should live in `data/processed/`. These CSVs should have stable column names and no spaces in column names.

## Final CSV Schemas For Shiny

### `sample_info.csv`

Used by:
- Sample Information Exploration tab.
- Filtering all other tabs by diagnosis, sex, age, etc.

Columns:

```text
sample_id
geo_accession
diagnosis
tissue_region
sex
age_at_death
pmi_hours
rin
vonsattel_grade
cag_repeat
source
notes
```

Example rows:

```csv
sample_id,geo_accession,diagnosis,tissue_region,sex,age_at_death,pmi_hours,rin,vonsattel_grade,cag_repeat,source,notes
C_0002,GSM1580869,Control,BA9 prefrontal cortex,F,72,18.5,7.2,NA,NA,GSE64810,
H_0001,GSM1580918,HD,BA9 prefrontal cortex,M,59,20.0,6.9,3,44,GSE64810,
```

### `counts_wide.csv`

Used by:
- Counts Matrix Exploration tab.
- Heatmaps, PCA, sample correlation, gene lookup.

Columns:

```text
gene_id
gene_symbol
gene_name
C_0002
C_0003
...
H_0750
```

Example:

```csv
gene_id,gene_symbol,gene_name,C_0002,C_0003,H_0001,H_0002
ENSG00000197386,HTT,huntingtin,1234.5,1188.2,1401.7,1320.1
ENSG00000141510,TP53,tumor protein p53,51.2,47.9,63.0,58.4
```

### `counts_long.csv`

Used by:
- Interactive filtering and plotting.
- Individual Gene Expression Visualization tab.

Columns:

```text
gene_id
gene_symbol
sample_id
normalized_count
log2_normalized_count
diagnosis
```

Example:

```csv
gene_id,gene_symbol,sample_id,normalized_count,log2_normalized_count,diagnosis
ENSG00000197386,HTT,C_0002,1234.5,10.27,Control
ENSG00000197386,HTT,H_0001,1401.7,10.45,HD
```

### `counts_summary_by_sample.csv`

Used by:
- Counts Matrix Exploration tab.
- QC-style sample overview plots.

Columns:

```text
sample_id
diagnosis
total_normalized_count
median_normalized_count
detected_genes
percent_zero_or_missing
```

Example:

```csv
sample_id,diagnosis,total_normalized_count,median_normalized_count,detected_genes,percent_zero_or_missing
C_0002,Control,25433210.4,18.2,22801,4.6
H_0001,HD,24992114.8,17.8,22640,5.1
```

### `counts_summary_by_gene.csv`

Used by:
- Counts Matrix Exploration tab.
- Gene filtering controls.

Columns:

```text
gene_id
gene_symbol
mean_control
mean_hd
median_control
median_hd
mean_all
detected_samples
detected_percent
```

Example:

```csv
gene_id,gene_symbol,mean_control,mean_hd,median_control,median_hd,mean_all,detected_samples,detected_percent
ENSG00000197386,HTT,1210.4,1355.8,1198.7,1341.2,1252.6,69,100
```

### `differential_expression.csv`

Used by:
- Differential Expression tab.
- Volcano plot, sortable DE table, top genes.

Columns:

```text
gene_id
gene_symbol
gene_name
baseMean
log2FoldChange
lfcSE
stat
pvalue
padj
neg_log10_padj
significance
direction
rank
```

Rules:
- `significance`: `FDR<0.05` or `Not significant`.
- `direction`: `Up in HD`, `Down in HD`, or `No significant change`.
- `rank`: sort by `padj` ascending, then absolute `log2FoldChange` descending.

Example:

```csv
gene_id,gene_symbol,gene_name,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,neg_log10_padj,significance,direction,rank
ENSG00000197386,HTT,huntingtin,1280.2,0.21,0.08,2.62,0.0088,0.041,1.39,FDR<0.05,Up in HD,152
```

### `gene_expression_long.csv`

Used by:
- Optional Individual Gene Expression Visualization tab.
- Boxplots, jitter plots, and per-gene summary statistics.

This can be the same as `counts_long.csv`, or a filtered version containing genes present in the DE table.

Columns:

```text
gene_id
gene_symbol
sample_id
diagnosis
normalized_count
log2_normalized_count
age_at_death
sex
vonsattel_grade
cag_repeat
```

### `app_data_dictionary.csv`

Used by:
- Developer documentation.
- Optional app help text or validation.

Columns:

```text
file_name
column_name
description
data_type
allowed_values
missing_value_rule
```

Example:

```csv
file_name,column_name,description,data_type,allowed_values,missing_value_rule
sample_info.csv,diagnosis,Disease group,categorical,Control|HD,Required
differential_expression.csv,padj,Benjamini-Hochberg adjusted p-value,numeric,0-1,May be NA for low-information genes
```

## Shiny Tab To CSV Mapping

Sample Information Exploration:
- `sample_info.csv`

Counts Matrix Exploration:
- `counts_wide.csv`
- `counts_summary_by_sample.csv`
- `counts_summary_by_gene.csv`

Differential Expression:
- `differential_expression.csv`

Optional Individual Gene Expression Visualization:
- `gene_expression_long.csv`
- `sample_info.csv`

## Notes For The Final Writeup

Important limitations to state:
- This is post-mortem bulk tissue, so expression reflects mixed cell types and disease-associated tissue composition changes.
- The first app version uses processed GEO data rather than raw FASTQ alignment.
- Differential expression should be interpreted as association with HD diagnosis, not causal effect.
- Metadata completeness may limit covariate-adjusted modeling.

