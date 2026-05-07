# Prepare app-ready CSV files from GSE64810 raw GEO files.
#
# Run from the project root:
#   Rscript scripts/02_prepare_processed_csvs.R
#
# Outputs:
#   data/processed/sample_info.csv
#   data/processed/counts_wide.csv
#   data/processed/differential_expression.csv
#   data/processed/gene_expression_long.csv
#   data/processed/app_data_dictionary.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

raw_dir <- file.path("data", "raw")
metadata_dir <- file.path("data", "metadata")
processed_dir <- file.path("data", "processed")

dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

series_matrix_file <- file.path(raw_dir, "GSE64810_series_matrix.txt.gz")
counts_file <- file.path(raw_dir, "GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz")
de_file <- file.path(raw_dir, "GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz")

required_raw <- c(series_matrix_file, counts_file, de_file)
missing_raw <- required_raw[!file.exists(required_raw)]
if (length(missing_raw) > 0) {
  stop(
    "Missing raw file(s):\n",
    paste(missing_raw, collapse = "\n"),
    "\nRun scripts/01_download_data.R first."
  )
}

clean_geo_value <- function(x) {
  x <- gsub('^"|"$', "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "na")] <- NA_character_
  x
}

to_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

standardize_diagnosis <- function(x) {
  dplyr::case_when(
    grepl("Huntington", x, ignore.case = TRUE) ~ "HD",
    grepl("normal|control", x, ignore.case = TRUE) ~ "Control",
    TRUE ~ NA_character_
  )
}

extract_characteristic <- function(characteristic_matrix, field_name) {
  pattern <- paste0("^", field_name, ":")
  matches <- apply(characteristic_matrix, 2, function(values) {
    hit <- values[grepl(pattern, values, ignore.case = TRUE)]
    if (length(hit) == 0) {
      return(NA_character_)
    }
    sub(pattern, "", hit[1], ignore.case = TRUE)
  })
  clean_geo_value(matches)
}

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_lines <- lines[grepl("^!Sample_", lines)]

  split_line <- function(line) {
    clean_geo_value(strsplit(line, "\t", fixed = TRUE)[[1]])
  }

  title_line <- split_line(sample_lines[grepl("^!Sample_title\t", sample_lines)][1])
  geo_line <- split_line(sample_lines[grepl("^!Sample_geo_accession\t", sample_lines)][1])
  source_line <- split_line(sample_lines[grepl("^!Sample_source_name_ch1\t", sample_lines)][1])
  characteristic_lines <- sample_lines[grepl("^!Sample_characteristics_ch1\t", sample_lines)]

  sample_id <- title_line[-1]
  geo_accession <- geo_line[-1]
  source_name <- source_line[-1]

  characteristic_matrix <- do.call(
    rbind,
    lapply(characteristic_lines, function(line) split_line(line)[-1])
  )

  tissue_region <- extract_characteristic(characteristic_matrix, "tissue")
  diagnosis_raw <- extract_characteristic(characteristic_matrix, "diagnosis")

  sample_info <- tibble(
    sample_id = sample_id,
    geo_accession = geo_accession,
    diagnosis = standardize_diagnosis(diagnosis_raw),
    diagnosis_raw = diagnosis_raw,
    tissue_region = tissue_region,
    sex = NA_character_,
    age_at_death = to_number(extract_characteristic(characteristic_matrix, "age of death")),
    pmi_hours = to_number(extract_characteristic(characteristic_matrix, "pmi")),
    rin = to_number(extract_characteristic(characteristic_matrix, "rin")),
    mrna_seq_reads = to_number(extract_characteristic(characteristic_matrix, "mrna-seq reads")),
    age_of_onset = to_number(extract_characteristic(characteristic_matrix, "age of onset")),
    disease_duration_years = to_number(extract_characteristic(characteristic_matrix, "duration")),
    vonsattel_grade = to_number(extract_characteristic(characteristic_matrix, "vonsattel grade")),
    cag_repeat = to_number(extract_characteristic(characteristic_matrix, "cag")),
    hv_striatal_score = to_number(extract_characteristic(characteristic_matrix, "h-v striatal score")),
    hv_cortical_score = to_number(extract_characteristic(characteristic_matrix, "h-v cortical score")),
    source_name = source_name,
    source = "GSE64810",
    notes = NA_character_
  )

  sample_info
}

cat("Parsing sample metadata...\n")
sample_info <- parse_series_matrix(series_matrix_file)

write_csv(sample_info, file.path(metadata_dir, "samples_clean.csv"), na = "")
write_csv(sample_info, file.path(processed_dir, "sample_info.csv"), na = "")

cat("Reading differential expression table...\n")
de_raw <- read_tsv(de_file, show_col_types = FALSE, name_repair = "minimal")
names(de_raw)[1] <- "gene_id"

differential_expression <- de_raw %>%
  rename(gene_symbol = symbol) %>%
  mutate(
    gene_name = NA_character_,
    neg_log10_padj = if_else(!is.na(padj) & padj > 0, -log10(padj), NA_real_),
    significance = if_else(!is.na(padj) & padj < 0.05, "FDR<0.05", "Not significant"),
    direction = case_when(
      significance == "FDR<0.05" & log2FoldChange > 0 ~ "Up in HD",
      significance == "FDR<0.05" & log2FoldChange < 0 ~ "Down in HD",
      TRUE ~ "No significant change"
    )
  ) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  mutate(rank = row_number()) %>%
  select(
    gene_id,
    gene_symbol,
    gene_name,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,
    neg_log10_padj,
    significance,
    direction,
    rank,
    everything()
  )

write_csv(differential_expression, file.path(processed_dir, "differential_expression.csv"), na = "")

gene_lookup <- differential_expression %>%
  select(gene_id, gene_symbol, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)

cat("Reading normalized counts matrix...\n")
counts_raw <- read_tsv(counts_file, show_col_types = FALSE, name_repair = "minimal")
names(counts_raw)[1] <- "gene_id"

sample_ids <- sample_info$sample_id
count_sample_ids <- setdiff(names(counts_raw), "gene_id")
missing_in_counts <- setdiff(sample_ids, count_sample_ids)
extra_in_counts <- setdiff(count_sample_ids, sample_ids)

if (length(missing_in_counts) > 0 || length(extra_in_counts) > 0) {
  stop(
    "Sample IDs do not match between metadata and counts.\n",
    "Missing in counts: ", paste(missing_in_counts, collapse = ", "), "\n",
    "Extra in counts: ", paste(extra_in_counts, collapse = ", ")
  )
}

counts_wide <- counts_raw %>%
  mutate(across(all_of(sample_ids), as.numeric)) %>%
  left_join(gene_lookup, by = "gene_id") %>%
  mutate(
    gene_symbol = if_else(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    gene_name = if_else(is.na(gene_name), "", gene_name)
  ) %>%
  select(gene_id, gene_symbol, gene_name, all_of(sample_ids))

write_csv(counts_wide, file.path(processed_dir, "counts_wide.csv"), na = "")

cat("Creating long gene expression table. This may take a moment...\n")
gene_expression_long <- counts_wide %>%
  select(gene_id, gene_symbol, all_of(sample_ids)) %>%
  pivot_longer(
    cols = all_of(sample_ids),
    names_to = "sample_id",
    values_to = "normalized_count"
  ) %>%
  left_join(
    sample_info %>%
      select(sample_id, diagnosis, age_at_death, sex, vonsattel_grade, cag_repeat),
    by = "sample_id"
  ) %>%
  mutate(log2_normalized_count = log2(normalized_count + 1)) %>%
  select(
    gene_id,
    gene_symbol,
    sample_id,
    diagnosis,
    normalized_count,
    log2_normalized_count,
    age_at_death,
    sex,
    vonsattel_grade,
    cag_repeat
  )

write_csv(gene_expression_long, file.path(processed_dir, "gene_expression_long.csv"), na = "")

cat("Writing data dictionary...\n")
app_data_dictionary <- tribble(
  ~file_name, ~column_name, ~description, ~data_type, ~allowed_values, ~missing_value_rule,
  "sample_info.csv", "sample_id", "Study sample identifier used in counts matrix", "character", "C_#### or H_####", "Required",
  "sample_info.csv", "geo_accession", "GEO sample accession", "character", "GSM accession", "Required",
  "sample_info.csv", "diagnosis", "Simplified disease group", "categorical", "Control|HD", "Required",
  "sample_info.csv", "diagnosis_raw", "Original diagnosis text from GEO", "character", "Neurologically normal|Huntington's Disease", "Required",
  "sample_info.csv", "tissue_region", "Tissue described in GEO metadata", "categorical", "BA9 prefrontal cortex", "Required",
  "sample_info.csv", "sex", "Biological sex if available", "categorical", "Not available in GEO series matrix", "May be missing",
  "sample_info.csv", "age_at_death", "Age at death in years", "numeric", "Positive number", "May be missing",
  "sample_info.csv", "pmi_hours", "Post-mortem interval in hours", "numeric", "Positive number", "May be missing",
  "sample_info.csv", "rin", "RNA integrity number", "numeric", "0-10", "May be missing",
  "sample_info.csv", "mrna_seq_reads", "Number of mRNA-seq reads reported in GEO", "numeric", "Non-negative integer", "May be missing",
  "sample_info.csv", "age_of_onset", "HD age of onset in years", "numeric", "Positive number", "Expected missing for controls",
  "sample_info.csv", "disease_duration_years", "HD disease duration in years", "numeric", "Positive number", "Expected missing for controls",
  "sample_info.csv", "vonsattel_grade", "HD neuropathology grade", "numeric", "3|4 in this dataset", "Expected missing for controls",
  "sample_info.csv", "cag_repeat", "HTT CAG repeat count", "numeric", "Positive number", "Expected missing for controls",
  "counts_wide.csv", "gene_id", "Ensembl gene identifier with version", "character", "ENSG...", "Required",
  "counts_wide.csv", "gene_symbol", "Gene symbol from GEO DESeq2 table when available", "character", "Gene symbol or gene_id fallback", "Required",
  "counts_wide.csv", "gene_name", "Gene description if available", "character", "Not available from provided GEO files", "May be missing",
  "counts_wide.csv", "<sample columns>", "DESeq2 normalized count for each sample", "numeric", "Non-negative number", "Required",
  "differential_expression.csv", "gene_id", "Ensembl gene identifier with version", "character", "ENSG...", "Required",
  "differential_expression.csv", "gene_symbol", "Gene symbol", "character", "Gene symbol", "May be missing",
  "differential_expression.csv", "baseMean", "Mean normalized expression across samples", "numeric", "Non-negative number", "May be missing",
  "differential_expression.csv", "log2FoldChange", "DESeq2 log2 fold change for HD vs Control", "numeric", "Any real number", "Required",
  "differential_expression.csv", "pvalue", "DESeq2 Wald test p-value", "numeric", "0-1", "Required",
  "differential_expression.csv", "padj", "Benjamini-Hochberg adjusted p-value", "numeric", "0-1", "Required, but may be NA for low-information genes",
  "differential_expression.csv", "significance", "FDR significance label using padj < 0.05", "categorical", "FDR<0.05|Not significant", "Required",
  "differential_expression.csv", "direction", "Direction among significant genes", "categorical", "Up in HD|Down in HD|No significant change", "Required",
  "gene_expression_long.csv", "gene_id", "Ensembl gene identifier with version", "character", "ENSG...", "Required",
  "gene_expression_long.csv", "gene_symbol", "Gene symbol from GEO DESeq2 table when available", "character", "Gene symbol or gene_id fallback", "Required",
  "gene_expression_long.csv", "sample_id", "Study sample identifier", "character", "Must match sample_info.csv", "Required",
  "gene_expression_long.csv", "diagnosis", "Simplified disease group", "categorical", "Control|HD", "Required",
  "gene_expression_long.csv", "normalized_count", "DESeq2 normalized count", "numeric", "Non-negative number", "Required",
  "gene_expression_long.csv", "log2_normalized_count", "log2(normalized_count + 1)", "numeric", "Non-negative number", "Required"
)

write_csv(app_data_dictionary, file.path(processed_dir, "app_data_dictionary.csv"), na = "")

cat("\nWrote processed files to", processed_dir, "\n")
cat("- sample_info.csv\n")
cat("- counts_wide.csv\n")
cat("- differential_expression.csv\n")
cat("- gene_expression_long.csv\n")
cat("- app_data_dictionary.csv\n")
