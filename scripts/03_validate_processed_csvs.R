# Validate processed CSV files before they are used by the Shiny app.
#
# Run from the project root:
#   Rscript scripts/03_validate_processed_csvs.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

processed_dir <- file.path("data", "processed")

sample_info_file <- file.path(processed_dir, "sample_info.csv")
counts_wide_file <- file.path(processed_dir, "counts_wide.csv")
de_file <- file.path(processed_dir, "differential_expression.csv")
gene_expression_file <- file.path(processed_dir, "gene_expression_long.csv")
dictionary_file <- file.path(processed_dir, "app_data_dictionary.csv")

required_files <- c(
  sample_info_file,
  counts_wide_file,
  de_file,
  gene_expression_file,
  dictionary_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing processed file(s):\n", paste(missing_files, collapse = "\n"))
}

fail_if <- function(condition, message) {
  if (isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

cat("Reading processed CSV files...\n")
sample_info <- read_csv(sample_info_file, show_col_types = FALSE)
counts_wide <- read_csv(counts_wide_file, show_col_types = FALSE)
differential_expression <- read_csv(de_file, show_col_types = FALSE)
gene_expression_long <- read_csv(gene_expression_file, show_col_types = FALSE)

cat("Validating sample_info.csv...\n")
fail_if(!all(c("sample_id", "diagnosis") %in% names(sample_info)),
        "sample_info.csv must contain sample_id and diagnosis columns.")
fail_if(any(is.na(sample_info$sample_id) | sample_info$sample_id == ""),
        "sample_info.csv has missing sample_id values.")
fail_if(any(is.na(sample_info$diagnosis) | sample_info$diagnosis == ""),
        "sample_info.csv has missing diagnosis labels.")
fail_if(!all(sample_info$diagnosis %in% c("Control", "HD")),
        "sample_info.csv diagnosis values must be Control or HD.")
fail_if(anyDuplicated(sample_info$sample_id) > 0,
        "sample_info.csv has duplicated sample_id values.")

cat("Validating counts_wide.csv...\n")
fail_if(!"gene_id" %in% names(counts_wide),
        "counts_wide.csv must contain a gene_id column.")
fail_if(ncol(counts_wide) <= 4,
        "counts_wide.csv should have genes as rows and sample columns after gene metadata.")
fail_if(anyDuplicated(counts_wide$gene_id) > 0,
        "counts_wide.csv has duplicated gene_id values.")

metadata_columns <- intersect(c("gene_id", "gene_symbol", "gene_name"), names(counts_wide))
count_sample_columns <- setdiff(names(counts_wide), metadata_columns)
fail_if(length(count_sample_columns) == 0,
        "counts_wide.csv has no sample count columns.")

non_numeric_count_columns <- count_sample_columns[
  !vapply(counts_wide[count_sample_columns], is.numeric, logical(1))
]
fail_if(length(non_numeric_count_columns) > 0,
        paste("These count columns are not numeric:", paste(non_numeric_count_columns, collapse = ", ")))

cat("Validating differential_expression.csv...\n")
required_de_columns <- c("gene_id", "log2FoldChange", "pvalue", "padj")
fail_if(!all(required_de_columns %in% names(differential_expression)),
        paste("differential_expression.csv must contain:", paste(required_de_columns, collapse = ", ")))

cat("Validating sample ID matching between metadata and counts...\n")
missing_in_counts <- setdiff(sample_info$sample_id, count_sample_columns)
extra_in_counts <- setdiff(count_sample_columns, sample_info$sample_id)
fail_if(length(missing_in_counts) > 0,
        paste("Sample IDs missing from counts_wide.csv:", paste(missing_in_counts, collapse = ", ")))
fail_if(length(extra_in_counts) > 0,
        paste("Extra sample IDs in counts_wide.csv:", paste(extra_in_counts, collapse = ", ")))

cat("Validating gene_expression_long.csv...\n")
required_long_columns <- c(
  "gene_id",
  "gene_symbol",
  "sample_id",
  "diagnosis",
  "normalized_count",
  "log2_normalized_count"
)
fail_if(!all(required_long_columns %in% names(gene_expression_long)),
        paste("gene_expression_long.csv must contain:", paste(required_long_columns, collapse = ", ")))
fail_if(!all(unique(gene_expression_long$sample_id) %in% sample_info$sample_id),
        "gene_expression_long.csv contains sample IDs not found in sample_info.csv.")

cat("\nValidation passed.\n")
cat("Samples:", nrow(sample_info), "\n")
cat("Diagnosis counts:\n")
print(table(sample_info$diagnosis))
cat("Genes in counts_wide.csv:", nrow(counts_wide), "\n")
cat("Sample columns in counts_wide.csv:", length(count_sample_columns), "\n")
cat("Rows in differential_expression.csv:", nrow(differential_expression), "\n")
cat("Rows in gene_expression_long.csv:", nrow(gene_expression_long), "\n")

