# Run the full data preparation workflow.
#
# Run from the project root:
#   Rscript scripts/run_all_data_prep.R

cat("Step 1: Download/check raw data\n")
source(file.path("scripts", "01_download_data.R"))

cat("\nStep 2: Prepare processed CSV files\n")
source(file.path("scripts", "02_prepare_processed_csvs.R"))

cat("\nStep 3: Validate processed CSV files\n")
source(file.path("scripts", "03_validate_processed_csvs.R"))

cat("\nData preparation complete.\n")

