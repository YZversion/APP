# Scripts

Run these from the project root in the RStudio terminal.

```bash
Rscript scripts/run_all_data_prep.R
```

Individual scripts:

1. `01_download_data.R`: downloads/checks the GEO series matrix and supplementary files.
2. `02_prepare_processed_csvs.R`: creates app-ready CSV files from raw GEO files.
3. `03_validate_processed_csvs.R`: validates that required CSV files, columns, and sample IDs match.
4. `run_all_data_prep.R`: runs all three steps in order.

Processed outputs:

- `data/processed/sample_info.csv`
- `data/processed/counts_wide.csv`
- `data/processed/differential_expression.csv`
- `data/processed/gene_expression_long.csv`
- `data/processed/app_data_dictionary.csv`
