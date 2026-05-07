# BF591 Final Project: Huntington's Disease BA9 RNA-seq

This project analyzes post-mortem human prefrontal cortex Brodmann Area 9 RNA-seq data from Huntington's disease and neurologically normal control subjects.

Primary dataset:
- GEO: GSE64810
- SRA: SRP051844
- Tissue: post-mortem human prefrontal cortex / BA9
- Samples: 20 Huntington's disease, 49 controls
- Platform: Illumina HiSeq 2000

The Shiny app will use prepared CSV files in `data/processed/`. The app itself is intentionally not implemented yet.

## Folder Structure

```text
.
├── app/                  # Future R Shiny app files
├── data/
│   ├── raw/              # Downloaded GEO/SRA files, unchanged
│   ├── metadata/         # Cleaned sample annotations and lookup files
│   └── processed/        # Final CSV files consumed by Shiny
├── docs/                 # Project notes and workflow documentation
├── results/              # Plots/tables generated during analysis
└── scripts/              # Data download, cleaning, preprocessing scripts
```

## First Milestone

1. Download the GEO processed files and sample metadata.
2. Clean sample information into a single sample-level CSV.
3. Convert counts into Shiny-friendly long and wide CSV formats.
4. Run or import differential expression results.
5. Save final app-ready CSV files under `data/processed/`.

See `docs/data_preparation_plan.md` for the detailed workflow and CSV schemas.

