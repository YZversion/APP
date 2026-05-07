# Development Difficulties and Fixes

This file summarizes the main issues encountered while finishing and testing the HD RNA-seq Explorer Shiny app.

## Large CSV Upload Limit

**Problem:** Uploading `counts_wide.csv` failed because the file is about 30 MB, larger than Shiny's default upload limit.

**Fix:** Added a larger request size limit near the top of `app/app.R`:

```r
options(shiny.maxRequestSize = 200 * 1024^2)
```

This allows uploads up to 200 MB, which covers `counts_wide.csv` and other generated processed files used during local testing.

## Tab 4 Gene Selector Freezing

**Problem:** After uploading `counts_wide.csv` in Tab 4, the app became slow or appeared unresponsive. The Plot button looked like it could not be clicked.

**Cause:** The Tab 4 gene selector initially sent all gene symbols to the browser at once. With tens of thousands of genes, the client-side `selectizeInput` became too heavy and froze the page.

**Fix:** Changed the gene selector to use server-side selectize:

```r
selectizeInput("tab4_gene", "Select gene:", choices = NULL)

observeEvent(tab4_counts(), {
  genes <- rownames(tab4_counts())
  updateSelectizeInput(session, "tab4_gene",
                       choices = genes,
                       selected = genes[1],
                       server = TRUE)
})
```

This keeps the large gene list on the server and only sends matching choices to the browser during search.

## Tab 4 Functional Validation

**Problem:** Tab 4 still needed final local validation after the other three tabs were confirmed.

**Checks performed:**

- Uploaded `counts_wide.csv` and `sample_info.csv`.
- Searched for `HTT`.
- Grouped expression by `diagnosis`.
- Confirmed the joined table had 69 rows, matching all samples.
- Confirmed both `Control` and `HD` groups were present.
- Verified all four plot types worked: boxplot, violin, bar, and beeswarm.
- Ran `tests/test_tab4.R`, which passed all 21 assertions.

## Full Test Suite Verification

**Problem:** After changing upload behavior and Tab 4 selector behavior, the full helper test suite needed to be rechecked.

**Result:** All four test files passed:

- `tests/test_tab1.R`: 13 assertions
- `tests/test_tab2.R`: 30 assertions
- `tests/test_tab3.R`: 17 assertions
- `tests/test_tab4.R`: 21 assertions

Total: 81 passing assertions.

## README Cleanup

**Problem:** The README still had unfinished wording from earlier drafts, including a citation TODO and unclear testing language.

**Fix:** Replaced the TODO with the full Labadorf et al. PLOS ONE citation and clarified the testing section so it accurately describes the 81-test assertion suite.

## RStudio Project File

**Problem:** The repository did not initially have an `.Rproj` file, making the project less convenient to open consistently in RStudio.

**Fix:** Created `APP.Rproj` so the app can be opened directly as an RStudio project.

## Git Hygiene

**Problem:** Several local-only files and generated data directories appeared during testing, including:

- `.claude/`
- `data/raw/`
- `data/processed/`
- local Shiny log files

**Resolution:** Only app code changes were committed when pushing the upload-size fix. Generated data remains untracked because it is reproducible from the preprocessing scripts and is excluded by project policy.
