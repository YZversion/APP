library(testthat)
source(here::here("app/app.R"))

# ── Test 1: get_gene_counts — returns correct named numeric vector ─────────────
test_that("get_gene_counts returns a named numeric vector for a valid gene symbol", {
  mat <- matrix(
    c(10, 20, 30, 40, 50, 60),
    nrow = 2,
    dimnames = list(c("HTT", "MAPT"), c("S1", "S2", "S3"))
  )
  result <- get_gene_counts(mat, "HTT")

  expect_true(is.numeric(result))
  expect_equal(length(result), 3L)
  expect_equal(names(result), c("S1", "S2", "S3"))
  expect_equal(as.numeric(result), c(10, 30, 50))  # matrix fills by column
})

# ── Test 2: get_gene_counts — stops with a "not found" error for absent gene ──
test_that("get_gene_counts stops with an informative error when gene is absent", {
  mat <- matrix(1:6, nrow = 2,
                dimnames = list(c("HTT", "MAPT"), c("S1", "S2", "S3")))
  expect_error(get_gene_counts(mat, "NONEXISTENT_GENE"), regexp = "not found")
})

# ── Test 3: pivot_gene_to_long — result has correct columns and row count ──────
test_that("pivot_gene_to_long returns sample_id, normalized_count, and group_col", {
  count_vec <- c(S1 = 10.5, S2 = 20.0, S3 = 5.5)
  sample_df <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    diagnosis = c("HD", "Control", "HD"),
    stringsAsFactors = FALSE
  )
  result <- pivot_gene_to_long(count_vec, sample_df, "diagnosis")

  expect_equal(nrow(result), 3L)
  expect_true(all(c("sample_id", "normalized_count", "diagnosis") %in% names(result)))
})

# ── Test 4: pivot_gene_to_long — values join correctly by sample_id ───────────
test_that("pivot_gene_to_long correctly joins count values to metadata by sample_id", {
  count_vec <- c(S1 = 10.5, S2 = 20.0, S3 = 5.5)
  sample_df <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    diagnosis = c("HD", "Control", "HD"),
    stringsAsFactors = FALSE
  )
  result <- pivot_gene_to_long(count_vec, sample_df, "diagnosis")

  s1 <- result[result$sample_id == "S1", ]
  s2 <- result[result$sample_id == "S2", ]
  expect_equal(s1$normalized_count, 10.5)
  expect_equal(s1$diagnosis, "HD")
  expect_equal(s2$normalized_count, 20.0)
  expect_equal(s2$diagnosis, "Control")
})

# ── Test 5: pivot_gene_to_long — all samples present even when order differs ──
test_that("pivot_gene_to_long includes all samples regardless of input order", {
  count_vec <- c(S3 = 5.5, S1 = 10.5, S2 = 20.0)  # deliberately shuffled
  sample_df <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    diagnosis = c("HD", "Control", "HD"),
    stringsAsFactors = FALSE
  )
  result <- pivot_gene_to_long(count_vec, sample_df, "diagnosis")

  expect_equal(nrow(result), 3L)
  expect_true(all(c("S1", "S2", "S3") %in% result$sample_id))
})

# ── Test 6: make_gene_plot — returns a gg object for all four plot types ───────
test_that("make_gene_plot returns a gg object for box, violin, bar, and beeswarm", {
  long_df <- data.frame(
    sample_id        = paste0("S", 1:6),
    normalized_count = c(10, 12, 15, 20, 22, 25),
    diagnosis        = c("HD", "HD", "HD", "Control", "Control", "Control"),
    stringsAsFactors = FALSE
  )
  for (pt in c("box", "violin", "bar", "beeswarm")) {
    p <- make_gene_plot(long_df, "HTT", "diagnosis", pt)
    expect_s3_class(p, "gg")
  }
})

# ── Test 7: make_gene_plot — axis and title labels match gene_name/group_col ──
test_that("make_gene_plot sets x, y, and title labels from its arguments", {
  long_df <- data.frame(
    sample_id        = paste0("S", 1:4),
    normalized_count = c(10, 15, 20, 25),
    diagnosis        = c("HD", "HD", "Control", "Control"),
    stringsAsFactors = FALSE
  )
  p <- make_gene_plot(long_df, "HTT", "diagnosis", "violin")

  expect_equal(p$labels$x, "diagnosis")
  expect_equal(p$labels$y, "Normalized Count")
  expect_true(grepl("HTT",       p$labels$title))
  expect_true(grepl("diagnosis", p$labels$title))
})
