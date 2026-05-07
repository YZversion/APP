library(testthat)
source(here::here("app/app.R"))

# ── Test 1: strip_count_matrix removes cols 1-3 and returns a numeric matrix ──
test_that("strip_count_matrix returns a numeric matrix with only sample columns", {
  df <- data.frame(
    gene_id     = c("g1", "g2", "g3"),
    gene_symbol = c("A",  "B",  "C"),
    gene_name   = c("x",  "y",  "z"),
    S1 = c(1.0,  2.0,  3.0),
    S2 = c(4.0,  5.0,  6.0),
    S3 = c(7.0,  8.0,  9.0),
    S4 = c(10.0, 11.0, 12.0),
    S5 = c(13.0, 14.0, 15.0),
    stringsAsFactors = FALSE
  )
  mat <- strip_count_matrix(df)

  expect_true(is.matrix(mat))
  expect_true(is.numeric(mat))
  expect_equal(nrow(mat), 3)
  expect_equal(ncol(mat), 5)
  expect_equal(colnames(mat), c("S1", "S2", "S3", "S4", "S5"))
})

# ── Test 2: pass_var_filter — exactly 3 high-variance genes pass at 70 pct ────
test_that("pass_var_filter selects the 3 high-variance genes at the 70th percentile", {
  mat <- matrix(0, nrow = 10, ncol = 6)
  rownames(mat) <- paste0("gene", 1:10)
  # Genes 1-7: low variance (alternating 1/2, var ≈ 0.3)
  for (i in 1:7)  mat[i, ] <- c(1, 2, 1, 2, 1, 2)
  # Genes 8-10: high variance (alternating 1/100, var ≈ 2940)
  for (i in 8:10) mat[i, ] <- c(1, 100, 1, 100, 1, 100)

  pass <- pass_var_filter(mat, var_pct = 70)

  expect_equal(sum(pass), 3L)
  expect_true(all(pass[8:10]))
  expect_false(any(pass[1:7]))
})

# ── Test 3: pass_nonzero_filter — exactly 4 genes meet the min non-zero cutoff ─
test_that("pass_nonzero_filter identifies genes with enough non-zero samples", {
  mat <- matrix(0, nrow = 10, ncol = 6)
  rownames(mat) <- paste0("gene", 1:10)
  # Genes 1-4: >= 4 non-zero samples
  mat[1, ] <- c(1, 2, 3, 4, 0, 0)   # 4 non-zero
  mat[2, ] <- c(1, 2, 3, 4, 5, 0)   # 5 non-zero
  mat[3, ] <- c(1, 2, 3, 4, 5, 6)   # 6 non-zero
  mat[4, ] <- c(0, 2, 3, 4, 5, 0)   # 4 non-zero
  # Genes 5-10: < 4 non-zero samples (0, 1, 2, or 3)
  mat[5, ] <- c(1, 2, 3, 0, 0, 0)   # 3
  mat[6, ] <- c(1, 2, 0, 0, 0, 0)   # 2
  mat[7, ] <- c(1, 0, 0, 0, 0, 0)   # 1
  mat[8, ] <- c(0, 0, 0, 0, 0, 0)   # 0
  mat[9, ] <- c(1, 2, 3, 0, 0, 0)   # 3
  mat[10,] <- c(0, 0, 1, 0, 0, 0)   # 1

  pass <- pass_nonzero_filter(mat, min_nonzero = 4)

  expect_equal(sum(pass), 4L)
  expect_true(all(pass[1:4]))
  expect_false(any(pass[5:10]))
})

# ── Test 4: combined filter is the intersection of both individual filters ─────
test_that("combined filter requires a gene to pass both var and nonzero filters", {
  mat <- rbind(
    c(1, 500, 1, 500, 1, 500),   # g1: highest var, 6 non-zero  → pass both
    c(1, 200, 1, 200, 1, 200),   # g2: high var,    6 non-zero  → pass both
    c(0, 100, 0, 100, 0,   0),   # g3: medium var,  2 non-zero  → pass var, fail nz
    c(1,   2, 1,   2, 1,   2),   # g4: low var,     6 non-zero  → fail var, pass nz
    c(1,   1, 1,   1, 1,   1),   # g5: zero var,    6 non-zero  → fail var, pass nz
    c(0,   0, 0,   0, 0,   0)    # g6: zero var,    0 non-zero  → fail both
  )
  rownames(mat) <- paste0("g", 1:6)

  pv  <- pass_var_filter(mat, var_pct = 70)
  pnz <- pass_nonzero_filter(mat, min_nonzero = 4)

  # g1 and g2 are the only ones with high enough variance
  expect_true(pv[1] && pv[2])
  expect_false(any(pv[3:6]))

  # g1, g2, g4, g5 have enough non-zero samples; g3 and g6 do not
  expect_true(all(pnz[c(1, 2, 4, 5)]))
  expect_false(pnz[3] || pnz[6])

  # Combined: only g1 and g2 pass both
  combined <- pv & pnz
  expect_equal(sum(combined), 2L)
  expect_true(combined[1] && combined[2])
  expect_false(any(combined[3:6]))
})

# ── Test 5: compute_gene_stats — log10 values match hand-computed expected ─────
test_that("compute_gene_stats returns numerically correct log10 median and variance", {
  # Single gene, 6 samples with known median and variance
  # median(c(0,1,4,9,16,25)) = (4+9)/2 = 6.5
  # var(c(0,1,4,9,16,25)) = 94.9667 (hand-computed)
  x   <- c(0, 1, 4, 9, 16, 25)
  mat <- matrix(x, nrow = 1, dimnames = list("g1", paste0("S", 1:6)))

  expected_log_med <- log10(median(x) + 1)
  expected_log_var <- log10(var(x)    + 1)

  stats <- compute_gene_stats(mat, pass = TRUE)

  expect_equal(stats$log_median[1], expected_log_med, tolerance = 1e-6)
  expect_equal(stats$log_var[1],    expected_log_var, tolerance = 1e-6)
  expect_equal(stats$n_zeros[1],    1L)                       # one zero in the vector
  expect_equal(as.character(stats$status[1]), "Pass")

  # Also verify a Fail gene is labelled correctly
  stats_fail <- compute_gene_stats(mat, pass = FALSE)
  expect_equal(as.character(stats_fail$status[1]), "Fail")
})

# ── Test 6: cap_heatmap_genes — returns 500 rows, all with top variance ────────
test_that("cap_heatmap_genes returns exactly 500 rows ordered by decreasing variance", {
  set.seed(42)
  # 600-gene × 6-sample matrix; multiply row i by i to give row i variance ∝ i²
  base_mat <- matrix(runif(600 * 6, min = 0, max = 10), nrow = 600, ncol = 6)
  mat_big  <- base_mat * (1:600)   # row-wise scaling (R recycles down columns)
  rownames(mat_big) <- paste0("gene", 1:600)

  capped <- cap_heatmap_genes(mat_big, n_max = 500)

  # Exactly 500 rows returned
  expect_equal(nrow(capped), 500L)

  # Every returned gene has variance >= the 501st-highest variance in the full set
  vars_full  <- apply(mat_big, 1, var)
  thresh_501 <- sort(vars_full, decreasing = TRUE)[501]
  vars_capped <- apply(capped, 1, var)
  expect_true(all(vars_capped >= thresh_501))

  # The overall top gene (highest variance) must be present; the bottom must not
  top_gene <- rownames(mat_big)[which.max(vars_full)]
  bot_gene <- rownames(mat_big)[which.min(vars_full)]
  expect_true(top_gene  %in% rownames(capped))
  expect_false(bot_gene %in% rownames(capped))
})

# ── Test 7: PCA sanity — correct number of PCs and variance sums to 1 ─────────
test_that("prcomp on a 6-sample x 20-gene matrix has correct PC count and variance", {
  set.seed(99)
  mat_pca <- matrix(runif(20 * 6), nrow = 20, ncol = 6)
  rownames(mat_pca) <- paste0("gene",   1:20)
  colnames(mat_pca) <- paste0("sample", 1:6)

  # t(mat_pca) is 6 samples × 20 genes — same call used in the app
  pca <- prcomp(t(mat_pca), scale. = TRUE)

  # prcomp returns min(n, p) components regardless of centering — centering
  # reduces the *effective* rank by 1 but prcomp still emits min(n,p) columns
  # in x (the last will have near-zero variance).  For a 6×20 matrix: min(6,20)=6.
  expected_pcs <- min(nrow(t(mat_pca)), ncol(t(mat_pca)))
  expect_equal(ncol(pca$x), expected_pcs)

  # Proportion of variance across all PCs must sum to 1
  pct_var <- pca$sdev^2 / sum(pca$sdev^2)
  expect_equal(sum(pct_var), 1.0, tolerance = 1e-6)

  # Cross-check via summary()$importance row "Proportion of Variance"
  importance <- summary(pca)$importance
  expect_equal(sum(importance["Proportion of Variance", ]), 1.0, tolerance = 1e-6)
})
