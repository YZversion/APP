library(testthat)
source(here::here("app/app.R"))

# ── Test 1: filter_de_table — empty search returns the full data frame ────────
test_that("filter_de_table with empty search returns all rows unchanged", {
  de <- data.frame(
    gene_symbol = c("HTT", "MAPT", "SNCA"),
    padj        = c(0.001, 0.05, 0.5),
    stringsAsFactors = FALSE
  )
  result <- filter_de_table(de, "")

  expect_equal(nrow(result), 3L)
  expect_equal(result, de)
})

# ── Test 2: filter_de_table — whitespace-only search treated as empty ─────────
test_that("filter_de_table treats whitespace-only search as empty and returns all rows", {
  de <- data.frame(
    gene_symbol = c("HTT", "MAPT"),
    padj        = c(0.001, 0.5),
    stringsAsFactors = FALSE
  )
  result <- filter_de_table(de, "   ")

  expect_equal(nrow(result), 2L)
})

# ── Test 3: filter_de_table — case-insensitive partial match ──────────────────
test_that("filter_de_table matches gene symbols case-insensitively and partially", {
  de <- data.frame(
    gene_symbol = c("MAPT", "MAP2", "SNCA", "map3k1"),
    padj        = c(0.001, 0.01, 0.5, 0.02),
    stringsAsFactors = FALSE
  )
  result <- filter_de_table(de, "map")

  # MAPT, MAP2, map3k1 all contain "map" (case-insensitive); SNCA does not
  expect_equal(nrow(result), 3L)
  expect_true(all(c("MAPT", "MAP2", "map3k1") %in% result$gene_symbol))
  expect_false("SNCA" %in% result$gene_symbol)
})

# ── Test 4: filter_de_table — no match returns zero rows ─────────────────────
test_that("filter_de_table returns 0 rows when no gene symbol matches", {
  de <- data.frame(
    gene_symbol = c("HTT", "MAPT", "SNCA"),
    padj        = c(0.001, 0.05, 0.5),
    stringsAsFactors = FALSE
  )
  result <- filter_de_table(de, "XYZ_NOMATCH_9999")

  expect_equal(nrow(result), 0L)
})

# ── Test 5: make_volcano — returns a ggplot with correct axis labels ──────────
test_that("make_volcano returns a gg object with x/y labels matching the chosen columns", {
  de <- data.frame(
    log2FoldChange = c(-3, -1,  0,  2,  4),
    neg_log10_padj = c( 8,  1,  0,  2, 10),
    pvalue         = c(0.001, 0.1, 1, 0.05, 0.0001),
    stringsAsFactors = FALSE
  )
  p <- make_volcano(de,
                    x_col         = "log2FoldChange",
                    y_col         = "neg_log10_padj",
                    base_col      = "#999999",
                    highlight_col = "#d7191c",
                    padj_thresh   = 5)

  expect_s3_class(p, "gg")
  expect_equal(p$labels$x, "log2FoldChange")
  expect_equal(p$labels$y, "neg_log10_padj")
})

# ── Test 6: make_volcano — highlighted flag matches the padj threshold exactly ─
test_that("make_volcano flags exactly the genes whose neg_log10_padj exceeds the threshold", {
  # neg_log10_padj values: 8, 1, 0, 2, 10
  # threshold = 5  →  genes at rows 1 (8>5) and 5 (10>5) should be TRUE; others FALSE
  de <- data.frame(
    log2FoldChange = c(-3, -1,  0,  2,  4),
    neg_log10_padj = c( 8,  1,  0,  2, 10),
    gene_symbol    = c("A", "B", "C", "D", "E"),
    stringsAsFactors = FALSE
  )
  p <- make_volcano(de,
                    x_col         = "log2FoldChange",
                    y_col         = "neg_log10_padj",
                    base_col      = "#999999",
                    highlight_col = "#d7191c",
                    padj_thresh   = 5)

  plot_data <- p$data
  expect_true("highlighted" %in% names(plot_data))

  # Exactly 2 genes are highlighted (rows 1 and 5)
  expect_equal(sum(plot_data$highlighted), 2L)
  expect_true( plot_data$highlighted[plot_data$gene_symbol == "A"])
  expect_true( plot_data$highlighted[plot_data$gene_symbol == "E"])
  expect_false(plot_data$highlighted[plot_data$gene_symbol == "B"])
  expect_false(plot_data$highlighted[plot_data$gene_symbol == "C"])
  expect_false(plot_data$highlighted[plot_data$gene_symbol == "D"])
})
