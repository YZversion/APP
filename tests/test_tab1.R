library(testthat)
source(here::here("app/app.R"))

# ── Test 1: valid CSV returns a data frame with correct dimensions ─────────────
test_that("load_csv_from_path returns a data frame with correct dimensions", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(
    data.frame(x = 1:3, y = c("a", "b", "c"), z = c(1.1, 2.2, 3.3)),
    tmp, row.names = FALSE
  )
  df <- load_csv_from_path(tmp)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_equal(ncol(df), 3)
})

# ── Test 2: wrong extension errors with a message containing "csv" ────────────
test_that("load_csv_from_path errors on a .txt file with a message mentioning csv", {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  write.csv(data.frame(a = 1:3, b = 4:6, c = 7:9), tmp, row.names = FALSE)
  expect_error(load_csv_from_path(tmp), regexp = "csv", ignore.case = TRUE)
})

# ── Test 3: single data row errors ────────────────────────────────────────────
test_that("load_csv_from_path errors when the file has only 1 data row", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(data.frame(a = 1, b = 2, c = 3), tmp, row.names = FALSE)
  expect_error(load_csv_from_path(tmp), regexp = "more than 1 data row")
})

# ── Test 4: build_col_summary returns correct Type and mean (sd) format ────────
test_that("build_col_summary Type column and mean (sd) format are correct", {
  df <- data.frame(
    score = c(10, 20, 30),
    label = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  summ <- build_col_summary(df)

  num_row  <- summ[summ[["Column Name"]] == "score", ]
  char_row <- summ[summ[["Column Name"]] == "label", ]

  # class() can return "numeric" or "integer" depending on how the column was created
  expect_true(num_row$Type %in% c("numeric", "integer"))
  expect_equal(char_row$Type, "character")

  # Numeric summary must match sprintf("%.2f (%.2f)", mean, sd)
  expected <- sprintf("%.2f (%.2f)", mean(c(10, 20, 30)), sd(c(10, 20, 30)))
  expect_equal(num_row[["Mean (sd) or Distinct Values"]], expected)
})

# ── Test 5: get_categorical_cols excludes all-NA columns ──────────────────────
test_that("get_categorical_cols excludes all-NA columns but keeps valid categoricals", {
  df <- data.frame(
    diagnosis = c("HD", "Control", "HD"),
    score     = c(1.0, 2.0, 3.0),
    sex       = c(NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  cat_cols <- get_categorical_cols(df)

  expect_true("diagnosis" %in% cat_cols)   # valid categorical — must appear
  expect_false("score" %in% cat_cols)      # numeric — must not appear
  expect_false("sex" %in% cat_cols)        # all-NA — must not appear
})

# ── Test 6: get_numeric_cols returns only numeric columns ──────────────────────
test_that("get_numeric_cols returns only numeric columns from a mixed data frame", {
  df <- data.frame(
    age   = c(30L, 40L, 50L),
    group = c("A", "B", "C"),
    value = c(1.1, 2.2, 3.3),
    stringsAsFactors = FALSE
  )
  num_cols <- get_numeric_cols(df)

  expect_setequal(num_cols, c("age", "value"))
  expect_false("group" %in% num_cols)
})
