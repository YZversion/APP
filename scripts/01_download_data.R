# Download raw GEO files for the Huntington's Disease BA9 RNA-seq project.
#
# Run from the project root:
#   Rscript scripts/01_download_data.R

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

files_to_download <- data.frame(
  file_name = c(
    "GSE64810_series_matrix.txt.gz",
    "GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz",
    "GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz"
  ),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64810/matrix/GSE64810_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64810/suppl/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64810/suppl/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz"
  ),
  stringsAsFactors = FALSE
)

cat("Checking required raw files in", raw_dir, "\n\n")

for (i in seq_len(nrow(files_to_download))) {
  destination <- file.path(raw_dir, files_to_download$file_name[i])

  if (file.exists(destination) && file.info(destination)$size > 0) {
    cat("Already exists:", destination, "\n")
    next
  }

  cat("Downloading:", files_to_download$file_name[i], "\n")
  download_ok <- tryCatch(
    {
      download.file(
        url = files_to_download$url[i],
        destfile = destination,
        mode = "wb",
        quiet = FALSE
      )
      TRUE
    },
    error = function(e) {
      message("Download failed for ", files_to_download$file_name[i])
      message("Reason: ", conditionMessage(e))
      FALSE
    }
  )

  if (!download_ok || !file.exists(destination) || file.info(destination)$size == 0) {
    cat("\nManual download needed.\n")
    cat("Open this URL in a browser:\n")
    cat(files_to_download$url[i], "\n")
    cat("Save the file as:\n")
    cat(normalizePath(destination, winslash = "/", mustWork = FALSE), "\n\n")
  }
}

missing_files <- files_to_download$file_name[
  !file.exists(file.path(raw_dir, files_to_download$file_name)) |
    file.info(file.path(raw_dir, files_to_download$file_name))$size == 0
]

if (length(missing_files) > 0) {
  cat("\nThe following files are still missing:\n")
  cat(paste("-", missing_files), sep = "\n")
  stop("\nDownload the missing files into data/raw/ before running preprocessing.")
}

cat("\nAll required raw files are present.\n")

