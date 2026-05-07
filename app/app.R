library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggbeeswarm)
library(colourpicker)

# ── Pure helpers (no Shiny context — testable standalone) ──────────────────────

# Reads and validates a CSV from a file path.
# Raises stop() with a plain message on any failure so tests can use expect_error().
# skip_ext_check = TRUE when the caller has already validated the original filename
# (e.g. from fileInput$name in the Shiny wrapper below).
load_csv_from_path <- function(filepath, skip_ext_check = FALSE) {
  if (!skip_ext_check) {
    ext <- tolower(tools::file_ext(filepath))
    if (ext != "csv") stop(paste0("Expected a .csv file; got .", ext, "."))
  }
  df <- tryCatch(
    read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Could not parse the file. Check that it is a valid CSV.")
  )
  if (nrow(df) <= 1) stop("File must have more than 1 data row.")
  if (ncol(df) <= 1) stop("File must have more than 1 column.")
  df
}

# Shiny wrapper around load_csv_from_path.
# Validates the original filename (file_input$name) before delegating to the
# path helper, since datapath is a temp file that may lack the original extension.
load_csv <- function(file_input) {
  validate(need(!is.null(file_input), "Upload a CSV file to begin."))
  validate(need(
    tolower(tools::file_ext(file_input$name)) == "csv",
    paste0("Expected a .csv file; got .", tools::file_ext(file_input$name), ".")
  ))
  tryCatch(
    load_csv_from_path(file_input$datapath, skip_ext_check = TRUE),
    error = function(e) validate(need(FALSE, conditionMessage(e)))
  )
}

# Returns names of numeric columns in df.
get_numeric_cols <- function(df) {
  names(df)[sapply(df, is.numeric)]
}

# Returns names of non-numeric columns, excluding all-NA columns.
# All-NA columns have no valid grouping values and would produce empty dropdowns.
get_categorical_cols <- function(df) {
  is_non_numeric <- !sapply(df, is.numeric)
  is_not_all_na  <- !sapply(df, function(x) all(is.na(x)))
  names(df)[is_non_numeric & is_not_all_na]
}

# Builds the per-column summary data frame for the Summary sub-tab.
# Numeric columns: "mean (sd)". Character columns: up to 5 distinct values.
build_col_summary <- function(df) {
  data.frame(
    `Column Name` = names(df),
    Type = sapply(df, function(x) class(x)[1]),
    `Mean (sd) or Distinct Values` = sapply(df, function(x) {
      if (is.numeric(x)) {
        sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
      } else {
        vals <- unique(x[!is.na(x) & nzchar(x)])
        if (length(vals) <= 5) {
          paste(vals, collapse = ", ")
        } else {
          paste0(paste(vals[1:5], collapse = ", "), " ... (", length(vals), " total)")
        }
      }
    }),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# --- tab2 helpers (also used by tests) ---------------------------------------

# Strips gene_id / gene_symbol / gene_name (cols 1-3) from counts_wide.csv
# and coerces the remaining sample columns to a plain numeric matrix.
strip_count_matrix <- function(df) {
  mat <- as.matrix(df[, -(1:3)])
  storage.mode(mat) <- "numeric"
  mat
}

# Returns a logical vector (length = nrow(mat)): TRUE where per-gene variance
# is at or above the var_pct-th percentile of all gene variances.
pass_var_filter <- function(mat, var_pct) {
  gene_vars <- apply(mat, 1, var, na.rm = TRUE)
  thresh    <- quantile(gene_vars, probs = var_pct / 100, na.rm = TRUE)
  gene_vars >= thresh
}

# Returns a logical vector: TRUE where a gene has >= min_nonzero samples with
# a non-zero count.
pass_nonzero_filter <- function(mat, min_nonzero) {
  rowSums(mat > 0, na.rm = TRUE) >= min_nonzero
}

# Builds the per-gene stats data frame used by both diagnostic scatter plots.
# pass is a logical vector of length nrow(mat).
compute_gene_stats <- function(mat, pass) {
  data.frame(
    log_median = log10(apply(mat, 1, median, na.rm = TRUE) + 1),
    log_var    = log10(apply(mat, 1, var,    na.rm = TRUE) + 1),
    n_zeros    = rowSums(mat == 0, na.rm = TRUE),
    # Factor levels set so Fail rows sort first (drawn behind Pass points)
    status     = factor(ifelse(pass, "Pass", "Fail"), levels = c("Fail", "Pass"))
  )
}

# Subsets a passed-filter matrix to the top n_max genes by variance.
# Called before pheatmap to keep rendering responsive on large datasets.
cap_heatmap_genes <- function(mat_pass, n_max = 500) {
  if (nrow(mat_pass) <= n_max) return(mat_pass)
  top_idx <- order(apply(mat_pass, 1, var, na.rm = TRUE), decreasing = TRUE)[seq_len(n_max)]
  mat_pass[top_idx, , drop = FALSE]
}

# --- tab3 helpers (also used by tests) ---------------------------------------

# Filters a DE data frame to rows whose gene_symbol contains gene_search.
# Empty or whitespace-only search returns the full data frame unchanged.
filter_de_table <- function(de_df, gene_search = "") {
  gene_search <- trimws(gene_search)
  if (!nzchar(gene_search)) return(de_df)
  de_df[grepl(gene_search, de_df$gene_symbol, ignore.case = TRUE), ]
}

# Builds a volcano ggplot from a DE data frame using the .data[[col]] pattern.
# Points where neg_log10_padj > padj_thresh are colored with highlight_col;
# all other points use base_col.
make_volcano <- function(de_df, x_col, y_col, base_col, highlight_col, padj_thresh) {
  de_df$highlighted <- de_df$neg_log10_padj > padj_thresh
  ggplot(de_df, aes(x = .data[[x_col]], y = .data[[y_col]],
                    color = highlighted)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(
      values = c("FALSE" = base_col, "TRUE" = highlight_col),
      name   = paste0("-log10(padj) > ", padj_thresh)
    ) +
    labs(x = x_col, y = y_col,
         title = paste("Volcano Plot:", x_col, "vs", y_col)) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# --- tab4 helpers (also used by tests) ---------------------------------------

# Extracts one gene's row from a counts matrix (rownames = gene symbols).
# Returns a named numeric vector; names are sample IDs (column names of the matrix).
# Stops with a clear message if the gene symbol is not present.
get_gene_counts <- function(counts_mat, gene_symbol) {
  if (!(gene_symbol %in% rownames(counts_mat))) {
    stop(paste0("Gene '", gene_symbol, "' not found in the counts matrix."))
  }
  counts_mat[gene_symbol, ]
}

# Joins a named count vector (names = sample IDs) to sample_df and returns a
# tidy data frame with columns: sample_id, normalized_count, and group_col.
pivot_gene_to_long <- function(count_vec, sample_df, group_col) {
  long_df <- data.frame(
    sample_id        = names(count_vec),
    normalized_count = as.numeric(count_vec),
    stringsAsFactors = FALSE
  )
  meta <- sample_df[, c("sample_id", group_col), drop = FALSE]
  merge(long_df, meta, by = "sample_id", all.x = TRUE)
}

# Builds a ggplot for one gene's expression split by group_col.
# Dispatches on plot_type: "box" | "violin" | "bar" | "beeswarm"
make_gene_plot <- function(long_df, gene_name, group_col, plot_type) {
  p <- switch(plot_type,
    box = ggplot(long_df, aes(x = .data[[group_col]], y = normalized_count,
                              fill = .data[[group_col]])) +
          geom_boxplot(alpha = 0.7, outlier.shape = 16),

    violin = ggplot(long_df, aes(x = .data[[group_col]], y = normalized_count,
                                 fill = .data[[group_col]])) +
             geom_violin(trim = FALSE, alpha = 0.7) +
             geom_jitter(width = 0.15, size = 1.2, alpha = 0.5),

    bar = ggplot(long_df, aes(x = .data[[group_col]], y = normalized_count,
                              fill = .data[[group_col]])) +
          stat_summary(fun = mean, geom = "bar", alpha = 0.8) +
          stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2),

    beeswarm = ggplot(long_df, aes(x = .data[[group_col]], y = normalized_count,
                                   color = .data[[group_col]])) +
               geom_beeswarm(size = 2, alpha = 0.8, cex = 2)
  )
  p +
    labs(
      x     = group_col,
      y     = "Normalized Count",
      title = paste("Expression of", gene_name, "by", group_col)
    ) +
    theme_bw() +
    theme(legend.position = "none")
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- navbarPage(
  title = "HD RNA-seq Explorer",
  theme = bs_theme(bootswatch = "flatly"),

  # ── Tab 1: Sample Information ──────────────────────────────────────────────
  tabPanel(
    "Sample Information",
    sidebarLayout(
      sidebarPanel(
        fileInput("tab1_file", "Upload sample_info.csv",
                  accept = ".csv",
                  placeholder = "No file selected"),
        helpText("Upload the sample information CSV exported by the data pipeline."),
        hr(),
        # These selectors are rendered dynamically once a file is loaded
        uiOutput("tab1_col_ui"),
        uiOutput("tab1_group_ui")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Summary",
            br(),
            verbatimTextOutput("tab1_dims"),
            tableOutput("tab1_col_summary")
          ),
          tabPanel(
            "Data Table",
            br(),
            DTOutput("tab1_dt")
          ),
          tabPanel(
            "Distributions",
            br(),
            plotOutput("tab1_plot", height = "450px")
          )
        )
      )
    )
  ),

  # ── Tab 2: Counts Exploration ──────────────────────────────────────────────
  tabPanel(
    "Counts Exploration",
    sidebarLayout(
      sidebarPanel(
        fileInput("tab2_file", "Upload counts_wide.csv",
                  accept = ".csv",
                  placeholder = "No file selected"),
        helpText("Upload the normalized counts matrix (genes × samples). Columns 1–3 are gene metadata and will be stripped automatically."),
        hr(),
        sliderInput("tab2_var_pct",
                    "Minimum variance percentile:",
                    min = 0, max = 100, value = 50, step = 1,
                    post = "%"),
        helpText("Keep genes whose variance is at or above this percentile."),
        sliderInput("tab2_nonzero",
                    "Minimum non-zero samples:",
                    min = 1, max = 69, value = 10, step = 1),
        helpText("Keep genes with at least this many samples having a non-zero count."),
        hr(),
        checkboxInput("tab2_log_heatmap",
                      "Log₂-transform counts for heatmap",
                      value = TRUE),
        hr(),
        radioButtons("tab2_pca_type", "PCA display type:",
                     choices = c("Scatter (PC vs PC)"  = "scatter",
                                 "Beeswarm (top N PCs)" = "beeswarm"),
                     selected = "scatter"),
        # Show PC axis selectors only in scatter mode
        conditionalPanel(
          condition = "input.tab2_pca_type == 'scatter'",
          numericInput("tab2_pc_x", "X-axis PC:", value = 1, min = 1, max = 69),
          numericInput("tab2_pc_y", "Y-axis PC:", value = 2, min = 1, max = 69)
        ),
        # Show top-N selector only in beeswarm mode
        conditionalPanel(
          condition = "input.tab2_pca_type == 'beeswarm'",
          numericInput("tab2_pc_n", "Number of top PCs:", value = 5, min = 2, max = 20)
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Filter Summary",
            br(),
            tableOutput("tab2_filter_summary")
          ),
          tabPanel(
            "Diagnostic Plots",
            br(),
            plotOutput("tab2_plot_var",   height = "350px"),
            br(),
            plotOutput("tab2_plot_zeros", height = "350px")
          ),
          tabPanel(
            "Heatmap",
            br(),
            plotOutput("tab2_heatmap", height = "600px")
          ),
          tabPanel(
            "PCA",
            br(),
            plotOutput("tab2_pca", height = "500px")
          )
        )
      )
    )
  ),

  # ── Tab 3: Differential Expression ────────────────────────────────────────
  tabPanel(
    "Differential Expression",
    sidebarLayout(
      sidebarPanel(
        fileInput("tab3_file", "Upload differential_expression.csv",
                  accept = ".csv",
                  placeholder = "No file selected"),
        helpText("Upload the pre-computed DESeq2 differential expression results (HD vs Control)."),
        hr(),
        h5("Volcano Plot Controls"),
        radioButtons("tab3_x_col", "X-axis column:",
                     choices = c("log2FoldChange", "baseMean", "lfcSE",
                                 "stat", "pvalue", "padj"),
                     selected = "log2FoldChange"),
        radioButtons("tab3_y_col", "Y-axis column:",
                     choices = c("neg_log10_padj", "pvalue", "padj",
                                 "stat", "baseMean", "log2FoldChange"),
                     selected = "neg_log10_padj"),
        colourInput("tab3_base_col", "Base point color:", value = "#999999"),
        colourInput("tab3_highlight_col", "Highlighted point color:", value = "#d7191c"),
        sliderInput("tab3_padj_thresh",
                    "Highlight threshold (−log₁₀ padj ≥):",
                    min = 0, max = 300, value = 60, step = 1),
        helpText("Points where −log₁₀(padj) exceeds this value are drawn in the highlight color."),
        actionButton("tab3_plot_btn", "Update Volcano Plot", class = "btn-primary"),
        hr(),
        h5("Table Controls"),
        textInput("tab3_gene_search", "Search gene symbol:",
                  value = "", placeholder = "e.g. HTT, MAPT"),
        helpText("Filters the table to rows whose gene symbol contains this string (case-insensitive).")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Table",
            br(),
            DTOutput("tab3_dt")
          ),
          tabPanel(
            "Volcano Plot",
            br(),
            plotOutput("tab3_volcano", height = "550px")
          )
        )
      )
    )
  ),

  # ── Tab 4: Individual Gene Expression ─────────────────────────────────────
  tabPanel(
    "Gene Expression",
    sidebarLayout(
      sidebarPanel(
        fileInput("tab4_counts_file", "Upload Normalized Counts Matrix (CSV)",
                  accept = ".csv",
                  placeholder = "No file selected"),
        helpText("Upload counts_wide.csv — columns 1–3 are gene metadata (stripped automatically); remaining columns are one sample each."),
        fileInput("tab4_sample_file", "Upload Sample Information (CSV)",
                  accept = ".csv",
                  placeholder = "No file selected"),
        helpText("Upload sample_info.csv to supply grouping metadata for the plot."),
        hr(),
        uiOutput("tab4_gene_ui"),
        uiOutput("tab4_group_ui"),
        radioButtons("tab4_plot_type", "Plot type:",
                     choices = c("Boxplot"  = "box",
                                 "Violin"   = "violin",
                                 "Bar"      = "bar",
                                 "Beeswarm" = "beeswarm"),
                     selected = "violin"),
        helpText("Choose how to display the distribution of normalized counts across sample groups."),
        actionButton("tab4_plot_btn", "Plot", icon = icon("chart-bar"))
      ),
      mainPanel(
        plotOutput("tab4_plot", height = "500px"),
        br(),
        DTOutput("tab4_gene_table")
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Tab 1 ──────────────────────────────────────────────────────────────────

  tab1_data <- reactive({
    load_csv(input$tab1_file)
  })

  # Numeric column selector — only available after a valid file loads
  output$tab1_col_ui <- renderUI({
    df       <- tab1_data()
    num_cols <- get_numeric_cols(df)
    validate(need(length(num_cols) > 0, "No numeric columns found in this file."))
    selectInput("tab1_col", "Column to plot:", choices = num_cols, selected = num_cols[1])
  })

  # Categorical group-by selector — "None" maps to "" so it is truly optional
  output$tab1_group_ui <- renderUI({
    df       <- tab1_data()
    cat_cols <- get_categorical_cols(df)
    selectInput("tab1_group", "Group by (optional):",
                choices = c("None" = "", cat_cols), selected = "")
  })

  # Summary: row/column counts
  output$tab1_dims <- renderText({
    df <- tab1_data()
    paste0("Rows: ", nrow(df), "    Columns: ", ncol(df))
  })

  # Summary: per-column type + mean/sd or distinct values
  output$tab1_col_summary <- renderTable({
    build_col_summary(tab1_data())
  })

  # Data Table: sortable with column filters, 15 rows per page
  output$tab1_dt <- renderDT({
    datatable(
      tab1_data(),
      options = list(pageLength = 15, scrollX = TRUE),
      filter = "top"
    )
  })

  # Distributions: violin when a group-by is chosen, histogram otherwise.
  # input$tab1_col and input$tab1_group live inside uiOutputs and may be NULL
  # briefly during re-renders — req() holds rendering until they are ready.
  output$tab1_plot <- renderPlot({
    df  <- tab1_data()
    req(input$tab1_col)
    col <- input$tab1_col
    grp <- input$tab1_group  # "" when "None" is selected

    if (!is.null(grp) && nzchar(grp)) {
      ggplot(df, aes(x = .data[[grp]], y = .data[[col]], fill = .data[[grp]])) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
        labs(x = grp, y = col, title = paste(col, "by", grp)) +
        theme_bw() +
        theme(legend.position = "none")
    } else {
      ggplot(df, aes(x = .data[[col]])) +
        geom_histogram(bins = 30, fill = "#4e79a7", color = "white", alpha = 0.85) +
        labs(x = col, y = "Count", title = paste("Distribution of", col)) +
        theme_bw()
    }
  })

  # ── Tab 2 ──────────────────────────────────────────────────────────────────

  # Load counts_wide.csv, strip metadata cols 1-3, return numeric matrix.
  tab2_data <- reactive({
    strip_count_matrix(load_csv(input$tab2_file))
  })

  # Apply both slider filters; returns list(mat = full matrix, pass = logical vector).
  # Keeping the full matrix lets diagnostic plots colour filtered-out genes differently.
  tab2_filtered <- reactive({
    mat  <- tab2_data()
    pass <- pass_var_filter(mat, input$tab2_var_pct) &
            pass_nonzero_filter(mat, input$tab2_nonzero)
    list(mat = mat, pass = pass)
  })

  # Filter Summary ─────────────────────────────────────────────────────────────
  output$tab2_filter_summary <- renderTable({
    res     <- tab2_filtered()
    n_total <- nrow(res$mat)
    n_samp  <- ncol(res$mat)
    n_pass  <- sum(res$pass)
    n_fail  <- n_total - n_pass

    data.frame(
      Metric = c("Number of samples",
                 "Total genes",
                 "Genes passing filter",
                 "Genes not passing filter"),
      Value  = c(
        as.character(n_samp),
        as.character(n_total),
        paste0(n_pass, " (", round(100 * n_pass / n_total, 1), "%)"),
        paste0(n_fail, " (", round(100 * n_fail / n_total, 1), "%)")
      )
    )
  }, striped = TRUE, hover = TRUE)

  # Both diagnostic scatter plots share the same per-gene stats; compute once.
  tab2_gene_stats <- reactive({
    res <- tab2_filtered()
    compute_gene_stats(res$mat, res$pass)
  })

  # Diagnostic Plot 1: median count vs variance
  output$tab2_plot_var <- renderPlot({
    df <- tab2_gene_stats()
    # Sort so Fail (grey) renders behind Pass (blue)
    df <- df[order(df$status), ]
    ggplot(df, aes(x = log_median, y = log_var, color = status)) +
      geom_point(alpha = 0.3, size = 0.6) +
      scale_color_manual(values = c("Fail" = "#d7d7d7", "Pass" = "#2c7bb6"),
                         name = "Filter") +
      labs(x = "log10(median count + 1)",
           y = "log10(variance + 1)",
           title = "Median Count vs Variance") +
      theme_bw() +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  })

  # Diagnostic Plot 2: median count vs number of zero samples
  output$tab2_plot_zeros <- renderPlot({
    df <- tab2_gene_stats()
    df <- df[order(df$status), ]
    ggplot(df, aes(x = log_median, y = n_zeros, color = status)) +
      geom_point(alpha = 0.3, size = 0.6) +
      scale_color_manual(values = c("Fail" = "#d7d7d7", "Pass" = "#d7191c"),
                         name = "Filter") +
      labs(x = "log10(median count + 1)",
           y = "Number of zero samples",
           title = "Median Count vs Zero Samples") +
      theme_bw() +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  })

  # Heatmap ────────────────────────────────────────────────────────────────────
  output$tab2_heatmap <- renderPlot({
    res      <- tab2_filtered()
    mat_pass <- res$mat[res$pass, , drop = FALSE]
    validate(need(nrow(mat_pass) >= 2,
                  "Fewer than 2 genes pass the current filters. Adjust the sliders."))

    # Cap to top 500 genes by variance to keep pheatmap responsive
    mat_pass <- cap_heatmap_genes(mat_pass)

    if (input$tab2_log_heatmap) {
      mat_pass <- log2(mat_pass + 1)
    }

    pheatmap(
      mat_pass,
      show_rownames = FALSE,
      color         = rev(brewer.pal(11, "RdBu")),
      border_color  = NA,
      fontsize_col  = 7
    )
  })

  # PCA ────────────────────────────────────────────────────────────────────────
  output$tab2_pca <- renderPlot({
    res      <- tab2_filtered()
    mat_pass <- res$mat[res$pass, , drop = FALSE]
    validate(need(nrow(mat_pass) >= 2 & ncol(mat_pass) >= 2,
                  "Not enough genes or samples for PCA. Adjust the filters."))

    # prcomp expects samples as rows — transpose so rows = samples, cols = genes
    pca     <- prcomp(t(mat_pass), scale. = TRUE)
    # % variance explained by each PC
    pct_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    n_pcs   <- ncol(pca$x)  # number of available PCs (at most min(genes, samples)-1)

    if (input$tab2_pca_type == "scatter") {
      req(input$tab2_pc_x, input$tab2_pc_y)
      pc_x <- min(as.integer(input$tab2_pc_x), n_pcs)
      pc_y <- min(as.integer(input$tab2_pc_y), n_pcs)
      validate(need(pc_x != pc_y, "X and Y axes must be different PCs."))

      plot_df <- data.frame(
        x = pca$x[, pc_x],
        y = pca$x[, pc_y]
      )
      ggplot(plot_df, aes(x = x, y = y)) +
        geom_point(color = "#4e79a7", size = 2.5, alpha = 0.8) +
        labs(
          x     = sprintf("PC%d (%.1f%% variance)", pc_x, pct_var[pc_x]),
          y     = sprintf("PC%d (%.1f%% variance)", pc_y, pct_var[pc_y]),
          title = "PCA — Filtered Genes"
        ) +
        theme_bw()

    } else {
      # Beeswarm: pivot top N PCs into long format, one point per sample per PC
      req(input$tab2_pc_n)
      n_show  <- min(as.integer(input$tab2_pc_n), n_pcs)
      pc_mat  <- as.data.frame(pca$x[, 1:n_show, drop = FALSE])
      pc_mat$sample <- rownames(pc_mat)

      long_df <- pivot_longer(pc_mat, cols = -sample,
                              names_to = "PC", values_to = "score")

      # Attach % variance label; ordered factor keeps PCs in numeric order on axis
      long_df$PC_label <- factor(
        long_df$PC,
        levels = paste0("PC", seq_len(n_show)),
        labels = sprintf("PC%d\n(%.1f%%)", seq_len(n_show), pct_var[seq_len(n_show)])
      )

      ggplot(long_df, aes(x = PC_label, y = score)) +
        geom_beeswarm(color = "#4e79a7", size = 1.5, alpha = 0.7, cex = 1.5) +
        labs(x = "Principal Component",
             y = "PC Score",
             title = "Top PCs — Beeswarm") +
        theme_bw()
    }
  })

  # ── Tab 3 ──────────────────────────────────────────────────────────────────

  tab3_data <- reactive({
    load_csv(input$tab3_file)
  })

  output$tab3_dt <- renderDT({
    df      <- tab3_data()
    display <- filter_de_table(df, input$tab3_gene_search)
    datatable(
      display,
      options = list(pageLength = 15, scrollX = TRUE),
      filter  = "top"
    )
  })

  output$tab3_volcano <- renderPlot({
    # Re-render only when the button is clicked; isolate() prevents every
    # radio/slider change from immediately triggering a re-plot.
    input$tab3_plot_btn
    isolate({
      df <- tab3_data()
      make_volcano(
        de_df         = df,
        x_col         = input$tab3_x_col,
        y_col         = input$tab3_y_col,
        base_col      = input$tab3_base_col,
        highlight_col = input$tab3_highlight_col,
        padj_thresh   = input$tab3_padj_thresh
      )
    })
  })

  # ── Tab 4 ──────────────────────────────────────────────────────────────────

  # Load counts_wide.csv; strip gene metadata cols 1-3; store gene_symbol as rownames
  tab4_counts <- reactive({
    df  <- load_csv(input$tab4_counts_file)
    mat <- as.matrix(df[, -(1:3)])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- df$gene_symbol
    mat
  })

  tab4_samples <- reactive({
    load_csv(input$tab4_sample_file)
  })

  # Gene selector — searchable selectize, populated once the counts file loads
  output$tab4_gene_ui <- renderUI({
    mat   <- tab4_counts()
    genes <- rownames(mat)
    selectizeInput("tab4_gene", "Select gene:",
                   choices  = genes,
                   selected = genes[1],
                   options  = list(placeholder = "Search gene symbol..."))
  })

  # Group selector — categorical cols from sample_info, defaulting to 'diagnosis'
  output$tab4_group_ui <- renderUI({
    df       <- tab4_samples()
    cat_cols <- get_categorical_cols(df)
    validate(need(length(cat_cols) > 0, "No categorical columns found in sample info."))
    default_col <- if ("diagnosis" %in% cat_cols) "diagnosis" else cat_cols[1]
    selectInput("tab4_group", "Group by:",
                choices  = cat_cols,
                selected = default_col)
  })

  # Subset to one gene row, pivot to long format, join sample metadata.
  # req() ensures both file reactives and both UI-driven inputs are ready.
  tab4_long <- reactive({
    req(input$tab4_gene, input$tab4_group)
    count_vec <- get_gene_counts(tab4_counts(), input$tab4_gene)
    pivot_gene_to_long(count_vec, tab4_samples(), input$tab4_group)
  })

  # Plot re-renders only on button click; bindEvent holds rendering until clicked
  output$tab4_plot <- renderPlot({
    long_df <- tab4_long()
    make_gene_plot(long_df, input$tab4_gene, input$tab4_group, input$tab4_plot_type)
  }) |> bindEvent(input$tab4_plot_btn)

  # Table updates reactively whenever gene or group selection changes
  output$tab4_gene_table <- renderDT({
    long_df <- tab4_long()
    datatable(
      long_df,
      options = list(pageLength = 15, scrollX = TRUE),
      filter  = "top"
    )
  })
}

shinyApp(ui, server)
