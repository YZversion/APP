library(shiny)
library(bslib)
library(ggplot2)
library(DT)

# ── Shared helper ──────────────────────────────────────────────────────────────
# Reads and validates a CSV uploaded via fileInput.
# Uses validate() so errors appear as user-facing messages, not crashes.
# Reusable across all tabs — call from inside a reactive({}).
load_csv <- function(file_input) {
  validate(need(!is.null(file_input), "Upload a CSV file to begin."))
  validate(need(
    tolower(tools::file_ext(file_input$name)) == "csv",
    paste0("Expected a .csv file; got .", tools::file_ext(file_input$name), ".")
  ))
  df <- tryCatch(
    read.csv(file_input$datapath, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  validate(need(!is.null(df), "Could not parse the file. Check that it is a valid CSV."))
  validate(need(nrow(df) > 1, "File must have more than 1 data row."))
  validate(need(ncol(df) > 1, "File must have more than 1 column."))
  df
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

  # ── Tab 2 ──────────────────────────────────────────────────────────────────
  tabPanel(
    "Counts Exploration",
    fluidPage(
      textOutput("placeholder_tab2")
    )
  ),

  # ── Tab 3 ──────────────────────────────────────────────────────────────────
  tabPanel(
    "Differential Expression",
    fluidPage(
      textOutput("placeholder_tab3")
    )
  ),

  # ── Tab 4 (choose-your-own — not yet decided) ──────────────────────────────
  tabPanel(
    "Tab 4",
    fluidPage(
      textOutput("placeholder_tab4")
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
    df <- tab1_data()
    num_cols <- names(df)[sapply(df, is.numeric)]
    validate(need(length(num_cols) > 0, "No numeric columns found in this file."))
    selectInput("tab1_col", "Column to plot:", choices = num_cols, selected = num_cols[1])
  })

  # Categorical group-by selector — "None" maps to "" so it is truly optional
  output$tab1_group_ui <- renderUI({
    df <- tab1_data()
    cat_cols <- names(df)[!sapply(df, is.numeric)]
    selectInput("tab1_group", "Group by (optional):",
                choices = c("None" = "", cat_cols), selected = "")
  })

  # Summary: row/column counts
  output$tab1_dims <- renderText({
    df <- tab1_data()
    paste0("Rows: ", nrow(df), "    Columns: ", ncol(df))
  })

  # Summary: one row per column — numeric cols show mean (sd), others show distinct values
  output$tab1_col_summary <- renderTable({
    df <- tab1_data()
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

  # ── Tabs 2–4 placeholders ──────────────────────────────────────────────────
  output$placeholder_tab2 <- renderText("Tab 2 — Counts Exploration (not yet implemented)")
  output$placeholder_tab3 <- renderText("Tab 3 — Differential Expression (not yet implemented)")
  output$placeholder_tab4 <- renderText("Tab 4 — TBD (not yet implemented)")
}

shinyApp(ui, server)
