###############################################################################
# app.R -- hiGN expression browser
#
# Answers one question for a lab member: is my gene expressed in hiPSC-derived
# glutamatergic neurons, and how strongly, at DIV0 / DIV14 / DIV28 / DIV50?
#
# Reads data/app_data.rds, built by prepare_data.R. CRAN packages only.
###############################################################################

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

DATA_PATH <- "data/app_data.rds"
if (!file.exists(DATA_PATH))
  stop("data/app_data.rds not found. Run prepare_data.R once to build it.")

D         <- readRDS(DATA_PATH)
cts       <- D$counts
lib_size  <- D$lib_size
coldata   <- D$coldata
gene_map  <- D$gene_map
bg        <- D$bg_density
reference <- D$reference
P         <- D$params

TP_LEVELS    <- levels(coldata$timepoint)
MAX_GENES    <- 12
GENE_CHOICES <- setNames(gene_map$ensembl, gene_map$label)

ROLE_COLOURS <- c("Selected gene"      = "#1f6f8b",
                  "Positive control"   = "#2e7d32",
                  "Negative control"   = "#9e9e9e")

theme_set(theme_bw(base_size = 12))


## ---------------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------------

# Per-sample values for a set of Ensembl IDs. CPM uses the full-matrix library
# sizes stored at prep time, so these numbers match the analysis script.
get_long <- function(ens, roles) {
  m   <- cts[ens, , drop = FALSE]
  cpm <- t(t(m) / lib_size) * 1e6
  data.frame(
    ensembl = rep(rownames(m), times = ncol(m)),
    sample  = rep(colnames(m), each  = nrow(m)),
    count   = as.vector(m),
    cpm     = as.vector(cpm),
    stringsAsFactors = FALSE
  ) |>
    mutate(logcpm = log2(cpm + 1),
           role   = factor(unname(roles[ensembl]), levels = names(ROLE_COLOURS))) |>
    left_join(gene_map[, c("ensembl", "symbol", "label")], by = "ensembl") |>
    left_join(coldata, by = "sample") |>
    mutate(symbol = factor(symbol, levels = unique(gene_map$symbol[
      match(ens, gene_map$ensembl)])))
}

# Detection call per gene per timepoint: raw reads AND CPM must both clear
# their thresholds in at least MIN_REPS replicates. The CPM half matters
# because the DIV50 libraries are shallow (9.5-16M vs 18-39M elsewhere).
make_calls <- function(long) {
  long |>
    group_by(ensembl, symbol, role, timepoint) |>
    summarise(n           = n(),
              n_pass      = sum(count >= P$MIN_COUNT & cpm >= P$MIN_CPM),
              mean_count  = mean(count),
              mean_cpm    = mean(cpm),
              mean_logcpm = mean(logcpm),
              .groups     = "drop") |>
    mutate(call = ifelse(n_pass >= pmin(P$MIN_REPS, n), "Expressed", "Not detected"),
           timepoint = factor(as.character(timepoint), levels = TP_LEVELS)) |>
    arrange(symbol, timepoint)
}

# One plain sentence per gene -- the thing most people actually came for.
verdict_text <- function(cl) {
  yes  <- as.character(cl$timepoint[cl$call == "Expressed"])
  peak <- cl[which.max(cl$mean_cpm), ]
  if (sum(cl$mean_count) == 0)
    return("Zero reads in all 15 samples. The gene is in the annotation, so it was measured -- nothing mapped to it.")
  if (!length(yes))
    return(sprintf("Not detected at any timepoint. Highest mean was %.2f CPM at %s, below the %g CPM threshold.",
                   peak$mean_cpm, peak$timepoint, P$MIN_CPM))
  if (length(yes) == length(TP_LEVELS))
    return(sprintf("Expressed at every timepoint. Highest at %s (%.1f CPM).",
                   peak$timepoint, peak$mean_cpm))
  sprintf("Expressed at %s, not at %s. Highest at %s (%.1f CPM).",
          paste(yes, collapse = ", "),
          paste(setdiff(TP_LEVELS, yes), collapse = ", "),
          peak$timepoint, peak$mean_cpm)
}

verdict_card <- function(cl) {
  n_yes <- sum(cl$call == "Expressed")
  accent <- if (n_yes == 0) "#9e9e9e" else "#2e7d32"
  chips <- lapply(TP_LEVELS, function(tp) {
    row <- cl[cl$timepoint == tp, ]
    on  <- nrow(row) && row$call == "Expressed"
    tags$span(class = if (on) "chip chip-on" else "chip chip-off",
              tp, tags$b(sprintf("%.1f", if (nrow(row)) row$mean_cpm else NA)), "CPM")
  })
  div(class = "card", style = paste0("border-left-color:", accent),
      div(class = "card-gene", as.character(cl$symbol[1]),
          tags$span(class = "card-ens", cl$ensembl[1])),
      div(class = "card-verdict", verdict_text(cl)),
      div(class = "chips", chips))
}


## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

app_css <- HTML("
  body { background:#fbfbfa; }
  .card { background:#fff; border:1px solid #e6e4e0; border-left:4px solid #ccc;
          border-radius:4px; padding:14px 16px; margin-bottom:12px; }
  .card-gene { font-size:17px; font-weight:600; letter-spacing:.01em; }
  .card-ens { font-size:12px; font-weight:400; color:#8a857e; margin-left:8px;
              font-family:ui-monospace,Menlo,Consolas,monospace; }
  .card-verdict { margin:6px 0 10px; color:#3d3a35; }
  .chips { display:flex; flex-wrap:wrap; gap:6px; }
  .chip { font-size:12px; padding:3px 9px; border-radius:20px;
          border:1px solid transparent; }
  .chip-on  { background:#e8f3e9; color:#1f5c23; border-color:#c6e2c8; }
  .chip-off { background:#f2f1ef; color:#77726b; border-color:#e2e0dc; }
  .empty { color:#8a857e; padding:40px 10px; text-align:center; }
  .note { font-size:12px; color:#8a857e; }
  .well { background:#fff; border-color:#e6e4e0; }
")

ui <- fluidPage(
  tags$head(tags$style(app_css)),
  titlePanel("hiGN expression browser"),
  tags$p(class = "note",
         "hiPSC-derived glutamatergic neurons, bulk RNA-seq across differentiation.",
         "15 libraries: DIV0 (n=4), DIV14 (n=4), DIV28 (n=3), DIV50 (n=4)."),
  tags$hr(),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput("genes", "Genes", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Type a symbol, e.g. TLR3",
                                    maxItems = MAX_GENES)),
      helpText(class = "note",
               "Search by HGNC symbol or Ensembl ID. Up to", MAX_GENES, "at a time."),
      checkboxInput("show_ref", "Show control genes for scale", TRUE),
      radioButtons("yaxis", "Plot y-axis",
                   c("CPM (log scale)"    = "cpm",
                     "Raw reads (log scale)" = "count",
                     "log2(CPM + 1)"      = "logcpm")),
      tags$hr(),
      downloadButton("dl_calls",   "Download summary (CSV)", class = "btn-sm"),
      tags$br(), tags$br(),
      downloadButton("dl_samples", "Download per-sample values (CSV)", class = "btn-sm"),
      tags$br(), tags$br(),
      downloadButton("dl_plot",    "Download plot (PDF)", class = "btn-sm")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("Result",
                 br(),
                 uiOutput("cards"),
                 conditionalPanel("input.genes.length > 0",
                                  h4("Per timepoint"), DTOutput("calls_tbl"))),
        tabPanel("Expression plot", br(), plotOutput("plot", height = "620px")),
        tabPanel("Per-sample values", br(), DTOutput("samples_tbl")),
        tabPanel("Compare to all genes", br(),
                 plotOutput("bg_plot", height = "560px"),
                 tags$p(class = "note",
                        "Grey curve is the distribution of all", format(nrow(cts), big.mark = ","),
                        "genes in that timepoint. A line sitting in the left-hand spike is a gene",
                        "behaving like the silent majority of the transcriptome.")),
        tabPanel("How this works", br(), uiOutput("about"))
      )
    )
  )
)


## ---------------------------------------------------------------------------
## server
## ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # Server-side selectize: 58k choices would freeze the browser client-side.
  updateSelectizeInput(session, "genes", choices = GENE_CHOICES, server = TRUE)

  roles <- reactive({
    sel <- input$genes
    r <- setNames(rep("Selected gene", length(sel)), sel)
    if (isTRUE(input$show_ref)) {
      r <- c(r,
             setNames(rep("Positive control", length(reference$positive)), reference$positive),
             setNames(rep("Negative control", length(reference$negative)), reference$negative))
      r <- r[!duplicated(names(r))]
    }
    r
  })

  long <- reactive({
    req(length(input$genes) > 0)
    get_long(names(roles()), roles())
  })

  calls <- reactive(make_calls(long()))

  selected_calls <- reactive(filter(calls(), role == "Selected gene"))

  output$cards <- renderUI({
    if (!length(input$genes))
      return(div(class = "empty",
                 h4("Pick a gene to start"),
                 tags$p("Type a gene symbol in the box on the left.",
                        "You will get a yes/no call at each timepoint, the numbers behind it,",
                        "and a plot you can download.")))
    cl <- selected_calls()
    lapply(split(cl, droplevels(cl$symbol)), verdict_card)
  })

  output$calls_tbl <- renderDT({
    calls() |>
      transmute(Gene = symbol, Ensembl = ensembl, Role = role, Timepoint = timepoint,
                `Mean reads` = round(mean_count, 1),
                `Mean CPM`   = round(mean_cpm, 2),
                `log2(CPM+1)`= round(mean_logcpm, 2),
                `Replicates passing` = paste0(n_pass, " / ", n),
                Call = call) |>
      datatable(rownames = FALSE, options = list(pageLength = 12, dom = "tip")) |>
      formatStyle("Call", color = styleEqual(c("Expressed", "Not detected"),
                                             c("#1f5c23", "#8a857e")),
                  fontWeight = styleEqual("Expressed", "600", default = "400"))
  })

  output$samples_tbl <- renderDT({
    long() |>
      transmute(Gene = symbol, Sample = sample, Timepoint = timepoint,
                Reads = count, CPM = round(cpm, 2), `log2(CPM+1)` = round(logcpm, 2)) |>
      datatable(rownames = FALSE, options = list(pageLength = 15, dom = "tip"))
  })

  build_plot <- reactive({
    d <- long()
    spec <- switch(input$yaxis,
      cpm    = list(y = d$cpm + 0.01, lab = "CPM + 0.01", log = TRUE,  hline = P$MIN_CPM),
      count  = list(y = d$count + 1,  lab = "reads + 1",  log = TRUE,  hline = P$MIN_COUNT),
      logcpm = list(y = d$logcpm,     lab = "log2(CPM + 1)", log = FALSE, hline = NA))
    d$yval <- spec$y

    p <- ggplot(d, aes(timepoint, yval, group = symbol, colour = role)) +
      stat_summary(fun = mean, geom = "line", linewidth = 0.7) +
      geom_point(size = 1.5, alpha = 0.75) +
      facet_wrap(~ symbol) +
      scale_colour_manual(values = ROLE_COLOURS, drop = TRUE) +
      labs(x = NULL, y = spec$lab, colour = NULL,
           title = "Expression across differentiation",
           subtitle = if (!is.na(spec$hline))
             sprintf("Dashed line = detection threshold (%g). Points are replicates; line is the mean.",
                     spec$hline)
           else "Points are replicates; line is the mean.") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            strip.background = element_rect(fill = "#f2f1ef", colour = "#e6e4e0"))

    if (!is.na(spec$hline))
      p <- p + geom_hline(yintercept = spec$hline, linetype = "dashed", colour = "grey40")
    if (spec$log) p <- p + scale_y_log10()
    p
  })

  output$plot <- renderPlot(build_plot())

  output$bg_plot <- renderPlot({
    cl <- calls()
    ggplot(bg, aes(x, y)) +
      geom_area(fill = "#ececea", colour = NA) +
      geom_vline(data = cl, aes(xintercept = mean_logcpm, colour = role),
                 linewidth = 0.5) +
      facet_grid(timepoint ~ .) +
      scale_colour_manual(values = ROLE_COLOURS, drop = TRUE) +
      labs(title = "Where your genes sit in the whole transcriptome",
           x = "log2(CPM + 1)", y = "density", colour = NULL) +
      theme(legend.position = "bottom")
  })

  output$about <- renderUI({
    tagList(
      h4("What the call means"),
      tags$p(sprintf(
        "A gene is called Expressed at a timepoint when at least %d replicates have both \u2265 %d raw reads and \u2265 %g CPM. Both halves matter: the DIV50 libraries are shallower (9.5-16 M reads) than the rest (18-39 M), so a reads-only cutoff would under-call expression at exactly the timepoint of most interest, and a CPM-only cutoff would let noise through in the deep libraries.",
        P$MIN_REPS, P$MIN_COUNT, P$MIN_CPM)),
      h4("Normalisation"),
      tags$p("Counts are converted to CPM using each library's total mapped reads, then log2(CPM + 1) where a log scale is shown. This is computed on the unfiltered matrix, so a gene with zero reads keeps a value of 0 and still appears. That distinction matters: a missing row reads as \"not measured\", a zero row reads as \"not expressed\"."),
      tags$p("DESeq2 median-of-ratios size factors were checked against library-size factors on this dataset and agree to within 10% on every sample, so the simpler normalisation does not change any conclusion here."),
      h4("What this app cannot tell you"),
      tags$ul(
        tags$li("Whether a difference between timepoints is statistically significant. This is a detection and magnitude view, not a differential expression test."),
        tags$li("Protein-level expression."),
        tags$li("Which cells in the culture express the gene. Bulk RNA-seq averages over everything in the well.")),
      h4("Duplicate identifiers"),
      tags$p("Version suffixes were stripped from Ensembl IDs and duplicates collapsed to the ID with the most reads. Where one symbol maps to several Ensembl IDs, all remain searchable -- check the ID under the gene name if a result looks unexpected."),
      tags$hr(),
      tags$p(class = "note", "Data built on ", D$built_on, " from raw_counts_table.txt.")
    )
  })

  ## downloads
  output$dl_calls <- downloadHandler(
    filename = function() paste0("expression_summary_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(calls(), f, row.names = FALSE))

  output$dl_samples <- downloadHandler(
    filename = function() paste0("per_sample_values_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(long(), f, row.names = FALSE))

  output$dl_plot <- downloadHandler(
    filename = function() paste0("expression_plot_", Sys.Date(), ".pdf"),
    content  = function(f) ggsave(f, build_plot(), width = 11, height = 8))
}

shinyApp(ui, server)
