###############################################################################
# prepare_data.R -- run ONCE on your own machine, then commit data/app_data.rds
#
# Turns raw_counts_table.txt into a single small .rds that the Shiny app reads.
# Everything that needs Bioconductor (org.Hs.eg.db) happens HERE, so the app
# itself only needs CRAN packages. That is what makes it deployable anywhere.
#
# Normalisation follows hiGN_gene_check_from_raw.R exactly:
#   log2(CPM + 1) on the UNFILTERED matrix, so a gene that is off keeps a
#   value of 0 and still appears -- it is never silently dropped.
#
# Usage:  Rscript prepare_data.R
###############################################################################

COUNTS_FILE <- "raw_counts_table.txt"
OUT_FILE    <- "data/app_data.rds"

# Detection rule -- identical to the analysis script
MIN_COUNT <- 10
MIN_CPM   <- 1
MIN_REPS  <- 3     # DIV28 has n = 3, so 3 means "all replicates" there

# Reference genes shown alongside any gene the user picks, so a result can be
# read against a known positive and a known negative on the same axis.
REF_POS <- c("ACTB", "GAPDH")
REF_NEG <- c("CD19", "INS")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

dir.create("data", showWarnings = FALSE)


## -- 1. counts ---------------------------------------------------------------

cts <- read.delim(COUNTS_FILE, row.names = 1, check.names = FALSE) |> as.matrix()
stopifnot(all(cts == round(cts)), all(cts >= 0), !any(is.na(cts)))
mode(cts) <- "integer"
message(sprintf("Raw counts: %d genes x %d samples", nrow(cts), ncol(cts)))

rownames(cts) <- sub("\\..*$", "", rownames(cts))          # drop version suffix
if (anyDuplicated(rownames(cts))) {
  n_dup <- sum(duplicated(rownames(cts)))
  cts <- cts[order(rowSums(cts), decreasing = TRUE), , drop = FALSE]
  cts <- cts[!duplicated(rownames(cts)), , drop = FALSE]
  message(sprintf("Collapsed %d duplicate Ensembl IDs (kept highest total)", n_dup))
}

# Library sizes from the FULL matrix, before any filtering. If you compute
# these after subsetting, every CPM in the app is wrong.
lib_size <- colSums(cts)


## -- 2. sample metadata ------------------------------------------------------

coldata <- data.frame(row.names = colnames(cts), sample = colnames(cts)) |>
  mutate(div       = as.integer(sub("^DIV(\\d+)_.*$", "\\1", sample)),
         replicate = sub("^.*_(\\d+)$", "\\1", sample)) |>
  mutate(timepoint = factor(paste0("DIV", div),
                            levels = paste0("DIV", sort(unique(div)))))
stopifnot(identical(rownames(coldata), colnames(cts)))
print(dplyr::count(coldata, timepoint))                    # expect 4 / 4 / 3 / 4


## -- 3. Ensembl -> symbol for EVERY gene -------------------------------------
# Done here so lab members can type "TLR3" instead of an Ensembl ID.

sym <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = rownames(cts), keytype = "ENSEMBL",
  column = "SYMBOL", multiVals = "first"))

gene_map <- tibble(ensembl = rownames(cts),
                   symbol  = unname(sym),
                   total   = rowSums(cts)) |>
  mutate(symbol = ifelse(is.na(symbol) | symbol == "", ensembl, symbol))

# One symbol can map to several Ensembl IDs. Mark the one with the most reads
# as primary so it sorts first in the dropdown; keep the others searchable.
# Done by sorting rather than group_by(): there are ~40,000 groups here, and
# a grouped mutate over that many is slow enough to look like a hang.
gene_map <- gene_map[order(gene_map$symbol, -gene_map$total), ]
gene_map$primary <- !duplicated(gene_map$symbol)
gene_map$label   <- ifelse(gene_map$symbol == gene_map$ensembl,
                           gene_map$ensembl,
                           paste0(gene_map$symbol, "  \u00b7  ", gene_map$ensembl))

stopifnot(all(c("ensembl", "symbol", "total", "primary", "label") %in% names(gene_map)))

message(sprintf("Mapped %d / %d genes to a symbol",
                sum(gene_map$symbol != gene_map$ensembl), nrow(gene_map)))


## -- 4. background distribution for the calibration plot ---------------------
# Precomputed so the app never has to hold the full log-CPM matrix in memory.

cpm    <- t(t(cts) / lib_size) * 1e6
logcpm <- log2(cpm + 1)

bg_density <- lapply(levels(coldata$timepoint), function(tp) {
  v <- as.vector(logcpm[, coldata$sample[coldata$timepoint == tp], drop = FALSE])
  d <- density(v, n = 512)
  data.frame(timepoint = tp, x = d$x, y = d$y)
}) |> bind_rows() |>
  mutate(timepoint = factor(timepoint, levels = levels(coldata$timepoint)))

rm(cpm, logcpm)


## -- 5. resolve reference genes ---------------------------------------------

resolve_ref <- function(symbols) {
  hit <- gene_map |> filter(symbol %in% symbols, primary)
  missing <- setdiff(symbols, hit$symbol)
  if (length(missing))
    warning("Reference gene not found: ", paste(missing, collapse = ", "))
  setNames(hit$ensembl, hit$symbol)
}

reference <- list(positive = resolve_ref(REF_POS),
                  negative = resolve_ref(REF_NEG))
print(reference)


## -- 6. save -----------------------------------------------------------------

saveRDS(list(
  counts     = cts,
  lib_size   = lib_size,
  coldata    = coldata,
  gene_map   = as.data.frame(gene_map),
  bg_density = bg_density,
  reference  = reference,
  params     = list(MIN_COUNT = MIN_COUNT, MIN_CPM = MIN_CPM, MIN_REPS = MIN_REPS),
  built_on   = as.character(Sys.Date())
), OUT_FILE, compress = "xz")

message(sprintf("Wrote %s (%.2f MB)", OUT_FILE, file.size(OUT_FILE) / 1e6))
