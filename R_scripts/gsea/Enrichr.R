#!/usr/bin/env Rscript
# Enrichr.R — online gene-set enrichment via the enrichR package (Maayan lab Enrichr).
#
#   Rscript R_scripts/gsea/Enrichr.R --genes genes.txt --treatment demo
#   Rscript R_scripts/gsea/Enrichr.R --deseq2 merged_DESeq2.csv --treatment X --padj 0.1
#
# Requires: install.packages("enrichR")

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- list(
    genes = NULL, deseq2 = NULL, te_group = NULL,
    treatment = Sys.getenv("RIBO_SEQ_TREATMENT", unset = "run"),
    padj = 0.1, site = "Enrichr", outdir = NULL, plot = TRUE, help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      opts$help <- TRUE
      return(opts)
    }
    val <- if (i < length(args)) args[[i + 1L]] else NA_character_
    if (key == "--genes") opts$genes <- val
    else if (key == "--deseq2") opts$deseq2 <- val
    else if (key == "--te-group") opts$te_group <- val
    else if (key == "--treatment") opts$treatment <- val
    else if (key == "--padj") opts$padj <- as.numeric(val)
    else if (key == "--site") opts$site <- val
    else if (key == "--outdir") opts$outdir <- val
    else if (key == "--no-plot") {
      opts$plot <- FALSE
      i <- i + 1L
      next
    } else stop("Unknown argument: ", key, call. = FALSE)
    i <- i + 2L
  }
  opts
}

print_help <- function() {
  cat("Enrichr.R — Enrichr online enrichment\n")
  cat("  --genes FILE / --deseq2 FILE [--te-group STR]\n")
  cat("  --treatment STR  --padj NUM  --site STR  --outdir DIR  --no-plot\n")
  cat("Requires: install.packages(\"enrichR\")\n")
}

source_common <- function() {
  for (p in c("common_variables.R", "R_scripts/common_variables.R")) {
    if (file.exists(p)) {
      source(p)
      return(invisible(TRUE))
    }
  }
  parent_dir <<- Sys.getenv("RIBO_SEQ_PARENT_DIR")
  if (!nzchar(parent_dir)) parent_dir <<- "results"
  invisible(TRUE)
}

load_genes <- function(opts) {
  if (!is.null(opts$genes)) {
    genes <- trimws(readLines(opts$genes, warn = FALSE))
    genes <- genes[nzchar(genes) & !startsWith(genes, "#")]
    return(list(default = unique(genes)))
  }
  if (!is.null(opts$deseq2)) {
    df <- utils::read.csv(opts$deseq2, stringsAsFactors = FALSE)
    stopifnot("gene_sym" %in% names(df))
    if (!is.null(opts$te_group)) {
      stopifnot("TE_group" %in% names(df))
      g <- unique(df$gene_sym[df$TE_group == opts$te_group])
      g <- g[!is.na(g) & nzchar(as.character(g))]
      return(setNames(list(g), gsub("[^A-Za-z0-9]+", "_", opts$te_group)))
    }
    out <- list()
    if ("TE_group" %in% names(df)) {
      for (grp in c("TE up", "TE down")) {
        g <- unique(df$gene_sym[df$TE_group == grp])
        g <- g[!is.na(g) & nzchar(as.character(g))]
        if (length(g) > 0) out[[gsub("[^A-Za-z0-9]+", "_", grp)]] <- g
      }
    }
    if (!length(out)) {
      g <- unique(df$gene_sym)
      out$all_genes <- g[!is.na(g) & nzchar(as.character(g))]
    }
    return(out)
  }
  stop("Provide --genes or --deseq2", call. = FALSE)
}

DEFAULT_DBS <- c(
  "GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
  "KEGG_2021_Human", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures",
  "TRRUST_Transcription_Factors_2019", "Reactome_2022", "TargetScan_microRNA_2017"
)

main <- function() {
  opts <- parse_args()
  if (isTRUE(opts$help)) {
    print_help()
    quit(save = "no", status = 0)
  }
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop("Package 'enrichR' is required. install.packages(\"enrichR\")", call. = FALSE)
  }
  suppressPackageStartupMessages(library(enrichR))

  source_common()
  outdir <- if (!is.null(opts$outdir) && nzchar(opts$outdir)) opts$outdir else file.path(parent_dir, "Analysis", "Enrichr")
  plotdir <- file.path(parent_dir, "plots", "Enrichr")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (opts$plot) dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

  message("Connecting to Enrichr site: ", opts$site)
  setEnrichrSite(opts$site)

  gene_sets <- load_genes(opts)
  dbs_env <- Sys.getenv("RIBO_SEQ_ENRICHR_DBS")
  dbs_selected <- if (nzchar(dbs_env)) trimws(strsplit(dbs_env, ",", fixed = TRUE)[[1]]) else DEFAULT_DBS

  for (label in names(gene_sets)) {
    genes <- gene_sets[[label]]
    message("Enrichr query [", label, "]: ", length(genes), " genes")
    if (length(genes) < 3) {
      warning("Skipping ", label, ": fewer than 3 genes")
      next
    }
    enriched <- enrichr(genes, dbs_selected)
    prefix <- paste(opts$treatment, label, sep = "_")
    for (db in names(enriched)) {
      tab <- enriched[[db]]
      if (is.null(tab) || !nrow(tab)) next
      out_csv <- file.path(outdir, paste0(prefix, "_", db, ".csv"))
      utils::write.csv(tab, out_csv, row.names = FALSE)
      message("  wrote ", out_csv)
    }
    if (opts$plot && requireNamespace("ggplot2", quietly = TRUE) &&
        !is.null(enriched[["MSigDB_Hallmark_2020"]]) && nrow(enriched$MSigDB_Hallmark_2020) > 0) {
      hall <- enriched$MSigDB_Hallmark_2020
      hall <- hall[hall$Adjusted.P.value < opts$padj, , drop = FALSE]
      if (nrow(hall) > 0) {
        hall <- hall[order(hall$Odds.Ratio), , drop = FALSE]
        hall$Term <- factor(hall$Term, levels = hall$Term)
        p <- ggplot2::ggplot(hall, ggplot2::aes(x = Term, y = Odds.Ratio, fill = -log10(Adjusted.P.value))) +
          ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_classic() +
          ggplot2::labs(title = paste(opts$treatment, label, "Hallmark"), y = "Odds ratio", x = NULL)
        png_path <- file.path(plotdir, paste0(prefix, "_Hallmark.png"))
        ggplot2::ggsave(png_path, p, width = 8, height = max(4, 0.25 * nrow(hall)))
        message("  wrote ", png_path)
      }
    }
  }
  message("Enrichr done. Outputs under ", outdir)
}

if (sys.nframe() == 0) main()
