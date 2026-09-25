#!/usr/bin/env Rscript
# gProfiler.R — g:Profiler g:GOSt via gprofiler2
# Docs: https://biit.cs.ut.ee/gprofiler/page/r
#
#   Rscript R_scripts/gsea/gProfiler.R --genes genes.txt --organism hsapiens
#   Rscript R_scripts/gsea/gProfiler.R --deseq2 merged.csv --te-group "TE up" --organism hsapiens
#
# Requires: install.packages("gprofiler2")

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- list(
    genes = NULL, deseq2 = NULL, te_group = NULL,
    treatment = Sys.getenv("RIBO_SEQ_TREATMENT", unset = "run"),
    organism = Sys.getenv("RIBO_SEQ_ORGANISM", unset = "hsapiens"),
    sources = NULL, user_threshold = 0.05, correction = "g_SCS",
    outdir = NULL, plot = TRUE, help = FALSE
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
    else if (key == "--organism") opts$organism <- val
    else if (key == "--sources") opts$sources <- val
    else if (key == "--user-threshold") opts$user_threshold <- as.numeric(val)
    else if (key == "--correction") opts$correction <- val
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
  cat("gProfiler.R — g:Profiler g:GOSt (gprofiler2)\n")
  cat("  --genes FILE / --deseq2 FILE [--te-group STR]\n")
  cat("  --organism ID  --sources CSV  --user-threshold NUM  --correction METHOD\n")
  cat("  --treatment STR  --outdir DIR  --no-plot\n")
  cat("Requires: install.packages(\"gprofiler2\")\n")
  cat("Docs: https://biit.cs.ut.ee/gprofiler/page/r\n")
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

load_gene_lists <- function(opts) {
  if (!is.null(opts$genes)) {
    genes <- trimws(readLines(opts$genes, warn = FALSE))
    genes <- genes[nzchar(genes) & !startsWith(genes, "#")]
    return(list(default = unique(genes)))
  }
  if (!is.null(opts$deseq2)) {
    df <- utils::read.csv(opts$deseq2, stringsAsFactors = FALSE)
    stopifnot("gene_sym" %in% names(df))
    out <- list()
    if (!is.null(opts$te_group)) {
      stopifnot("TE_group" %in% names(df))
      g <- unique(df$gene_sym[df$TE_group == opts$te_group])
      g <- g[!is.na(g) & nzchar(as.character(g))]
      out[[gsub("[^A-Za-z0-9]+", "_", opts$te_group)]] <- g
      return(out)
    }
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

main <- function() {
  opts <- parse_args()
  if (isTRUE(opts$help)) {
    print_help()
    quit(save = "no", status = 0)
  }
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package 'gprofiler2' is required. install.packages(\"gprofiler2\")", call. = FALSE)
  }
  suppressPackageStartupMessages(library(gprofiler2))
  source_common()

  outdir <- if (!is.null(opts$outdir) && nzchar(opts$outdir)) opts$outdir else file.path(parent_dir, "Analysis", "gProfiler")
  plotdir <- file.path(parent_dir, "plots", "gProfiler")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (opts$plot) dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

  sources <- if (!is.null(opts$sources) && nzchar(opts$sources)) {
    trimws(strsplit(opts$sources, ",", fixed = TRUE)[[1]])
  } else {
    c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP")
  }

  gene_sets <- load_gene_lists(opts)
  for (label in names(gene_sets)) {
    genes <- as.character(gene_sets[[label]])
    message("g:Profiler query [", label, "]: ", length(genes), " genes; organism=", opts$organism)
    if (length(genes) < 3) {
      warning("Skipping ", label, ": fewer than 3 genes")
      next
    }
    gostres <- gost(
      query = genes, organism = opts$organism, sources = sources,
      user_threshold = opts$user_threshold, correction_method = opts$correction,
      significant = TRUE, evcodes = FALSE
    )
    prefix <- paste(opts$treatment, label, sep = "_")
    if (is.null(gostres) || is.null(gostres$result) || !nrow(gostres$result)) {
      message("  no significant terms for ", label)
      next
    }
    res <- gostres$result
    keep <- vapply(res, function(col) !is.list(col), logical(1))
    csv_path <- file.path(outdir, paste0(prefix, "_gost.csv"))
    utils::write.csv(res[, keep, drop = FALSE], csv_path, row.names = FALSE)
    message("  wrote ", csv_path)
    if (opts$plot) {
      png_path <- file.path(plotdir, paste0(prefix, "_gostplot.png"))
      tryCatch({
        p <- gostplot(gostres, capped = TRUE, interactive = FALSE)
        grDevices::png(png_path, width = 1000, height = 600)
        print(p)
        grDevices::dev.off()
        message("  wrote ", png_path)
      }, error = function(e) message("  plot skipped: ", conditionMessage(e)))
    }
  }
  message("gProfiler done. Outputs under ", outdir)
}

if (sys.nframe() == 0) main()
