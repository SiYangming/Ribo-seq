#Set the parent directory (same meaning as Shell_scripts/common_variables.sh)
parent_dir <- Sys.getenv("RIBO_SEQ_PARENT_DIR")
if (parent_dir == "") {
  parent_dir <- file.path(dirname(getwd()), "results")
  warning("RIBO_SEQ_PARENT_DIR not set. Using default: ", parent_dir)
}

# Project root (sibling of results/ when using default layout)
project_root <- Sys.getenv("RIBO_SEQ_PROJECT_ROOT")
if (project_root == "") {
  if (basename(parent_dir) == "results") {
    project_root <- dirname(parent_dir)
  } else {
    project_root <- parent_dir
  }
}

# Genome Version
genome_version <- Sys.getenv("GENOME_VERSION")
if (genome_version == "") {
  genome_version <- "v49"
}

# Optional subset suffix for test refs, e.g. "_chr20"
ref_suffix <- Sys.getenv("RIBO_SEQ_REF_SUFFIX")

# Fasta Directory
fasta_dir <- Sys.getenv("RIBO_SEQ_FASTA_DIR")
if (fasta_dir == "") {
  fasta_dir <- file.path(project_root, "reference")
}

resolve_region_lengths <- function() {
  env_path <- Sys.getenv("RIBO_SEQ_REGION_LENGTHS")
  if (nzchar(env_path)) {
    return(env_path)
  }
  info_dir <- file.path(fasta_dir, "GENCODE", genome_version, "transcript_info")
  candidates <- c(
    file.path(info_dir, paste0("gencode.", genome_version, ".pc_transcripts", ref_suffix, "_region_lengths.csv")),
    file.path(info_dir, paste0("gencode.", genome_version, ".pc_transcripts_region_lengths.csv")),
    file.path(info_dir, paste0("gencode.", genome_version, ".pc_transcripts_chr20_region_lengths.csv"))
  )
  for (p in candidates) {
    if (file.exists(p)) {
      return(p)
    }
  }
  candidates[[1]]
}

region_lengths_file <- resolve_region_lengths()

#set sample names
RPF_filenames_env <- Sys.getenv("RIBO_SEQ_RPF_FILENAMES")
if (RPF_filenames_env != "") {
  RPF_sample_names <- strsplit(RPF_filenames_env, " ")[[1]]
} else {
  RPF_sample_names <- character(0)
}

Totals_filenames_env <- Sys.getenv("RIBO_SEQ_TOTALS_FILENAMES")
if (Totals_filenames_env != "") {
  Total_sample_names <- strsplit(Totals_filenames_env, " ")[[1]]
} else {
  Total_sample_names <- character(0)
}

# Prefer info.csv next to project root (not hardcoded sample names)
info_candidates <- c(
  Sys.getenv("RIBO_SEQ_INFO_CSV"),
  file.path(project_root, "info.csv"),
  file.path(dirname(parent_dir), "info.csv")
)
info_file <- ""
for (cand in info_candidates) {
  if (nzchar(cand) && file.exists(cand)) {
    info_file <- cand
    break
  }
}

use_default_info <- TRUE

if (nzchar(info_file)) {
  info_data <- read.csv(info_file, stringsAsFactors = FALSE)

  rpf_rows <- info_data[info_data$type == "riboseq", , drop = FALSE]
  if (length(RPF_sample_names) == 0 && nrow(rpf_rows) > 0) {
    RPF_sample_names <- rpf_rows$sample
  }
  if (nrow(rpf_rows) == length(RPF_sample_names) && length(RPF_sample_names) > 0) {
    rpf_rows$replicate <- ave(seq_along(rpf_rows$treatment), rpf_rows$treatment, FUN = seq_along)
    RPF_sample_info <- data.frame(
      sample = RPF_sample_names,
      condition = rpf_rows$treatment,
      replicate = factor(rpf_rows$replicate)
    )
    use_default_info <- FALSE
  }

  total_rows <- info_data[info_data$type == "rnaseq", , drop = FALSE]
  if (length(Total_sample_names) == 0 && nrow(total_rows) > 0) {
    Total_sample_names <- total_rows$sample
  }
  if (nrow(total_rows) == length(Total_sample_names) && length(Total_sample_names) > 0) {
    total_rows$replicate <- ave(seq_along(total_rows$treatment), total_rows$treatment, FUN = seq_along)
    Total_sample_info <- data.frame(
      sample = Total_sample_names,
      condition = total_rows$treatment,
      replicate = factor(total_rows$replicate)
    )
  }
}

if (use_default_info) {
  if (length(RPF_sample_names) == 0) {
    RPF_sample_names <- c("Ctrl_RPFs_1", "Ctrl_RPFs_2", "Ctrl_RPFs_3", "Treatment_RPFs_1", "Treatment_RPFs_2", "Treatment_RPFs_3")
  }
  n <- length(RPF_sample_names)
  half <- max(1L, as.integer(n / 2))
  RPF_sample_info <- data.frame(
    sample = RPF_sample_names,
    condition = c(rep("Ctrl", half), rep("Treatment", n - half)),
    replicate = factor(rep_len(seq_len(half), n))
  )
  if (length(Total_sample_names) == 0) {
    Total_sample_names <- c("Ctrl_Totals_1", "Ctrl_Totals_2", "Ctrl_Totals_3", "Treatment_Totals_1", "Treatment_Totals_2", "Treatment_Totals_3")
  }
  n2 <- length(Total_sample_names)
  half2 <- max(1L, as.integer(n2 / 2))
  Total_sample_info <- data.frame(
    sample = Total_sample_names,
    condition = c(rep("Ctrl", half2), rep("Treatment", n2 - half2)),
    replicate = factor(rep_len(seq_len(half2), n2))
  )
}
