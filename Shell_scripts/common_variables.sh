#!/usr/bin/env bash

# Logging helpers
if ! command -v info &> /dev/null; then
    info() { echo "[INFO] $*"; }
fi

if ! command -v error &> /dev/null; then
    error() { echo "[ERROR] $*"; }
fi

###filenames
# Expected to be set by run.sh via prepare_data.py / env
if [ -z "${RPF_filenames:-}" ]; then
    echo "Warning: RPF_filenames not set. Defaulting to empty or check prepare_data.py"
    RPF_filenames=''
fi

if [ -z "${Totals_filenames:-}" ]; then
    echo "Warning: Totals_filenames not set. Defaulting to empty or check prepare_data.py"
    Totals_filenames=''
fi

### threads / memory / adaptors (overridable via env)
if [ -n "${RIBO_SEQ_THREADS:-}" ]; then
    threadN="$RIBO_SEQ_THREADS"
else
    threadN=4
fi

if [ -n "${RIBO_SEQ_BBMAP_MEMORY:-}" ]; then
    bbmap_memory="$RIBO_SEQ_BBMAP_MEMORY"
else
    bbmap_memory='Xmx=4g'
fi

if [ -n "${RIBO_SEQ_RPF_ADAPTOR:-}" ]; then
    RPF_adaptor="$RIBO_SEQ_RPF_ADAPTOR"
else
    RPF_adaptor='TGGAATTCTCGGGTGCCAAGG' # NextFlex small RNA
fi

if [ -n "${RIBO_SEQ_TOTALS_ADAPTOR:-}" ]; then
    Totals_adaptor="$RIBO_SEQ_TOTALS_ADAPTOR"
else
    Totals_adaptor='AGATCGGAAGAG' # LEXOGEN CORALL
fi

### paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPTS_DIR="${PROJECT_ROOT}/Python_scripts"

if [ -n "${RIBO_SEQ_PARENT_DIR:-}" ]; then
    parent_dir="$RIBO_SEQ_PARENT_DIR"
else
    parent_dir="${PROJECT_ROOT}/results"
fi

bcl_dir="${RIBO_SEQ_BCL_DIR:-Path/to/bcl/data}"

fastq_dir=${parent_dir}/fastq_files
fastqc_dir=${parent_dir}/fastQC_files
SAM_dir=${parent_dir}/SAM_files
BAM_dir=${parent_dir}/BAM_files
log_dir=${parent_dir}/logs
counts_dir=${parent_dir}/Counts_files
csv_counts_dir=${parent_dir}/Counts_files/csv_files
csv_R_objects=${parent_dir}/Counts_files/R_objects

STAR_dir=${parent_dir}/STAR
rsem_dir=${parent_dir}/rsem

analysis_dir=${parent_dir}/Analysis

region_counts_dir=${analysis_dir}/region_counts
spliced_counts_dir=${analysis_dir}/spliced_counts
periodicity_dir=${analysis_dir}/periodicity
cds_counts_dir=${analysis_dir}/CDS_counts
UTR5_counts_dir=$analysis_dir/UTR5_counts
codon_counts_dir=${analysis_dir}/codon_counts
most_abundant_transcripts_dir=${analysis_dir}/most_abundant_transcripts
DESeq2_dir=${analysis_dir}/DESeq2_output
reads_summary_dir=${analysis_dir}/reads_summary
fgsea_dir=${analysis_dir}/fgsea

plots_dir=${parent_dir}/plots

summed_counts_plots_dir=${plots_dir}/summed_counts
periodicity_plots_dir=${plots_dir}/periodicity
offset_plots_dir=${plots_dir}/offset
heatmaps_plots_dir=${plots_dir}/heatmaps
DE_analysis_dir=${plots_dir}/DE_analysis
PCA_dir=${plots_dir}/PCAs
Interactive_scatters_dir=${plots_dir}/Interactive_scatters
fgsea_plots_dir=${plots_dir}/fgsea
fgsea_scatters_dir=${plots_dir}/fgsea/scatters
fgsea_interactive_scatters_dir=${plots_dir}/fgsea/Interactive_scatters
read_counts_summary_dir=${plots_dir}/read_counts_summary
binned_plots_dir=${plots_dir}/binned_plots
single_transcript_binned_plots_dir=${plots_dir}/binned_plots/single_transcripts
normalisation_binned_plots_dir=${plots_dir}/binned_plots/normalisation

# Fastas / indices — all overridable; no machine-specific absolute paths
if [ -n "${RIBO_SEQ_FASTA_DIR:-}" ]; then
    fasta_dir="$RIBO_SEQ_FASTA_DIR"
else
    fasta_dir="${PROJECT_ROOT}/reference"
fi

if [ -n "${GENOME_VERSION:-}" ]; then
    genome_version="$GENOME_VERSION"
else
    genome_version='v49'
fi

# Optional subset tag for test/demo refs, e.g. "_chr20". Empty = full genome defaults.
# Example: export RIBO_SEQ_REF_SUFFIX=_chr20
ref_suffix="${RIBO_SEQ_REF_SUFFIX:-}"

if [ -n "${RIBO_SEQ_RRNA_FASTA:-}" ]; then
    rRNA_fasta="$RIBO_SEQ_RRNA_FASTA"
else
    rRNA_fasta=${fasta_dir}/rRNA/sortmerna_rrna.fasta
fi

if [ -n "${RIBO_SEQ_TRNA_FASTA:-}" ]; then
    tRNA_fasta="$RIBO_SEQ_TRNA_FASTA"
else
    tRNA_fasta=${fasta_dir}/tRNA/hg38-mature-tRNAs-dna.fasta
fi

if [ -n "${RIBO_SEQ_PC_FASTA:-}" ]; then
    pc_fasta="$RIBO_SEQ_PC_FASTA"
else
    pc_fasta=${fasta_dir}/GENCODE/${genome_version}/gencode.${genome_version}.pc_transcripts${ref_suffix}_filtered.fa
fi

if [ -n "${RIBO_SEQ_GENOME_FASTA:-}" ]; then
    genome_fasta="$RIBO_SEQ_GENOME_FASTA"
else
    genome_fasta=${fasta_dir}/GENCODE/${genome_version}/gencode.${genome_version}.dna${ref_suffix}.fa.gz
fi

if [ -n "${RIBO_SEQ_RSEM_INDEX:-}" ]; then
    rsem_index="$RIBO_SEQ_RSEM_INDEX"
else
    rsem_index=${fasta_dir}/GENCODE/${genome_version}/rsem_bowtie2_index/gencode.${genome_version}.pc_transcripts${ref_suffix}_filtered
fi

if [ -n "${RIBO_SEQ_STAR_INDEX:-}" ]; then
    STAR_index="$RIBO_SEQ_STAR_INDEX"
else
    STAR_index=${fasta_dir}/GENCODE/${genome_version}/STAR_index
fi

if [ -n "${RIBO_SEQ_STAR_GTF:-}" ]; then
    STAR_GTF="$RIBO_SEQ_STAR_GTF"
else
    STAR_GTF=${fasta_dir}/GENCODE/${genome_version}/gencode.${genome_version}.annotation${ref_suffix}.gtf
fi

most_abundant_fasta=$most_abundant_transcripts_dir/most_abundant_transcripts.fa

if [ -n "${RIBO_SEQ_REGION_LENGTHS:-}" ]; then
    region_lengths="$RIBO_SEQ_REGION_LENGTHS"
else
    region_lengths=${fasta_dir}/GENCODE/${genome_version}/transcript_info/gencode.${genome_version}.pc_transcripts${ref_suffix}_region_lengths.csv
fi
