#!/bin/bash

# Candida albicans RNAseq Analysis Pipeline
# This script performs:
# 1. Quality control with FastQC
# 2. Read alignment with STAR
# 3. Read counting with featureCounts
# 4. Differential expression analysis with DESeq2

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SAMPLES_FILE="data/samples.csv"
REFERENCE="data/candida_reference.fasta.gz"
ANNOTATION="data/candida_annotation.gtf"
OUTDIR="results"
STAR_INDEX="star_index"
THREADS=8
MEMORY="32G"

# Source conda/mamba initialization
if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source "/opt/conda/etc/profile.d/conda.sh"
fi
if [ -f "/opt/conda/etc/profile.d/mamba.sh" ]; then
    source "/opt/conda/etc/profile.d/mamba.sh"
fi

# Create output directories
mkdir -p "${OUTDIR}/fastqc"
mkdir -p "${OUTDIR}/star_index"
mkdir -p "${OUTDIR}/star_align"
mkdir -p "${OUTDIR}/featurecounts"
mkdir -p "${OUTDIR}/deseq2"
mkdir -p "${OUTDIR}/reports"

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to run commands with conda environment activation
run_in_env() {
    local env_name="$1"
    shift
    log "Activating $env_name environment..."
    conda activate "$env_name"
    "$@"
    local exit_code=$?
    conda deactivate
    return $exit_code
}

# Function to check if a conda environment exists
check_env() {
    if ! conda env list | grep -q "^$1\s"; then
        log "ERROR: Conda environment '$1' not found."
        exit 1
    fi
}

# Check required conda environments
log "Checking required conda environments..."
check_env "fastqc"
check_env "star"
check_env "subread"
check_env "samtools"
check_env "R"

# Parse samples file
log "Reading sample metadata from ${SAMPLES_FILE}"
declare -A SAMPLES
declare -A CONDITIONS

while IFS=, read -r sample_id condition reads; do
    # Skip header
    if [[ "$sample_id" == "sample_id" ]]; then
        continue
    fi
    
    # Remove quotes from reads paths
    reads=$(echo "$reads" | tr -d '"')
    
    SAMPLES["$sample_id"]="$reads"
    CONDITIONS["$sample_id"]="$condition"
    
done < "$SAMPLES_FILE"
log "Samples: ${SAMPLES[*]}, Condition: ${CONDITIONS[*]}"

# Step 1: Quality Control with FastQC
log "Starting quality control with FastQC..."
for sample_id in "${!SAMPLES[@]}"; do
    reads="${SAMPLES[$sample_id]}"
    log "Processing $sample_id: $reads"
    
    # Parse comma-separated read files
    IFS=',' read -ra READ_FILES <<< "$reads"
    for read_file in "${READ_FILES[@]}"; do
        run_in_env "fastqc" fastqc "$read_file" -o "${OUTDIR}/fastqc" -t $THREADS
    done
done

# Step 2: Build STAR index if it doesn't exist
if [[ ! -d "${STAR_INDEX}" ]]; then
    log "Building STAR index..."
    mkdir -p "${STAR_INDEX}"
    
    # Decompress reference and annotation if needed
    if [[ "$REFERENCE" == *.gz ]]; then
        REFERENCE_DECOMPRESSED="${REFERENCE%.gz}"
        if [[ ! -f "$REFERENCE_DECOMPRESSED" ]]; then
            log "Decompressing reference genome..."
            gunzip -c "$REFERENCE" > "$REFERENCE_DECOMPRESSED"
        fi
        REFERENCE="$REFERENCE_DECOMPRESSED"
    fi
    
    if [[ "$ANNOTATION" == *.gz ]]; then
        ANNOTATION_DECOMPRESSED="${ANNOTATION%.gz}"
        if [[ ! -f "$ANNOTATION_DECOMPRESSED" ]]; then
            log "Decompressing annotation file..."
            gunzip -c "$ANNOTATION" > "$ANNOTATION_DECOMPRESSED"
        fi
        ANNOTATION="$ANNOTATION_DECOMPRESSED"
    fi
    
    run_in_env "star" STAR --runMode genomeGenerate \
         --genomeDir "${STAR_INDEX}" \
         --genomeFastaFiles "$REFERENCE" \
         --sjdbGTFfile "$ANNOTATION" \
         --runThreadN $THREADS \
         --genomeSAindexNbases 12
else
    log "STAR index already exists, skipping index building..."
fi

# Step 3: STAR alignment
log "Starting STAR alignment..."
declare -A BAM_FILES

for sample_id in "${!SAMPLES[@]}"; do
    reads="${SAMPLES[$sample_id]}"
    log "Aligning $sample_id..."
    
    # Parse paired-end reads
    IFS=',' read -ra READ_FILES <<< "$reads"
    read1="${READ_FILES[0]}"
    read2="${READ_FILES[1]}"
    
    run_in_env "star" STAR --genomeDir "${STAR_INDEX}" \
         --readFilesIn "$read1" "$read2" \
         --runThreadN $THREADS \
         --outFileNamePrefix "${OUTDIR}/star_align/${sample_id}." \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard
    
    # Move and index BAM file
    mv "${OUTDIR}/star_align/${sample_id}.Aligned.sortedByCoord.out.bam" "${OUTDIR}/star_align/${sample_id}.bam"
    BAM_FILES["$sample_id"]="${OUTDIR}/star_align/${sample_id}.bam"
    
    log "Indexing BAM file for $sample_id..."
    run_in_env "samtools" samtools index "${OUTDIR}/star_align/${sample_id}.bam"
done

# Step 4: Read counting with featureCounts
log "Starting read counting with featureCounts..."
BAM_LIST=""
for bam_file in "${BAM_FILES[@]}"; do
    BAM_LIST="$BAM_LIST $bam_file"
done

run_in_env "subread" featureCounts -p -a "$ANNOTATION" \
              -o "${OUTDIR}/featurecounts/counts.txt" \
              -T $THREADS \
              $BAM_LIST

# Step 5: DESeq2 analysis
log "Starting DESeq2 differential expression analysis..."
run_in_env "R" Rscript workflow/deseq2_analysis.R \
    "${OUTDIR}/featurecounts/counts.txt" \
    "$SAMPLES_FILE" \
    "${OUTDIR}/deseq2/deseq2"

# Generate MultiQC report
log "Generating MultiQC report..."
if conda env list | grep -q "^multiqc\s"; then
    run_in_env "multiqc" multiqc "${OUTDIR}" -o "${OUTDIR}/reports"
else
    log "MultiQC environment not found, skipping report generation..."
fi

# Final summary
log "Pipeline completed successfully!"
log "Results available in: ${OUTDIR}"
log " - FastQC reports: ${OUTDIR}/fastqc"
log " - STAR alignments: ${OUTDIR}/star_align"
log " - Read counts: ${OUTDIR}/featurecounts/counts.txt"
log " - DESeq2 results: ${OUTDIR}/deseq2/"
log " - MultiQC report: ${OUTDIR}/reports/multiqc_report.html"

# Clean up temporary files
if [[ -f "${REFERENCE_DECOMPRESSED:-}" ]]; then
    rm "$REFERENCE_DECOMPRESSED"
fi
if [[ -f "${ANNOTATION_DECOMPRESSED:-}" ]]; then
    rm "$ANNOTATION_DECOMPRESSED"
fi

exit 0
