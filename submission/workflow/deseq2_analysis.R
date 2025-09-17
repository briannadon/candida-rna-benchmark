#!/usr/bin/env Rscript

# DESeq2 Analysis Script for Candida albicans RNAseq
# This script performs differential expression analysis

library(DESeq2)
library(readr)
library(ggplot2)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript deseq2_analysis.R <counts_file> <samples_file> <output_prefix>")
}

counts_file <- args[1]
samples_file <- args[2]
output_prefix <- ifelse(length(args) > 2, args[3], "deseq2")

# Read counts data
cat("Reading counts data from:", counts_file, "\n")
counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", skip=1, comment.char="")

# Remove annotation columns (Chr, Start, End, Strand, Length)
counts <- counts[, -c(1:5), drop=FALSE]

# Extract sample IDs from column names (remove file paths)
colnames(counts) <- gsub(".*\\.(SRR[0-9]+)\\.bam", "\\1", colnames(counts))

# Read sample metadata
cat("Reading sample metadata from:", samples_file, "\n")
samples <- read.csv(samples_file)

# Print debug information
cat("Number of columns in counts:", ncol(counts), "\n")
cat("Counts column names:", paste(colnames(counts), collapse=", "), "\n")
cat("Number of samples in metadata:", nrow(samples), "\n")
cat("Sample IDs:", paste(samples$sample_id, collapse=", "), "\n")

# Ensure sample order matches counts column order
samples <- samples[match(colnames(counts), samples$sample_id), ]
cat("Sample order after matching:", paste(samples$sample_id, collapse=", "), "\n")
cat("Number of samples after matching:", nrow(samples), "\n")

# Add type column (all samples are paired-end)
samples$type <- "paired-end"

# Set row names to sample IDs for colData
rownames(samples) <- samples$sample_id

# Remove the sample_id column as it's now the row names
samples <- samples[, -which(names(samples) == "sample_id"), drop=FALSE]

# Create DESeq2 object
cat("Creating DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~ condition)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Check if we have replicates for proper DESeq2 analysis
if (ncol(counts) <= 2) {
  cat("Warning: Only one replicate per condition detected.\n")
  cat("DESeq2 cannot perform proper differential expression analysis without biological replicates.\n")
  cat("Proceeding with exploratory analysis and normalized counts...\n")
  
  # Estimate size factors and get normalized counts
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized=TRUE)
  
  # Calculate simple fold changes (not statistically tested)
  control_mean <- rowMeans(norm_counts[, samples$condition == "control", drop=FALSE])
  treatment_mean <- rowMeans(norm_counts[, samples$condition == "treatment", drop=FALSE])
  
  # Avoid division by zero
  control_mean[control_mean == 0] <- 0.5
  treatment_mean[treatment_mean == 0] <- 0.5
  
  log2fc <- log2(treatment_mean / control_mean)
  
  # Create results dataframe
  res <- data.frame(
    baseMean = (control_mean + treatment_mean) / 2,
    log2FoldChange = log2fc,
    lfcSE = NA,
    stat = NA,
    pvalue = NA,
    padj = NA
  )
  
  res <- res[order(abs(res$log2FoldChange), decreasing=TRUE), ]
  
  cat("Note: Results are based on simple fold change calculations, not statistical testing.\n")
  cat("Statistical testing requires biological replicates.\n")
  
} else {
  # Run full DESeq2 analysis if we have replicates
  cat("Running DESeq2 analysis...\n")
  dds <- DESeq(dds)
  
  # Get results
  cat("Extracting results...\n")
  res <- results(dds, contrast=c("condition", "treatment", "control"))
  res <- res[order(res$padj), ]
}

# Write results
results_file <- paste0(output_prefix, "_results.csv")
cat("Writing results to:", results_file, "\n")
write.csv(as.data.frame(res), results_file)

# Write normalized counts
norm_counts_file <- paste0(output_prefix, "_normalized_counts.csv")
cat("Writing normalized counts to:", norm_counts_file, "\n")
norm_counts <- counts(dds, normalized=TRUE)
write.csv(norm_counts, norm_counts_file)

# Summary statistics
cat("\n=== Analysis Summary ===\n")
cat("Total genes analyzed:", nrow(res), "\n")
if (ncol(counts) > 2) {
  cat("Genes with padj < 0.05:", sum(res$padj < 0.05, na.rm=TRUE), "\n")
  cat("Genes with |log2FC| > 1 & padj < 0.05:", sum(abs(res$log2FoldChange) > 1 & res$padj < 0.05, na.rm=TRUE), "\n")
  cat("Up-regulated genes (log2FC > 1, padj < 0.05):", sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm=TRUE), "\n")
  cat("Down-regulated genes (log2FC < -1, padj < 0.05):", sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm=TRUE), "\n")
} else {
  cat("Genes with |log2FC| > 2:", sum(abs(res$log2FoldChange) > 2, na.rm=TRUE), "\n")
  cat("Up-regulated genes (log2FC > 2):", sum(res$log2FoldChange > 2, na.rm=TRUE), "\n")
  cat("Down-regulated genes (log2FC < -2):", sum(res$log2FoldChange < -2, na.rm=TRUE), "\n")
  cat("Note: These are exploratory results without statistical significance testing.\n")
}

cat("\nAnalysis completed!\n")
