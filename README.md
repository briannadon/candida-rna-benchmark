

# Candida albicans RNAseq Analysis

## 1. Data sources

This task is based on publicly available sequencing data from a study of **Candida albicans gene expression**. The dataset includes RNAseq samples and was originally sequenced using **Illumina sequencing platform**.

The FASTQ files are stored in `data/` and are used as the inputs for the workflow. They are subsampled down to 1M reads via `seqtk sample -s100 sample.fq.gz 1000000`.
- `data/SRR24303833_1.fastq.subsampled.gz`: Forward reads
- `data/SRR24303833_2.fastq.subsampled.gz`: Reverse reads
- `data/candida_reference.fasta.gz`: Reference genome
- `data/candida_annotation.gff3.gz`: Genome annotation

---

## 2. How to download

The data was downloaded from NCBI SRA using the SRA Toolkit. The SRA accession numbers are listed in `SRAs.txt`.

### Example using SRA Toolkit

```bash
# Install SRA Toolkit
conda install -c bioconda sra-tools

# Download specific SRA accession
prefetch SRR24303833
fasterq-dump SRR24303833 --split-files
```

---

## 3. Pre-processing / subsampling

The data has been pre-processed and is ready for analysis. No additional subsampling was performed as this is a single sample workflow skeleton.

---

## 4. How the workflow works

The workflow is built using Nextflow and performs RNAseq analysis through the following steps:

### Step 1 – Quality Control

**Purpose:** Assess read quality and generate QC reports
**Tools:** `FastQC`
**Inputs:** Raw FASTQ files from `data/`
**Outputs:** HTML and ZIP QC reports in `results/fastqc/`
**Command:**
```bash
fastqc data/SRR24303833_1.fastq data/SRR24303833_2.fastq
```

---

### Step 2 – Genome Indexing

**Purpose:** Build STAR index for alignment
**Tools:** `STAR`
**Inputs:** Reference genome (`data/candida_reference.fasta.gz`) and annotation (`data/candida_annotation.gff3.gz`)
**Outputs:** STAR genome index in `results/star_index/`
**Command:**
```bash
STAR --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles data/candida_reference.fasta.gz \
     --sjdbGTFfile data/candida_annotation.gff3.gz
```

---

### Step 3 – Read Alignment

**Purpose:** Map reads to reference genome
**Tools:** `STAR`
**Inputs:** FASTQ files and STAR index
**Outputs:** Aligned BAM files in `results/star_align/`
**Command:**
```bash
STAR --genomeDir star_index \
     --readFilesIn data/SRR24303833_1.fastq data/SRR24303833_2.fastq \
     --outSAMtype BAM SortedByCoordinate
```

---

### Step 4 – Read Counting

**Purpose:** Count reads per gene/feature
**Tools:** `featureCounts`
**Inputs:** Aligned BAM files and annotation
**Outputs:** Count matrix in `results/featurecounts/counts.txt`
**Command:**
```bash
featureCounts -a data/candida_annotation.gff3.gz \
              -o counts.txt \
              results/star_align/SRR24303833.Aligned.sortedByCoord.out.bam
```

---

### Step 5 – Differential Expression Analysis

**Purpose:** Identify differentially expressed genes (requires multiple samples)
**Tools:** `DESeq2` (R package)
**Inputs:** Count matrix and sample metadata
**Outputs:** DESeq2 results and plots in `results/deseq2/`
**Command:**
```R
# R script for DESeq2 analysis
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)
dds <- DESeq(dds)
results <- results(dds)
```

---

## 5. Running the workflow

The workflow is executed using Nextflow with Conda for dependency management:

```bash
# Test the workflow structure
./run_test.sh

# Run the complete pipeline
nextflow run main.nf -profile standard

# Run quick test (skip alignment and DESeq2)
nextflow run main.nf -profile test
```

Software dependencies are managed through Conda using the `environment.yml` file. Nextflow will automatically create and use the Conda environment.

