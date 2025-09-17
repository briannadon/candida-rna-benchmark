

# Candida albicans RNAseq Analysis

## 1. Data sources

This task is based on publicly available sequencing data from a study of **Candida albicans gene expression**. The dataset includes RNAseq samples and was originally sequenced using **Illumina sequencing platform**.

The FASTQ files are stored in `submission/data/` and are used as the inputs for the workflow. They are subsampled down to 1M reads via `seqtk sample -s100 sample.fq.gz 1000000`.
- `submission/data/SRR24303833_1.fastq.subsampled`: Forward reads (control sample)
- `submission/data/SRR24303833_2.fastq.subsampled`: Reverse reads (control sample)
- `submission/data/SRR28012743_1.fastq.subsampled`: Forward reads (treatment sample)
- `submission/data/SRR28012743_2.fastq.subsampled`: Reverse reads (treatment sample)
- `submission/data/candida_reference.fasta`: Reference genome
- `submission/data/candida_annotation.gtf`: Genome annotation (GTF format)
- `submission/data/candida_annotation.gff3`: Genome annotation (GFF3 format)
- `submission/data/samples.csv`: Sample metadata with conditions

---

## 2. How to download

The data was downloaded from NCBI SRA using the SRA Toolkit. The SRA accession numbers are listed in `submission/data/SRAs.txt`.

### Example using SRA Toolkit

```bash
# Install SRA Toolkit
conda install -c bioconda sra-tools

# Download specific SRA accessions
prefetch SRR24303833
prefetch SRR28012743
fasterq-dump SRR24303833 --split-files
fasterq-dump SRR28012743 --split-files
```

---

## 3. Pre-processing / subsampling

The data has been pre-processed and is ready for analysis. The workflow includes:

1. **FASTA header conversion**: The `convert_fasta_headers.py` script converts ENA format headers to GFF3-compatible format for proper alignment
2. **Subsampling**: FASTQ files are subsampled to 1M reads for faster processing
3. **Sample metadata**: The `samples.csv` file defines sample conditions (control vs treatment) for differential expression analysis

---

## 4. How the workflow works

The workflow is implemented as a bash script (`submission/workflow/rnaseq_pipeline.sh`) and performs RNAseq analysis through the following steps:

### Step 1 – Quality Control

**Purpose:** Assess read quality and generate QC reports
**Tools:** `FastQC`
**Inputs:** Raw FASTQ files from `submission/data/`
**Outputs:** HTML and ZIP QC reports in `submission/workflow/results/fastqc/`
**Command:**
```bash
fastqc submission/data/SRR24303833_1.fastq.subsampled submission/data/SRR24303833_2.fastq.subsampled
fastqc submission/data/SRR28012743_1.fastq.subsampled submission/data/SRR28012743_2.fastq.subsampled
```

---

### Step 2 – Genome Indexing

**Purpose:** Build STAR index for alignment
**Tools:** `STAR`
**Inputs:** Reference genome (`submission/data/candida_reference.fasta`) and annotation (`submission/data/candida_annotation.gtf`)
**Outputs:** STAR genome index in `submission/workflow/results/star_index/`
**Command:**
```bash
STAR --runMode genomeGenerate \
     --genomeDir submission/workflow/results/star_index \
     --genomeFastaFiles submission/data/candida_reference.fasta \
     --sjdbGTFfile submission/data/candida_annotation.gtf \
     --runThreadN 8 \
     --genomeSAindexNbases 12
```

---

### Step 3 – Read Alignment

**Purpose:** Map reads to reference genome
**Tools:** `STAR`
**Inputs:** FASTQ files and STAR index
**Outputs:** Aligned BAM files in `submission/workflow/results/star_align/`
**Command:**
```bash
STAR --genomeDir submission/workflow/results/star_index \
     --readFilesIn submission/data/SRR24303833_1.fastq.subsampled submission/data/SRR24303833_2.fastq.subsampled \
     --runThreadN 8 \
     --outFileNamePrefix submission/workflow/results/star_align/SRR24303833. \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard
```

---

### Step 4 – Read Counting

**Purpose:** Count reads per gene/feature
**Tools:** `featureCounts`
**Inputs:** Aligned BAM files and annotation
**Outputs:** Count matrix in `submission/workflow/results/featurecounts/counts.txt`
**Command:**
```bash
featureCounts -p -a submission/data/candida_annotation.gtf \
              -o submission/workflow/results/featurecounts/counts.txt \
              -T 8 \
              submission/workflow/results/star_align/SRR24303833.bam \
              submission/workflow/results/star_align/SRR28012743.bam
```

---

### Step 5 – Differential Expression Analysis

**Purpose:** Identify differentially expressed genes between control and treatment conditions
**Tools:** `DESeq2` (R package)
**Inputs:** Count matrix and sample metadata (`submission/data/samples.csv`)
**Outputs:** DESeq2 results and normalized counts in `submission/workflow/results/deseq2/`
**Command:**
```bash
Rscript submission/workflow/deseq2_analysis.R \
    submission/workflow/results/featurecounts/counts.txt \
    submission/data/samples.csv \
    submission/workflow/results/deseq2/deseq2
```

**Note:** The DESeq2 analysis script (`submission/workflow/deseq2_analysis.R`) handles both scenarios:
- **With replicates**: Performs full statistical analysis with p-values and adjusted p-values
- **Without replicates**: Provides exploratory analysis with fold change calculations (current dataset has only 2 samples)

---

## 5. Running the workflow

The workflow can be executed in two ways:

### Option 1: Direct execution (requires conda environments)

```bash
# Navigate to the submission directory
cd submission

# Run the complete pipeline
bash workflow/rnaseq_pipeline.sh
```

### Option 2: Docker execution (recommended)

```bash
# Build the Docker container (includes all dependencies)
docker build -t bioinf-bench .

# Run the complete pipeline in Docker
bash submission/workflow/run_full_pipeline_in_docker.sh
```

### Software Dependencies

The workflow requires the following conda environments:
- `fastqc`: For quality control
- `star`: For genome indexing and read alignment
- `subread`: For read counting with featureCounts
- `samtools`: For BAM file indexing
- `R`: For DESeq2 differential expression analysis
- `multiqc`: For generating summary reports (optional)

All dependencies are pre-installed in the Docker container (`bioinf-bench`).

