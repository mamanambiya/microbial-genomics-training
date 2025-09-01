# Day 2: Data Analysis Fundamentals

**Duration**: Full day (09:00-13:00)  
**Focus**: Quality control, read processing, and assembly introduction

## Overview

Day 2 introduces the fundamental concepts of genomic data analysis. We'll work with real sequencing data to understand quality metrics, perform essential preprocessing steps, and begin exploring genome assembly approaches.

## Learning Objectives

By the end of Day 2, you will be able to:

- Assess the quality of raw sequencing data
- Perform quality control and read processing
- Understand different sequencing technologies and their characteristics
- Apply basic filtering and trimming operations
- Recognize common data quality issues and their solutions

## Schedule

| Time | Topic | Instructor |
|------|-------|------------|
| 09:00-09:30 | Data Quality Overview | Sindiswa Lukhele |
| 09:30-10:30 | [Quality Control Analysis](quality-control.md) | Sindiswa Lukhele |
| 10:30-10:45 | *Coffee Break* | |
| 10:45-11:30 | [Read Processing & Filtering](read-processing.md) | Arash Iranzadeh |
| 11:30-12:15 | [Assembly Basics](assembly.md) | Ephifania Geza |
| 12:15-13:00 | Hands-on Practice | All Instructors |

## Key Topics

### 1. Understanding Sequencing Data
- FASTQ format and quality scores
- Illumina, Oxford Nanopore, and PacBio technologies
- Common artifacts and error patterns
- Quality metrics interpretation

### 2. Quality Control Analysis
- FastQC reports and interpretation
- MultiQC for batch analysis
- Identifying problematic samples
- Quality threshold decisions

### 3. Read Processing
- Adapter trimming strategies
- Quality filtering approaches
- Length filtering considerations
- Paired-end read handling

### 4. Assembly Introduction
- Assembly concepts and challenges
- Short-read vs long-read approaches
- Assembly quality metrics
- Choosing appropriate assemblers

## Datasets Used

### Primary Dataset: *E. coli* Outbreak Strains
- **Source**: Simulated outbreak investigation
- **Technology**: Illumina paired-end (2×150bp)
- **Coverage**: 50-100x
- **Samples**: 12 outbreak isolates + 3 controls

### Secondary Dataset: *M. tuberculosis* Clinical Isolates
- **Source**: Clinical surveillance study
- **Technology**: Illumina paired-end (2×150bp)
- **Coverage**: 80-120x
- **Samples**: 6 diverse strains

## Tools Introduced

### Quality Control
- **FastQC** - Individual sample QC
- **MultiQC** - Aggregate QC reporting
- **Qualimap** - Mapping quality assessment

### Read Processing  
- **Trimmomatic** - Adapter trimming and filtering
- **Cutadapt** - Flexible adapter removal
- **Fastp** - All-in-one preprocessing

### Assembly Tools
- **SPAdes** - De novo genome assembly
- **Unicycler** - Hybrid assembly approach
- **QUAST** - Assembly quality assessment

## Hands-on Exercises

### Exercise 1: Quality Assessment (45 minutes)
Analyze raw sequencing data quality and identify potential issues.

```bash
# Basic FastQC analysis
fastqc *.fastq.gz -o qc_reports/

# Generate MultiQC report
multiqc qc_reports/ -o multiqc_output/
```

### Exercise 2: Read Processing (60 minutes)
Apply appropriate trimming and filtering to improve data quality.

```bash
# Trimmomatic processing
trimmomatic PE sample_R1.fastq.gz sample_R2.fastq.gz \
    output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
    output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### Exercise 3: Assembly Comparison (45 minutes)
Compare assembly results from different approaches and parameters.

```bash
# SPAdes assembly
spades.py -1 processed_R1.fastq.gz -2 processed_R2.fastq.gz \
    -o spades_output/ --careful

# Assembly quality assessment
quast.py spades_output/scaffolds.fasta -o quast_results/
```

## Key Concepts

### Quality Scores (Phred Scale)
- Q20 = 1% error rate (99% accuracy)
- Q30 = 0.1% error rate (99.9% accuracy)  
- Q40 = 0.01% error rate (99.99% accuracy)

### Assembly Metrics
- **N50**: Scaffold length where 50% of assembly is in longer scaffolds
- **Coverage**: Average depth of sequencing reads
- **Contiguity**: Number and size distribution of contigs

### Common Issues
- **Adapter contamination**: Affects assembly quality
- **Low quality regions**: Lead to assembly errors
- **Coverage bias**: Uneven representation of genome regions

## Assessment Activities

### Practical Assessment
- Quality control report interpretation
- Successful read processing workflow
- Assembly quality evaluation
- Troubleshooting common problems

### Discussion Topics
- When to reject samples based on quality
- Optimal trimming parameters for different applications
- Assembly strategy selection criteria

## Common Problems & Solutions

### Low Quality Data
```bash
# More aggressive trimming
trimmomatic PE input_R1.fastq.gz input_R2.fastq.gz \
    output_R1.fastq.gz output_R1_unpaired.fastq.gz \
    output_R2.fastq.gz output_R2_unpaired.fastq.gz \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:50
```

### Adapter Contamination
```bash
# Custom adapter removal
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
    input_R1.fastq.gz input_R2.fastq.gz
```

### Assembly Fragmentation
```bash
# Try different k-mer sizes
spades.py -1 R1.fastq.gz -2 R2.fastq.gz \
    -o assembly_k55,77,99/ -k 55,77,99
```

## Resources

### Documentation
- [FastQC Manual](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [SPAdes Manual](https://github.com/ablab/spades)
- [Trimmomatic Documentation](http://www.usadellab.org/cms/?page=trimmomatic)

### Reference Papers
- Bolger et al. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data
- Bankevich et al. (2012). SPAdes: a new genome assembly algorithm
- Ewels et al. (2016). MultiQC: summarize analysis results for multiple tools

### Online Resources
- [Galaxy Training Materials](https://training.galaxyproject.org/)
- [Bioinformatics Workbook](https://bioinformaticsworkbook.org/)

## Looking Ahead

**Day 3 Preview**: Pathogen-specific genomic analysis including:
- *M. tuberculosis* strain characterization
- *V. cholerae* outbreak investigation
- Antimicrobial resistance detection

## Homework (Optional)

1. Process additional datasets provided on course shared folder
2. Explore different assembly parameters and compare results
3. Read assigned papers on assembly methods

---

**Key Takeaway**: Quality control is the foundation of all downstream analyses. Invest time in understanding your data quality before proceeding to interpretation!