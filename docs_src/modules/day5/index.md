# Day 5: Metagenomic Profiling

**Date**: September 5, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Metagenomic sequencing and microbiome analysis

## Overview

Day 5 introduces metagenomic analysis for studying complex microbial communities. We'll explore metagenomic sequencing principles, quality control methods specific to metagenomic data, and microbiome profiling techniques using R and QIIME2.

## Learning Objectives

By the end of Day 5, you will be able to:

- Understand principles of metagenomic sequencing and its applications
- Apply quality control methods specific to metagenomic data
- Perform microbiome profiling using computational tools
- Calculate and interpret microbial diversity metrics
- Analyze community composition and structure
- Compare microbiome profiles between samples

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Reproducible workflows with Nextflow and nf-core* | | Mamana Mbiyavanga |
| **10:30** | *Developing a Nextflow pipeline for QC, de novo assembly, quality assessment and annotation* | | Mamana Mbiyavanga |
| **11:30** | **Break** | | |
| **12:00** | *Developing a Nextflow pipeline for QC, de novo assembly, quality assessment and annotation* | | Mamana Mbiyavanga |

## Key Topics

### 1. Metagenomic Sequencing Principles
- Shotgun vs amplicon (16S rRNA) sequencing approaches
- Sample preparation and DNA extraction considerations
- Sequencing depth requirements for different analyses
- Bias sources in metagenomic studies

### 2. Quality Control in Metagenomics
- Host DNA contamination removal
- Adapter and quality trimming for metagenomic reads
- Contamination detection and removal
- Quality assessment metrics specific to metagenomics

### 3. Microbiome Profiling
- Taxonomic classification methods (Kraken2, MetaPhlAn)
- Abundance estimation and normalization
- Reference database selection and limitations
- Functional profiling approaches

### 4. Diversity Analysis
- Alpha diversity: richness, evenness, Shannon, Simpson indices
- Beta diversity: Bray-Curtis, Jaccard, UniFrac distances
- Ordination methods: PCoA, NMDS, correspondence analysis
- Statistical testing for community differences

## Tools and Software

### Metagenomic Analysis Tools
- **Kraken2** - Taxonomic classification
- **MetaPhlAn4** - Marker gene-based profiling  
- **HUMAnN** - Functional profiling
- **QIIME2** - Comprehensive microbiome analysis platform

### R Packages for Microbiome Analysis
- **phyloseq** - Analysis of microbiome data
- **vegan** - Community ecology analysis
- **DESeq2** - Differential abundance testing
- **ggplot2** - Data visualization

## Hands-on Exercises

### Exercise 1: Metagenomic QC (45 minutes)
Process raw metagenomic reads through quality control pipeline.

```bash
# Remove host DNA contamination
kneaddata --input sample_R1.fastq.gz --input sample_R2.fastq.gz \
    --reference-db human_genome --output kneaddata_output/

# Quality assessment of cleaned reads
fastqc kneaddata_output/*_paired_*.fastq
```

### Exercise 2: Taxonomic Profiling (60 minutes)
Classify metagenomic reads and generate abundance profiles.

```bash
# Run Kraken2 taxonomic classification
kraken2 --db kraken2_db --paired \
    sample_R1_paired.fastq sample_R2_paired.fastq \
    --output sample.kraken --report sample_kraken_report.txt

# Generate abundance table
bracken -d kraken2_db -i sample_kraken_report.txt \
    -o sample_bracken.txt -r 150 -l S
```

### Exercise 3: Diversity Analysis in R (75 minutes)
Calculate diversity metrics and visualize community structure.

```r
# Load required packages
library(phyloseq)
library(vegan)
library(ggplot2)

# Import abundance data
otu_table <- read.csv("abundance_table.csv", row.names = 1)
metadata <- read.csv("sample_metadata.csv", row.names = 1)

# Create phyloseq object
ps <- phyloseq(otu_table(otu_table, taxa_are_rows = TRUE),
               sample_data(metadata))

# Calculate alpha diversity
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Perform PCoA ordination
ord <- ordinate(ps, method = "PCoA", distance = "bray")
plot_ordination(ps, ord, color = "sample_type")
```

## Key Concepts

### Metagenomic Study Design
- **Sample collection**: Standardized protocols critical
- **Storage conditions**: -80°C storage, avoid freeze-thaw
- **Controls**: Negative controls, mock communities
- **Replication**: Biological vs technical replicates

### Taxonomic Databases
- **NCBI RefSeq**: Curated reference sequences
- **GTDB**: Genome Taxonomy Database with standardized taxonomy
- **Silva**: rRNA gene database for amplicon studies
- **Database selection**: Impacts classification accuracy

### Diversity Metrics Interpretation
| Metric | Measures | Interpretation |
|--------|----------|----------------|
| Richness | Number of species | Higher = more diverse |
| Shannon | Species diversity | Accounts for evenness |
| Simpson | Dominance | Lower = more diverse |
| Bray-Curtis | Community similarity | 0 = identical, 1 = completely different |

## Clinical Applications

### Diagnostic Metagenomics
- Pathogen identification in clinical samples
- Antimicrobial resistance gene detection
- Outbreak investigation and source tracking
- Culture-independent diagnosis

### Microbiome-Health Associations
- Disease biomarker discovery
- Treatment response prediction
- Personalized medicine approaches
- Host-microbiome interactions

## Assessment Activities

### Individual Analysis
- Complete quality control of metagenomic dataset
- Generate taxonomic abundance profiles
- Calculate diversity metrics for sample comparison
- Create publication-ready visualizations

### Group Discussion
- Compare results from different taxonomic classifiers
- Interpret diversity patterns in clinical context
- Discuss limitations of metagenomic approaches
- Evaluate study design considerations

## Common Challenges

### Low Biomass Samples
```r
# Detection of contamination
contamination_check <- is_contaminant(ps_neg_controls, method="prevalence")

# Remove contaminants
ps_clean <- prune_taxa(!contamination_check$contaminant, ps)
```

### Batch Effects
```r
# Visualize batch effects
plot_ordination(ps, ord, color = "batch", shape = "sample_type")

# Batch correction (if appropriate)
library(limma)
abundance_corrected <- removeBatchEffect(abundance_matrix, batch = metadata$batch)
```

### Statistical Testing
```r
# Test for differential abundance
library(DESeq2)
dds <- phyloseq_to_deseq2(ps, ~ condition)
dds <- DESeq(dds)
results <- results(dds, alpha = 0.05)
```

## Resources

### Essential Reading
- Knight et al. (2018). Best practices for analyzing microbiomes
- Pollock et al. (2018). The madness of microbiome: Attempting to find consensus
- Schloss (2020). Reintroducing mothur for 16S rRNA microbiome analysis

### Software Documentation
- [QIIME2 Tutorials](https://docs.qiime2.org/)
- [phyloseq Tutorial](https://joey711.github.io/phyloseq/)
- [Kraken2 Manual](https://ccb.jhu.edu/software/kraken2/)

### Databases
- [GTDB](https://gtdb.ecogenomic.org/) - Genome Taxonomy Database
- [Silva](https://www.arb-silva.de/) - rRNA gene database
- [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)

## Looking Ahead

**Day 6 Preview**: Co-infection and mutations including:
- Detecting rare pathogens in metagenomic data
- Co-infection analysis methods
- Community shifts and dynamic changes
- Variant calling and phenotype-genotype correlations

## Homework (Optional)

1. Analyze additional metagenomic datasets from different body sites
2. Compare taxonomic profiles from different classification methods
3. Practice diversity analysis with public microbiome datasets
4. Read case studies of clinical metagenomics applications

---

**Key Learning Outcome**: Metagenomic sequencing provides powerful insights into microbial community structure and function, but requires careful attention to experimental design, quality control, and statistical analysis for meaningful biological interpretation.