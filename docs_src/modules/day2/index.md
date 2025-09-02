# Day 2: Introduction to Commandline, High Performance Computing, & Quality Control

**Date**: September 2, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: High Performance Computing, command line proficiency, quality control

## Overview

Day 2 expands on the command line basics from Day 1, introduces high-performance computing resources, and covers essential quality control methods for genomic data. The day includes a special guest lecture on tuberculosis and co-infection.

## Learning Objectives

By the end of Day 2, you will be able to:

- Connect to and navigate high-performance computing systems
- Execute complex command line operations for bioinformatics
- Perform quality control on sequencing data
- Identify species from genomic data
- Understand M. tuberculosis genomics and co-infection patterns

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Introduction High Performance Computing (HPC) – Ilifu* | [Slides](#) • [HPC Guide](#) | Mamana Mbiyavanga |
| **09:30** | *Recap: Introduction to command line interface* | [Tutorial](#) • [Practice](#) | Arash Iranzadeh |
| **10:30** | *Quality checking and control, as well as species identification* | [Slides](#) • [Practical](#) | Arash Iranzadeh |
| **11:30** | **Break** | | |
| **12:00** | *Guest talk: MtB and co-infection* | [Slides](#) | Bethlehem Adnew |

## Key Topics

### 1. High Performance Computing (HPC) Introduction
- Ilifu cluster overview and capabilities
- SSH connections and authentication
- Job submission and queue systems
- Resource allocation and management
- File systems and storage
- Module loading and software access

### 2. Advanced Command Line Operations
- File manipulation and text processing
- Pipes and command chaining
- Shell scripting basics
- Environment variables and PATH
- Process management and monitoring
- Batch job submission

### 3. Quality Control in Genomics
- FastQC analysis and interpretation
- Quality score distributions
- Adapter contamination detection
- Per-base and per-sequence quality metrics
- Identifying low-quality data
- Decision making for data processing

### 4. Species Identification Methods
- Kraken2 taxonomic classification
- Database selection and limitations
- Contamination detection
- Result interpretation and validation
- Confidence assessment
- Troubleshooting classification issues

### 5. M. tuberculosis and Co-infection
- TB genomics and strain typing
- Co-infection patterns and detection
- Clinical implications
- Molecular epidemiology approaches
- Drug resistance mechanisms
- Public health applications

## Tools and Software

### HPC Environment
- **Ilifu cluster** - High-performance computing infrastructure
- **SLURM** - Job scheduling system
- **Module system** - Software environment management
- **SSH clients** - Remote connection tools

### Command Line Tools
- **Advanced bash** - Shell scripting and automation
- **GNU utilities** - Text processing and file manipulation
- **tmux/screen** - Terminal multiplexing
- **rsync** - File synchronization

### Quality Control Tools
- **FastQC** - Sequencing data quality assessment
- **MultiQC** - Aggregate quality reporting
- **Trimmomatic** - Read trimming and filtering
- **Kraken2** - Taxonomic sequence classification

## Hands-on Exercises

### Exercise 1: HPC Connection and Navigation (45 minutes)
Connect to Ilifu cluster and explore the computing environment.

```bash
# Connect to Ilifu via SSH
ssh username@slurm.ilifu.ac.za

# Load required modules
module load fastqc/0.11.9
module load kraken2/2.1.2

# Check available resources
sinfo
squeue
```

### Exercise 2: Command Line Practice (60 minutes)
Practice advanced command line operations for bioinformatics.

```bash
# File processing with pipes
cat sample.fastq | grep -c "^@"  # Count reads
zcat sample.fastq.gz | head -20   # View compressed file

# Batch processing with loops
for file in *.fastq.gz; do
    echo "Processing $file"
    fastqc "$file" -o qc_output/
done

# Text processing and filtering
awk 'NR%4==2{print length($0)}' sample.fastq | sort -n | uniq -c
```

### Exercise 3: Quality Control Analysis (45 minutes)
Analyze sequencing data quality and identify issues.

```bash
# Run FastQC on multiple samples
fastqc *.fastq.gz -o qc_reports/

# Species identification with Kraken2
kraken2 --db minikraken2_v2 --paired sample_R1.fastq sample_R2.fastq \
    --output sample.kraken --report sample_report.txt

# Interpret results
less sample_report.txt
```

## Key Concepts

### HPC Computing Principles
- **Parallel processing**: Utilizing multiple CPU cores
- **Job scheduling**: SLURM queue management
- **Resource sharing**: Fair allocation among users
- **Module system**: Dynamic software environment

### Quality Metrics Understanding
| Metric | Good | Moderate | Poor | Action |
|--------|------|----------|------|---------|
| Per-base quality | >Q30 | Q20-Q30 | <Q20 | Trim/filter |
| Per-sequence quality | >Q25 | Q15-Q25 | <Q15 | Remove reads |
| Adapter content | <5% | 5-10% | >10% | Trim adapters |
| Duplication level | <20% | 20-50% | >50% | Consider PCR bias |

### Species Identification Workflow
1. **Database selection**: Choose appropriate reference
2. **Classification**: Run taxonomic classifier
3. **Result filtering**: Apply confidence thresholds
4. **Validation**: Cross-check with expected species
5. **Contamination check**: Identify unexpected organisms

## Assessment Activities

### Individual Tasks
- Successfully connect to HPC system
- Execute quality control pipeline
- Interpret FastQC reports correctly
- Perform species identification
- Document analysis workflow

### Group Discussion
- Compare quality control results
- Discuss HPC resource management strategies
- Share troubleshooting experiences
- Evaluate species identification confidence

## Common Challenges

### HPC Connection Issues
```bash
# SSH key problems
ssh-keygen -t ed25519 -C "your_email@example.com"
ssh-copy-id username@slurm.ilifu.ac.za

# Module loading issues
module avail    # List available modules
module list     # Show loaded modules
module purge    # Clear all modules
```

### Quality Control Interpretation
- **Low quality at read ends**: Normal for Illumina, trim accordingly
- **Adapter contamination**: Choose appropriate trimming parameters
- **High duplication**: May be biological or technical
- **Unusual GC content**: Could indicate contamination

### Species Identification Problems
```bash
# Database issues
kraken2-build --download-taxonomy --db custom_db
kraken2-build --download-library bacteria --db custom_db

# Low classification rates
# Try different databases or confidence thresholds
kraken2 --confidence 0.1 --db database sample.fastq
```

## Resources

### HPC Documentation
- [Ilifu User Guide](https://docs.ilifu.ac.za/)
- [SLURM Quick Start](https://slurm.schedmd.com/quickstart.html)
- [SSH Key Management](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)

### Quality Control Resources
- [FastQC Manual](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Kraken2 Documentation](https://ccb.jhu.edu/software/kraken2/)
- [Quality Control Best Practices](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4766705/)

### M. tuberculosis Resources
- [TB-Profiler](https://github.com/jodyphelan/TBProfiler)
- [ReSeqTB](https://platform.reseqtb.org/)
- [TBDB](http://genome.tbdb.org/)

## Guest Lecture: MtB and Co-infection

### Speaker: Bethlehem Adnew

#### Key Topics Covered
- **M. tuberculosis genomics**: Strain diversity and typing methods
- **Co-infection dynamics**: TB-HIV and other respiratory pathogens
- **Diagnostic challenges**: Molecular detection in complex samples
- **Treatment implications**: Drug resistance in co-infected patients
- **Epidemiological insights**: Transmission patterns and control strategies

#### Interactive Discussion Points
- Current challenges in TB diagnosis
- Role of genomics in outbreak investigation
- Future directions in TB research
- Integration of genomic and clinical data

## Looking Ahead

**Day 3 Preview**: Genomic Characterization including:
- Genome assembly, quality assessment and annotation
- Multi-locus sequence typing and serotyping
- Antimicrobial resistance gene detection
- Role of mobile genetic elements in AMR spread

---

**Key Learning Outcome**: HPC fundamentals, command line proficiency, and quality control methods provide the essential computational foundation for all subsequent genomic analyses in the course.