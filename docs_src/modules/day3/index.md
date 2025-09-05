# Day 3: Accelerating Bioinformatics: HPC, QC, and Species Identification Essentials

**Date**: September 3, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: HPC infrastructure, quality control, species identification

## Overview

Day 3 continues building computational skills with an introduction to the Ilifu HPC infrastructure, then introduces essential genomic characterization techniques including quality control, species identification. These foundational skills are critical for all downstream genomic analyses.

## Learning Objectives

By the end of Day 3, you will be able to:

- Connect to and navigate the Ilifu high-performance computing cluster
- Submit and manage jobs using the SLURM scheduler
- Understand resource allocation and job queue management on HPC systems
- Perform quality checking and control on sequencing data using FastQC
- Identify species from genomic data using Kraken2 and other tools

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Introduction High Performance Computing (HPC) – Ilifu* | [Notes](../../day2/hpc-ilifu-training.md) • [Practical 1](../../day2/slurm-practical-tutorial.md) • [Practical 2](../../day2/unix-commands-pathogen-examples.md) | Mamana Mbiyavanga |
| **11:00** | **Break** | | |
| **11:30** | *Quality checking and control, as well as species identification* | [Practical](https://github.com/Arash-Iranzadeh/Microbial-Genomics/blob/main/scripts/qc_cleaning_detection.sh) | Arash Iranzadeh |



## Key Topics

### 1. QC and Species Identification Essentials
- Quality checking, and adapter and low read quality filtering
- Contamination detection and removal
- Species identification or confirmation

### 2. High Performance Computing (HPC) - Ilifu Infrastructure
- Introduction to cluster computing concepts
- Ilifu cluster architecture and capabilities
- SSH connections and authentication
- SLURM job scheduling system
- Resource allocation (CPU, memory, time)
- Module system for software management

### 3. SLURM Job Management
- Writing and submitting batch scripts
- Interactive vs batch jobs
- Job monitoring and queue management
- Resource specification and optimization
- Output and error file handling
- Best practices for efficient HPC usage

## Tools and Software

### HPC Environment
- **Ilifu cluster** - High-performance computing infrastructure
- **SLURM** - Job scheduling system
- **Module system** - Software environment management
- **SSH clients** - Remote connection tools

### Quality Control Tools
- **FASTQC** - Quality checking
- **MULTIQC** - Quality checking and amalgation of reports
- **Trimmomatic** - Filter adapters and low quality reads
- **Fastp** -  Filter adapters and low quality reads
- **KRAKEN2** - Species Identification
  
## Hands-on Exercises

### Exercise 2: Ilifu HPC Connection and Setup (30 minutes)
Connect to the Ilifu cluster and set up your working environment.

```bash
# Connect to Ilifu via SSH
ssh username@slurm.ilifu.ac.za

# Explore the HPC environment
pwd                    # Check current directory
module avail          # List available software
module load python    # Load software module

# Check cluster resources
sinfo                 # View cluster partitions
squeue               # Check job queue
```

### Exercise 3: SLURM Job Submission (60 minutes)
Learn to submit and manage jobs on the HPC cluster.

```bash
# Create a simple batch script
cat > my_first_job.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=test_job
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

echo "Hello from HPC!"
echo "Running on node: $HOSTNAME"
date
EOF

# Submit the job
sbatch my_first_job.sh

# Monitor job progress
squeue -u $USER
```

### Exercise 1: Species Identification (60 minutes)
Species identification and contamination screening

# Contamination screening
kraken2 --db minikraken2_v2 assembly_output/scaffolds.fasta --report contamination_check.txt
```

## Key Concepts

### HPC Computing Principles
- **Cluster architecture**: Login nodes vs compute nodes
- **Job scheduling**: SLURM queue management and priority
- **Resource allocation**: CPU, memory, and time specifications
- **Module system**: Managing software environments
- **Parallel processing**: Utilizing multiple cores efficiently

### SLURM Job Management
| Component | Description | Example |
|-----------|-------------|---------|
| Partition | Compute resource group | `main`, `gpu`, `bigmem` |
| Job State | Current job status | `PD` (pending), `R` (running) |
| Resources | CPU/Memory/Time | `--cpus-per-task=4 --mem=8G` |
| Output | Job results and logs | `slurm-jobid.out` |

## Assessment Activities

### Individual Analysis
- Write and submit a SLURM batch script
- Monitor job status and retrieve results
- Complete genome assembly workflow
- Perform quality assessment and interpretation

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

### SLURM Job Troubleshooting
```bash
# Job stuck in pending
scontrol show job <jobid>  # Check job details
squeue -j <jobid>          # Check specific job

# Resource issues
sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed
```

## Resources

### HPC Documentation
- [Ilifu User Guide](https://docs.ilifu.ac.za/)
- [SLURM Quick Start](https://slurm.schedmd.com/quickstart.html)
- [SSH Key Management](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)

### Assembly Issues
```bash
# Low coverage assemblies
spades.py --careful --cov-cutoff 5 -1 R1.fastq -2 R2.fastq -o low_cov_assembly/

# Contamination removal
# Remove contaminant contigs based on taxonomy
seqtk subseq scaffolds.fasta clean_contigs.txt > clean_assembly.fasta
```


## Clinical Applications

### Routine Surveillance
- Mantaining data quality control standards
- Rapid species identification and typing

## Looking Ahead

**Day 4 Preview**: 
- Genome assembly and
- Genome quality assessment
- Genome annotation with Prokka


---

**Key Learning Outcome**: Quality genomes to increase our confidence in our characterization capabilities essential for clinical genomics and public health surveillance.
