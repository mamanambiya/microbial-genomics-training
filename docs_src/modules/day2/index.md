# Day 2: Introduction to Commandline, High Performance Computing, & Quality Control

**Date**: September 2, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Command line proficiency, High Performance Computing introduction, M. tuberculosis genomics

## Overview

Day 2 focuses on building strong command line skills essential for bioinformatics work, introduces the Ilifu high-performance computing infrastructure, and includes a special guest lecture on M. tuberculosis and co-infection. This day provides the computational foundation needed for all subsequent genomic analyses in the course.

## Learning Objectives

By the end of Day 2, you will be able to:

- Master essential Unix/Linux command line operations for bioinformatics workflows
- Connect to and navigate the Ilifu high-performance computing cluster
- Submit and manage jobs using the SLURM scheduler
- Understand resource allocation and job queue management on HPC systems
- Apply command line skills to pathogen genomics data processing
- Understand M. tuberculosis genomics and co-infection patterns

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Introduction to command line interface* | [Practical](https://github.com/Arash-Iranzadeh/Microbial-Genomics/blob/main/scripts/practice_unix_commands.sh) | Arash Iranzadeh |
| **11:00** | *Introduction High Performance Computing (HPC) – Ilifu* | [Notes](../../day2/hpc-ilifu-training.md) • [Practical 1](../../day2/slurm-practical-tutorial.md) • [Practical 2](../../day2/unix-commands-pathogen-examples.md) | Mamana Mbiyavanga |
| **11:30** | **Break** | | |
| **12:00** | *Guest talk: MtB and co-infection* | | Bethlehem Adnew |

## Key Topics

### 1. Command Line Interface Fundamentals
- Unix/Linux file system navigation
- Essential commands for bioinformatics (grep, awk, sed)
- File manipulation and text processing
- Pipes and command chaining
- Working with compressed files (gzip, tar)
- Shell scripting basics for automation

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

### 4. Practical Applications in Pathogen Genomics
- Processing FASTQ files on HPC
- Running bioinformatics tools in parallel
- Managing large-scale genomic datasets
- Automating analysis pipelines
- Data transfer and storage strategies

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
- **Bash shell** - Command line interface and scripting
- **GNU coreutils** - Essential Unix utilities (ls, cd, grep, etc.)
- **Text processing** - awk, sed, cut, sort, uniq
- **File compression** - gzip, tar, zip
- **tmux/screen** - Terminal session management
- **rsync/scp** - File transfer and synchronization

## Hands-on Exercises

### Exercise 1: Command Line Fundamentals (90 minutes)
Master essential Unix commands for bioinformatics through practical exercises.

```bash
# Navigate file systems and manipulate files
cd ~/data
ls -la
mkdir analysis_output

# Process text files with Unix tools
grep "^>" sequences.fasta | wc -l  # Count sequences
cat sample.fastq | head -20         # View file contents

# Work with compressed files
gzip large_file.txt
gunzip -c compressed.gz | head

# Use pipes and redirection
cat data.txt | sort | uniq > unique_values.txt
```

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

## Key Concepts

### Command Line Essentials
- **File system navigation**: Understanding directory structure and paths
- **Text processing**: Using grep, sed, awk for data manipulation
- **Pipes and redirection**: Chaining commands for complex operations
- **Shell scripting**: Automating repetitive tasks
- **Regular expressions**: Pattern matching in bioinformatics data

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

### Individual Tasks
- Successfully connect to Ilifu HPC system
- Navigate Unix file system and manipulate files
- Write and submit a SLURM batch script
- Monitor job status and retrieve results
- Complete command line exercises for pathogen genomics data

### Group Discussion
- Share command line tips and tricks
- Discuss HPC resource management strategies
- Troubleshoot connection and job submission issues
- Compare different approaches to batch processing

## Common Challenges

### Command Line Challenges
```bash
# Permission denied errors
chmod +x script.sh    # Make script executable
ls -la                # Check file permissions

# Path issues
echo $PATH            # Check current PATH
export PATH=$PATH:/new/path  # Add to PATH
```

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

### Command Line Resources
- [Unix for Bioinformatics](https://bioinformatics.ucdavis.edu/research-computing/documentation/unix-basics/)
- [Bash Scripting Guide](https://www.gnu.org/software/bash/manual/)
- [Command Line for Genomics](https://datacarpentry.org/shell-genomics/)

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
- Quality checking and control with FastQC
- Species identification using Kraken2
- Genome assembly and quality assessment
- Genome annotation with Prokka

---

**Key Learning Outcome**: Mastery of command line operations and HPC infrastructure usage provides the essential computational foundation for all subsequent genomic analyses in the course. Understanding of M. tuberculosis genomics adds clinical context to the technical skills.