# Day 6: Nextflow Foundations & Core Concepts

**Date**: September 8, 2025
**Duration**: 09:00-13:00 CAT
**Focus**: Introduction to workflow management, Nextflow fundamentals, and first pipelines

## Overview

Day 6 introduces participants to workflow management systems and Nextflow fundamentals. This comprehensive session covers the theoretical foundations of reproducible workflows, core Nextflow concepts, and hands-on development of basic pipelines. Participants will understand why workflow management is crucial for bioinformatics and gain practical experience with Nextflow's core components.

## Learning Objectives

By the end of Day 6, you will be able to:

- Understand the challenges in bioinformatics reproducibility and benefits of workflow management systems
- Explain Nextflow's core features and architecture
- Identify the main components of a Nextflow script (processes, channels, workflows)
- Write and execute basic Nextflow processes and workflows
- Use channels to manage data flow between processes
- Configure Nextflow for different execution environments
- Debug common Nextflow issues and understand error messages
- Apply best practices for pipeline development

## Schedule

| Time (CAT) | Topic | Duration | Trainer |
|------------|-------|----------|---------|
| **09:00** | *Foundation: Command line, Git, Containers review* | 30 min | Mamana Mbiyavanga |
| **09:30** | *Introduction to Workflow Management Systems* | 45 min | Mamana Mbiyavanga |
| **10:15** | *Nextflow Basics: Core concepts and first pipelines* | 75 min | Mamana Mbiyavanga |
| **11:30** | **Break** | 15 min | |
| **11:45** | *Hands-on: Building your first Nextflow pipeline* | 75 min | Mamana Mbiyavanga |

## Key Topics

### 1. Foundation Review (30 minutes)

- Command line proficiency check
- Git basics and version control concepts
- Container technologies (Docker/Singularity) overview
- Setting up the development environment

### 2. Introduction to Workflow Management (45 minutes)

- The challenge of complex genomics analyses
- Problems with traditional scripting approaches
- Benefits of workflow management systems
- Nextflow vs other systems (Snakemake, CWL, WDL)
- Reproducibility, portability, and scalability

### 3. Nextflow Core Concepts (75 minutes)

- Nextflow architecture and execution model
- Processes: encapsulated tasks with inputs, outputs, and scripts
- Channels: asynchronous data streams connecting processes
- Workflows: orchestrating process execution and data flow
- The work directory structure and caching mechanism
- Executors and execution platforms

### 4. Hands-on Pipeline Development (75 minutes)

- Writing your first Nextflow process
- Creating channels and managing data flow
- Building a simple QC workflow
- Testing and debugging pipelines
- Understanding the work directory

## Tools and Software

### Core Requirements

- **Nextflow** (version 20.10.0 or later) - Workflow orchestration system
- **Java** (version 11 or later) - Required for Nextflow execution
- **Docker** or **Singularity** - Container platforms for reproducibility
- **Text editor** - VS Code with Nextflow extension recommended

### Bioinformatics Tools

- **FastQC** - Read quality control assessment
- **MultiQC** - Aggregate quality control reports
- **Trimmomatic** - Read trimming and filtering
- **SPAdes** - Genome assembly (for later exercises)

### Development Environment

- **Git** - Version control system
- **GitHub account** - For sharing and collaboration
- **Terminal/Command line** - For running Nextflow commands

## Part 1: The Challenge of Complex Genomics Analyses

### Why Workflow Management Matters

Consider analyzing 100 bacterial genomes without workflow management:

```bash
# Manual approach - tedious and error-prone
for sample in sample1 sample2 sample3 ... sample100; do
    fastqc ${sample}_R1.fastq ${sample}_R2.fastq
    if [ $? -ne 0 ]; then echo "FastQC failed"; exit 1; fi

    trimmomatic PE ${sample}_R1.fastq ${sample}_R2.fastq \
        ${sample}_R1_trimmed.fastq ${sample}_R1_unpaired.fastq \
        ${sample}_R2_trimmed.fastq ${sample}_R2_unpaired.fastq \
        SLIDINGWINDOW:4:20
    if [ $? -ne 0 ]; then echo "Trimming failed"; exit 1; fi

    spades.py -1 ${sample}_R1_trimmed.fastq -2 ${sample}_R2_trimmed.fastq \
        -o ${sample}_assembly
    if [ $? -ne 0 ]; then echo "Assembly failed"; exit 1; fi

    # What if step 3 fails for sample 67?
    # How do you restart from where it failed?
    # How do you run samples in parallel efficiently?
    # How do you ensure reproducibility across different systems?
done
```

### Problems with Traditional Approaches

1. **Error Handling**: Manual error checking is verbose and error-prone
2. **Parallelization**: Difficult to efficiently use multiple cores/nodes
3. **Resumability**: No easy way to restart from failed steps
4. **Reproducibility**: Hard to ensure same results across different systems
5. **Scalability**: Doesn't scale well from laptop to HPC to cloud
6. **Dependency Management**: Software installation and version conflicts
7. **Resource Management**: No automatic optimization of CPU/memory usage

### The Workflow Management Solution

With Nextflow, you define the workflow once and it handles:

- **Automatic parallelization** of all 100 samples
- **Intelligent resource management** (memory, CPUs)
- **Automatic retry** of failed tasks with different resources
- **Resume capability** from the last successful step
- **Container integration** for reproducibility
- **Detailed execution reports** and monitoring
- **Platform portability** (laptop → HPC → cloud)

## Part 2: Nextflow Architecture and Core Concepts

### Nextflow's Key Components

#### 1. **Nextflow Engine**

The core runtime that interprets and executes your pipeline:

- Parses the workflow script
- Manages task scheduling and execution
- Handles data flow between processes
- Provides caching and resume capabilities

#### 2. **Work Directory**

Where Nextflow stores intermediate files and task execution:

```
work/
├── 12/
│   └── 3456789abcdef.../
│       ├── .command.sh      # The actual script executed
│       ├── .command.run     # Wrapper script
│       ├── .command.out     # Standard output
│       ├── .command.err     # Standard error
│       ├── .command.log     # Execution log
│       ├── .exitcode       # Exit status
│       └── input_file.fastq # Staged input files
└── ab/
    └── cdef123456789.../
        └── ...
```

#### 3. **Executors**

Interface with different computing platforms:

- **Local**: Run on your laptop/desktop
- **SLURM**: Submit jobs to HPC clusters
- **AWS Batch**: Execute on Amazon cloud
- **Kubernetes**: Run on container orchestration platforms

### Core Nextflow Components

#### **Process**

A process defines a task to be executed. It's the basic building block of a Nextflow pipeline:

```nextflow
process FASTQC {
    // Process directives
    tag "$sample_id"
    container 'biocontainers/fastqc:v0.11.9_cv8'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.{html,zip}"), emit: reports

    script:
    """
    fastqc ${reads}
    """
}
```

**Key Elements:**

- **Directives**: Configure how the process runs (container, resources, etc.)
- **Input**: Define what data the process expects
- **Output**: Define what data the process produces
- **Script**: The actual command(s) to execute

#### **Channel**

Channels are asynchronous data streams that connect processes:

```nextflow
// Create channel from file pairs
reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")

// Create channel from a list
samples_ch = Channel.from(['sample1', 'sample2', 'sample3'])

// Create channel from a file
reference_ch = Channel.fromPath("reference.fasta")
```

**Channel Types:**

- **Queue channels**: Can be consumed only once
- **Value channels**: Can be consumed multiple times
- **File channels**: Handle file paths and staging

#### **Workflow**

The workflow block orchestrates process execution:

```nextflow
workflow {
    // Define input channels
    reads_ch = Channel.fromFilePairs(params.reads)

    // Execute processes
    FASTQC(reads_ch)

    // Chain processes together
    TRIMMOMATIC(reads_ch)
    SPADES(TRIMMOMATIC.out.trimmed)

    // Access outputs
    FASTQC.out.reports.view()
}
```

## Part 3: Hands-on Exercises

### Exercise 1: Installation and Setup (15 minutes)

**Objective**: Install Nextflow and verify the environment

```bash
# Check Java version (must be 11 or later)
java -version

# Install Nextflow
curl -s https://get.nextflow.io | bash

# Make executable and add to PATH
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow info

# Test with hello world
nextflow run hello
```

### Exercise 2: Your First Nextflow Script (30 minutes)

**Objective**: Create and run a simple Nextflow pipeline

Create a file called `word_count.nf`:

```nextflow
#!/usr/bin/env nextflow

// Pipeline parameters
params.input = "data/yeast/reads/ref1_1.fq.gz"

// Input channel
input_ch = Channel.fromPath(params.input)

// Main workflow
workflow {
    NUM_LINES(input_ch)
    NUM_LINES.out.view()
}

// Process definition
process NUM_LINES {
    input:
    path read

    output:
    stdout

    script:
    """
    printf '${read}\\t'
    gunzip -c ${read} | wc -l
    """
}
```

**Run the pipeline:**

```bash
# Create test data directory
mkdir -p data/yeast/reads

# Download test file (or create a dummy file)
echo -e "@read1\nACGT\n+\nIIII" | gzip > data/yeast/reads/ref1_1.fq.gz

# Run the pipeline
nextflow run word_count.nf

# Examine the work directory
ls -la work/
```

### Exercise 3: Understanding Channels (20 minutes)

**Objective**: Learn different ways to create and manipulate channels

Create `channel_examples.nf`:

```nextflow
#!/usr/bin/env nextflow

workflow {
    // Channel from file pairs
    reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    reads_ch.view { sample, files -> "Sample: $sample, Files: $files" }

    // Channel from list
    samples_ch = Channel.from(['sample1', 'sample2', 'sample3'])
    samples_ch.view { "Processing: $it" }

    // Channel from path pattern
    ref_ch = Channel.fromPath("*.fasta")
    ref_ch.view { "Reference: $it" }
}
```

Save and push to GitHub:
```bash
git add README.md
git commit -m "Add README documentation"
git push
```

## Key Concepts

### Version Control Benefits
- **Track Changes**: See what changed, when, and why
- **Backup**: Your code is safe on GitHub
- **Collaboration**: Work with others easily
- **Reproducibility**: Others can use your exact pipeline version

### Git Workflow
```
Working Directory → Staging Area → Local Repository → Remote Repository
     (edit)           (git add)      (git commit)        (git push)
```

### Commit Best Practices
- Write clear, descriptive messages
- Commit logical units of change
- Commit frequently
- Don't commit generated files (use .gitignore)

## Common Issues and Solutions

### Git Configuration
```bash
# Set up your identity (first time only)
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### Authentication Issues
```bash
# If password authentication fails, use personal access token:
# GitHub → Settings → Developer settings → Personal access tokens
# Generate new token with 'repo' permissions
# Use token instead of password when pushing
```

### Common Git Commands
| Command | Purpose | Example |
|---------|---------|---------|
| `git status` | Check repository state | `git status` |
| `git diff` | See changes | `git diff main.nf` |
| `git add .` | Stage all changes | `git add .` |
| `git commit -m` | Save changes | `git commit -m "Fix bug"` |
| `git push` | Upload to GitHub | `git push origin main` |
| `git pull` | Download from GitHub | `git pull origin main` |

## Assessment Activities

### Individual Tasks
- Successfully complete and run the Nextflow pipeline
- Initialize a Git repository for your pipeline
- Create a GitHub account and repository
- Push your pipeline to GitHub
- Create a comprehensive README file

### Group Discussion
- Share GitHub repository URLs with the class
- Discuss pipeline design choices
- Review each other's README files
- Troubleshoot Git/GitHub issues together

## Resources

### Nextflow Resources
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nf-core pipelines](https://nf-co.re/)

### Git/GitHub Resources
- [GitHub Guides](https://guides.github.com/)
- [Git Handbook](https://guides.github.com/introduction/git-handbook/)
- [GitHub Learning Lab](https://lab.github.com/)
- [Pro Git Book (free)](https://git-scm.com/book)

## Looking Ahead

**Day 7 Preview**: Applied Genomics & Advanced Topics

- MTB analysis pipeline development
- Genome assembly workflows
- Advanced Nextflow features and optimization
- Pipeline deployment and best practices

### Exercise 4: Building a QC Process (30 minutes)

**Objective**: Create a real bioinformatics process

Create `fastqc_pipeline.nf`:

```nextflow
#!/usr/bin/env nextflow

// Parameters
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = "results"

// Main workflow
workflow {
    // Create channel from paired reads
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Run FastQC
    FASTQC(reads_ch)

    // View results
    FASTQC.out.view { sample, reports ->
        "FastQC completed for $sample: $reports"
    }
}

// FastQC process
process FASTQC {
    tag "$sample_id"
    container 'biocontainers/fastqc:v0.11.9_cv8'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.{html,zip}")

    script:
    """
    fastqc ${reads}
    """
}
```

**Test the pipeline:**

```bash
# Create test data
mkdir -p data
echo -e "@read1\nACGTACGT\n+\nIIIIIIII" > data/sample1_R1.fastq
echo -e "@read1\nTGCATGCA\n+\nIIIIIIII" > data/sample1_R2.fastq

# Run pipeline
nextflow run fastqc_pipeline.nf

# Check results
ls -la results/fastqc/
```

## Troubleshooting Guide

### Installation Issues

```bash
# Java version problems
java -version  # Must be 11 or later

# Nextflow not found
echo $PATH
which nextflow

# Permission issues
chmod +x nextflow
```

### Pipeline Debugging

```bash
# Verbose output
nextflow run pipeline.nf -with-trace -with-report -with-timeline

# Check work directory
ls -la work/

# Resume from failure
nextflow run pipeline.nf -resume
```

---

**Key Learning Outcome**: Understanding workflow management fundamentals and Nextflow core concepts provides the foundation for building reproducible, scalable bioinformatics pipelines.
