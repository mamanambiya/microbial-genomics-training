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

### Why This Approach is "Tedious and Error-Prone"

**Major Problems with Traditional Shell Scripting:**

1. **No Parallelization**
   - Processes samples sequentially (one after another)
   - Wastes computational resources on multi-core systems
   - Takes unnecessarily long time

2. **Poor Error Recovery & Resumability**
   - If one sample fails, entire pipeline stops
   - No way to resume from failure point
   - Must restart from beginning
   - Manual error checking is verbose and error-prone

3. **Resource Management Issues**
   - No control over CPU/memory usage
   - Can overwhelm system or underutilize resources
   - No queue management for HPC systems
   - No automatic optimization of resource allocation

4. **Lack of Reproducibility**
   - Hard to track software versions
   - Environment dependencies not managed
   - Difficult to share and reproduce results across different systems
   - Software installation and version conflicts

5. **Poor Scalability**
   - Doesn't scale well from laptop to HPC to cloud
   - No automatic adaptation to different computing environments
   - Limited ability to handle varying data volumes

6. **Maintenance Nightmare**
   - Adding new steps requires modifying the entire script
   - Parameter changes need manual editing throughout
   - No modular design for reusable components
   - Difficult to test individual components

7. **No Progress Tracking**
   - Can't easily see which samples completed
   - No reporting or logging mechanisms
   - Difficult to debug failures
   - No visibility into pipeline performance

## The Workflow Management Solution

### Overview of Workflow Management Systems

**Workflow management systems (WMS)** are specialized programming languages and frameworks designed specifically to address the challenges of complex, multi-step computational pipelines. They provide a higher-level abstraction that automatically handles the tedious and error-prone aspects of traditional shell scripting.

#### How Workflow Management Systems Solve Traditional Problems:

**Automatic Parallelization**
- Analyze task dependencies and run independent steps simultaneously
- Efficiently utilize all available CPU cores and computing nodes
- Scale from single machines to massive HPC clusters and cloud environments

**Built-in Error Recovery**
- Automatic retry mechanisms for failed tasks
- Resume functionality to restart from failure points
- Intelligent caching to avoid re-running successful steps

**Resource Management**
- Automatic CPU and memory allocation based on task requirements
- Integration with job schedulers (SLURM, PBS, SGE)
- Dynamic scaling in cloud environments

**Reproducibility by Design**
- Container integration (Docker, Singularity) for consistent environments
- Version tracking for all software dependencies
- Portable execution across different computing platforms

**Progress Monitoring**
- Real-time pipeline execution tracking
- Detailed logging and reporting
- Performance metrics and resource usage statistics

**Modular Architecture**
- Reusable workflow components
- Easy parameter configuration
- Clean separation of logic and execution

### Comparison of Popular Workflow Languages

The bioinformatics community has developed several powerful workflow management systems, each with unique strengths and design philosophies:

#### 1. **Nextflow**
- **Language Base**: Groovy (JVM-based)
- **Philosophy**: Dataflow programming with reactive streams
- **Strengths**: Excellent parallelization, cloud-native, strong container support
- **Community**: Large bioinformatics community, nf-core ecosystem

#### 2. **Snakemake**
- **Language Base**: Python
- **Philosophy**: Rule-based workflow definition inspired by GNU Make
- **Strengths**: Pythonic syntax, excellent for Python developers, strong academic adoption
- **Community**: Very active in computational biology and data science

#### 3. **Common Workflow Language (CWL)**
- **Language Base**: YAML/JSON
- **Philosophy**: Vendor-neutral, standards-based approach
- **Strengths**: Platform independence, strong metadata support, scientific reproducibility focus
- **Community**: Broad industry and academic support across multiple domains

#### 4. **Workflow Description Language (WDL)**
- **Language Base**: Custom domain-specific language
- **Philosophy**: Human-readable workflow descriptions with strong typing
- **Strengths**: Excellent cloud integration, strong at Broad Institute and genomics centers
- **Community**: Strong in genomics, particularly for large-scale sequencing projects

### Feature Comparison Table

| Feature | Nextflow | Snakemake | CWL | WDL |
|---------|----------|-----------|-----|-----|
| **Syntax Base** | Groovy | Python | YAML/JSON | Custom DSL |
| **Learning Curve** | Moderate | Easy (for Python users) | Steep | Moderate |
| **Parallelization** | Excellent (automatic) | Excellent | Good | Excellent |
| **Container Support** | Native (Docker/Singularity) | Native | Native | Native |
| **Cloud Integration** | Excellent (AWS, GCP, Azure) | Good | Good | Excellent |
| **HPC Support** | Excellent (SLURM, PBS, etc.) | Excellent | Good | Good |
| **Resume Capability** | Excellent | Excellent | Limited | Good |
| **Community Size** | Large (bioinformatics) | Large (data science) | Medium | Medium |
| **Package Ecosystem** | nf-core (500+ pipelines) | Snakemake Wrappers | Limited | Limited |
| **Debugging Tools** | Good (Tower, reports) | Excellent | Limited | Good |
| **Best Use Cases** | Multi-omics, clinical pipelines | Data analysis, research | Standards compliance | Large-scale genomics |
| **Industry Adoption** | High (pharma, biotech) | High (academia) | Growing | High (genomics centers) |

### Simple Code Examples

Let's see how the same basic task - running FastQC on multiple samples - would be implemented in different workflow languages:

#### **Traditional Shell Script** (for comparison)
```bash
# Manual approach - sequential processing
for sample in sample1 sample2 sample3; do
    fastqc ${sample}_R1.fastq ${sample}_R2.fastq -o results/
    if [ $? -ne 0 ]; then echo "FastQC failed for $sample"; exit 1; fi
done
```

#### **Nextflow Implementation**
```groovy
#!/usr/bin/env nextflow

// Define input channel
Channel
    .fromFilePairs("data/*_{R1,R2}.fastq")
    .set { read_pairs_ch }

// FastQC process
process fastqc {
    container 'biocontainers/fastqc:v0.11.9'
    publishDir 'results/', mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc ${reads} -t ${task.cpus}
    """
}

// Run the workflow
workflow {
    fastqc(read_pairs_ch)
}
```

#### **Snakemake Implementation**
```python
# Snakefile
SAMPLES = ["sample1", "sample2", "sample3"]

rule all:
    input:
        expand("results/{sample}_{read}_fastqc.html",
               sample=SAMPLES, read=["R1", "R2"])

rule fastqc:
    input:
        "data/{sample}_{read}.fastq"
    output:
        html="results/{sample}_{read}_fastqc.html",
        zip="results/{sample}_{read}_fastqc.zip"
    container:
        "docker://biocontainers/fastqc:v0.11.9"
    shell:
        "fastqc {input} -o results/"
```

#### **CWL Implementation**
```yaml
# fastqc-workflow.cwl
cwlVersion: v1.2
class: Workflow

inputs:
  fastq_files:
    type: File[]

outputs:
  fastqc_reports:
    type: File[]
    outputSource: fastqc/html_report

steps:
  fastqc:
    run: fastqc-tool.cwl
    scatter: fastq_file
    in:
      fastq_file: fastq_files
    out: [html_report, zip_report]

# fastqc-tool.cwl
cwlVersion: v1.2
class: CommandLineTool

baseCommand: fastqc

inputs:
  fastq_file:
    type: File
    inputBinding:
      position: 1

outputs:
  html_report:
    type: File
    outputBinding:
      glob: "*_fastqc.html"
  zip_report:
    type: File
    outputBinding:
      glob: "*_fastqc.zip"

requirements:
  DockerRequirement:
    dockerPull: biocontainers/fastqc:v0.11.9
```

#### **Key Differences in Syntax:**

**Nextflow**: Uses Groovy syntax with channels for data flow, processes define computational steps
**Snakemake**: Python-based with rules that define input/output relationships, uses wildcards for pattern matching
**CWL**: YAML-based with explicit input/output definitions, requires separate tool and workflow files
**WDL**: Custom syntax with strong typing, task-based approach with explicit variable declarations

### Why Nextflow for This Course

This course focuses on **Nextflow** for several compelling reasons that make it particularly well-suited for microbial genomics workflows:

#### **1. Bioinformatics Community Adoption**
- **nf-core ecosystem**: Over 500 community-curated pipelines specifically for bioinformatics
- **Industry standard**: Widely adopted by pharmaceutical companies, biotech firms, and genomics centers
- **Active development**: Strong community support with regular updates and improvements

#### **2. Excellent Parallelization for Genomics**
- **Automatic scaling**: Seamlessly scales from single samples to thousands of genomes
- **Dataflow programming**: Natural fit for genomics pipelines with complex dependencies
- **Resource optimization**: Intelligent task scheduling maximizes computational efficiency

#### **3. Clinical and Production Ready**
- **Robust error handling**: Critical for clinical pipelines where reliability is essential
- **Comprehensive logging**: Detailed audit trails required for regulatory compliance
- **Resume capability**: Minimizes computational waste in long-running genomic analyses

#### **4. Multi-Platform Flexibility**
- **HPC integration**: Native support for SLURM, PBS, and other job schedulers common in genomics
- **Cloud-native**: Excellent support for AWS, Google Cloud, and Azure for scalable genomics
- **Container support**: Seamless Docker and Singularity integration for reproducible environments

#### **5. Microbial Genomics Specific Advantages**
- **Pathogen surveillance pipelines**: Many nf-core pipelines designed for bacterial genomics
- **AMR analysis workflows**: Established patterns for antimicrobial resistance detection
- **Outbreak investigation**: Scalable phylogenetic analysis capabilities
- **Metagenomics support**: Robust handling of complex metagenomic datasets

#### **6. Learning and Career Benefits**
- **Industry relevance**: Skills directly transferable to genomics industry positions
- **Growing demand**: Increasing adoption means more job opportunities
- **Comprehensive ecosystem**: Learning Nextflow provides access to hundreds of ready-to-use pipelines

The combination of these factors makes Nextflow an ideal choice for training the next generation of microbial genomics researchers and practitioners. Its balance of power, usability, and industry adoption ensures that skills learned in this course will be immediately applicable in real-world genomics applications.
```

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

```text
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
```text
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
