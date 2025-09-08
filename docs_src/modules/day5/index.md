# Day 5: Tracking Threats: Genomic Detection of AMR, Virulence, and Plasmid Mobility

**Date**: September 5, 2025
**Duration**: 09:00-13:00 CAT
**Focus**: Genome quality and functional gene annotation fundamentals, AMR and virulence factors and plasmid detection

## Overview

Day 5 introduces Nextflow, a powerful workflow management system for creating reproducible and scalable bioinformatics pipelines. We'll explore the fundamentals of Nextflow, the nf-core community standards, and begin developing a pipeline for genomic analysis including QC, assembly, quality assessment, and annotation.

## Learning Objectives

By the end of Day 5, you will be able to:

- Understand the principles of reproducible computational workflows
- Write basic Nextflow scripts with processes and channels
- Utilize nf-core tools and community pipelines
- Design workflow architecture for genomic analysis
- Implement data flow using Nextflow channels
- Begin developing a pipeline for QC, assembly, and annotation

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Reproducible workflows with Nextflow and nf-core* | | Mamana Mbiyavanga |
| **10:30** | *Developing a Nextflow pipeline for QC, de novo assembly, quality assessment and annotation* | | Mamana Mbiyavanga |
| **11:30** | **Break** | | |
| **12:00** | *Developing a Nextflow pipeline for QC, de novo assembly, quality assessment and annotation* | | Mamana Mbiyavanga |

## Key Topics

### 1. Introduction to Workflow Management
- Challenges in bioinformatics reproducibility
- Benefits of workflow management systems
- Nextflow vs other workflow systems (Snakemake, CWL, WDL)
- Container technologies (Docker, Singularity)

### 2. Nextflow Fundamentals
- Nextflow architecture and concepts
- Processes, channels, and operators
- Configuration files and profiles
- Resource management and executors
- Error handling and resume capabilities

### 3. nf-core Community and Standards
- nf-core pipeline structure
- Community guidelines and best practices
- Using nf-core tools
- Available nf-core pipelines for genomics
- Contributing to nf-core

### 4. Building a Genomic Analysis Pipeline
- Pipeline design and planning
- Implementing QC processes (FastQC, MultiQC)
- Assembly process integration (SPAdes)
- Quality assessment steps (QUAST)
- Annotation process (Prokka)

### 5. Nextflow Scripting
- Writing process definitions
- Channel operations and data flow
- Parameter handling
- Conditional execution
- Module organization

## Tools and Software

### Workflow Management
- **Nextflow** - Workflow orchestration system
- **nf-core tools** - Pipeline development framework
- **Tower** - Workflow monitoring platform

### Containerization
- **Docker** - Container platform
- **Singularity** - HPC-friendly containers
- **Conda** - Package management

### Pipeline Components
- **FastQC** - Read quality control
- **MultiQC** - Aggregate reporting
- **SPAdes** - Genome assembly
- **QUAST** - Assembly assessment
- **Prokka** - Genome annotation

## Hands-on Exercises

### Exercise 1: First Nextflow Script (30 minutes)
Create and run a simple Nextflow pipeline.

```groovy
#!/usr/bin/env nextflow

// Define parameters
params.input = "data/*.fastq"
params.outdir = "results"

// Create a channel from input files
Channel
    .fromPath(params.input)
    .set { fastq_ch }

// Define a process
process countReads {
    input:
    path fastq from fastq_ch
    
    output:
    path "*.count" into counts_ch
    
    script:
    """
    echo "Processing ${fastq}"
    wc -l ${fastq} > ${fastq.baseName}.count
    """
}

// View the results
counts_ch.view()
```

### Exercise 2: Building a QC Pipeline (60 minutes)
Implement quality control with FastQC and MultiQC.

```groovy
process fastqc {
    container 'biocontainers/fastqc:v0.11.9'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{zip,html}" into fastqc_results
    
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

process multiqc {
    publishDir params.outdir, mode: 'copy'
    container 'ewels/multiqc:latest'
    
    input:
    path '*' from fastqc_results.collect()
    
    output:
    path 'multiqc_report.html'
    
    script:
    """
    multiqc .
    """
}
```

### Exercise 3: Integrating Assembly (90 minutes)
Add genome assembly to the pipeline.

```groovy
process spades_assembly {
    container 'staphb/spades:latest'
    cpus 4
    memory '8 GB'
    
    input:
    tuple val(sample_id), path(reads1), path(reads2)
    
    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")
    
    script:
    """
    spades.py \
        -1 ${reads1} \
        -2 ${reads2} \
        -o spades_output \
        -t ${task.cpus} \
        --careful
    
    cp spades_output/contigs.fasta ${sample_id}_contigs.fasta
    """
}
```

## Key Concepts

### Workflow Principles
- **Reproducibility**: Same input → same output
- **Portability**: Run anywhere (laptop, HPC, cloud)
- **Scalability**: Handle any data volume
- **Resumability**: Restart from failure points

### Nextflow Components
| Component | Description | Example |
|-----------|-------------|---------|
| Process | Computational step | `process fastqc { ... }` |
| Channel | Data flow connection | `Channel.fromPath()` |
| Operator | Channel transformation | `.map()`, `.filter()` |
| Directive | Process configuration | `cpus 4` |

### Best Practices
1. **Use containers**: Ensure environment reproducibility
2. **Parameterize everything**: Make pipelines flexible
3. **Version control**: Track pipeline changes
4. **Document thoroughly**: Help users and future self
5. **Test incrementally**: Build and test step by step

## Assessment Activities

### Individual Tasks
- Create a basic Nextflow script with at least 2 processes
- Successfully run a pipeline with test data
- Modify pipeline parameters and observe changes
- Debug a pipeline with intentional errors
- Document pipeline usage

### Group Discussion
- Compare Nextflow with traditional shell scripting
- Discuss reproducibility challenges and solutions
- Share pipeline design strategies
- Explore nf-core pipeline catalog

## Common Challenges

### Installation Issues
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
./nextflow run hello

# Set up environment
export PATH=$PATH:$PWD
export NXF_VER=23.10.0
```

### Channel Operations
```groovy
// Common channel patterns
Channel
    .fromFilePairs(params.reads)
    .ifEmpty { error "No read files found!" }
    .set { read_pairs_ch }

// Combining channels
fastqc_ch
    .join(assembly_ch)
    .map { sample, qc, assembly -> 
        [sample, qc, assembly]
    }
```

### Resource Management
```groovy
process memory_intensive {
    memory { 2.GB * task.attempt }
    maxRetries 3
    errorStrategy 'retry'
    
    script:
    """
    # Your command here
    """
}
```

## Resources

### Documentation
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nf-core Website](https://nf-co.re/)
- [Nextflow Training](https://training.nextflow.io/)

### Tutorials
- [Nextflow Tutorial](https://www.nextflow.io/docs/latest/getstarted.html)
- [nf-core Tutorials](https://nf-co.re/usage/tutorials)
- [Seqera Labs Training](https://training.seqera.io/)

### Community
- [Nextflow Slack](https://www.nextflow.io/slack-invite.html)
- [nf-core Slack](https://nf-co.re/join)
- [GitHub Discussions](https://github.com/nextflow-io/nextflow/discussions)

## Looking Ahead

**Day 6 Preview**: Nextflow Pipeline Development
- Continue building the genomic analysis pipeline
- Advanced Nextflow features and optimization
- Pipeline testing and validation
- Deployment strategies

---

**Key Learning Outcome**: Understanding workflow management principles and gaining hands-on experience with Nextflow enables creation of reproducible, scalable bioinformatics pipelines essential for modern genomic analysis.