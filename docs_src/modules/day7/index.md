# Day 7: Nextflow

**Date**: September 9, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Reproducible workflows with Nextflow and pipeline development

## Overview

Day 7 introduces Nextflow as a workflow management system for creating reproducible, scalable bioinformatics pipelines. We'll explore the nf-core community and develop a complete pipeline for quality control, de novo assembly, quality assessment, and annotation.

## Learning Objectives

By the end of Day 7, you will be able to:

- Understand the principles of reproducible computational workflows
- Use Nextflow for bioinformatics pipeline development
- Leverage nf-core community pipelines and resources
- Develop custom Nextflow processes and workflows
- Implement best practices for pipeline development
- Create modular, reusable pipeline components

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Reproducible workflows with Nextflow and nf-core* | [Slides](#) • [Nextflow Guide](#) | Mamana Mbiyavanga |
| **10:30** | *Developing a Nextflow pipeline for QC, de novo assembly, quality assessment and annotation* | [Tutorial](#) • [Template](#) | Mamana Mbiyavanga |
| **11:30** | **Break** | | |
| **12:00** | *Developing a Nextflow pipeline for QC, de novo assembly, quality assessment and annotation* | [Practical](#) • [Scripts](#) | Mamana Mbiyavanga |

## Key Topics

### 1. Workflow Management Systems
- Challenges in bioinformatics workflows
- Benefits of workflow managers vs shell scripts
- Comparison of workflow systems (Nextflow, Snakemake, CWL)
- Reproducibility and scalability considerations

### 2. Nextflow Fundamentals
- Nextflow language syntax and concepts
- Processes, channels, and operators
- Configuration files and profiles
- Container integration (Docker/Singularity)

### 3. nf-core Community
- Pre-built, community-maintained pipelines
- Pipeline standards and best practices
- Tools for pipeline development and testing
- Contributing to the nf-core community

### 4. Pipeline Development
- Modular pipeline design principles
- Error handling and debugging
- Testing and validation strategies
- Documentation and metadata

## Tools and Resources

### Core Nextflow Tools
- **Nextflow** - Workflow management system
- **nf-core tools** - Pipeline development utilities
- **Docker/Singularity** - Container technologies
- **Conda** - Environment management

### Development Environment
- **Visual Studio Code** with Nextflow extension
- **Git** for version control
- **GitHub Actions** for continuous integration
- **nf-test** for pipeline testing

## Hands-on Exercises

### Exercise 1: Nextflow Basics (45 minutes)
Learn fundamental Nextflow concepts through simple examples.

```nextflow
#!/usr/bin/env nextflow

// Simple process definition
process FASTQC {
    container 'biocontainers/fastqc:v0.11.9_cv8'
    
    input:
    path reads
    
    output:
    path "*.html"
    path "*.zip"
    
    script:
    """
    fastqc ${reads}
    """
}

// Workflow definition
workflow {
    // Create channel from input files
    reads_ch = Channel.fromPath(params.reads)
    
    // Run FastQC process
    FASTQC(reads_ch)
}
```

### Exercise 2: Building QC-Assembly Pipeline (90 minutes)
Develop a complete pipeline for microbial genome analysis.

```nextflow
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process definitions
process FASTQC {
    tag "$sample_id"
    container 'biocontainers/fastqc:v0.11.9_cv8'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_fastqc.{html,zip}")
    
    script:
    """
    fastqc ${reads}
    """
}

process TRIMMOMATIC {
    tag "$sample_id"
    container 'biocontainers/trimmomatic:v0.39-1'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_paired.fastq.gz")
    
    script:
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \
        ${sample_id}_R1_paired.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \
        ${sample_id}_R2_paired.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process SPADES {
    tag "$sample_id"
    container 'biocontainers/spades:v3.15.3-1'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta")
    
    script:
    """
    spades.py --careful -1 ${reads[0]} -2 ${reads[1]} -o spades_output
    cp spades_output/scaffolds.fasta ${sample_id}_assembly.fasta
    """
}

process QUAST {
    tag "$sample_id"
    container 'biocontainers/quast:v5.0.2-1'
    publishDir "${params.outdir}/quast", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "${sample_id}_quast"
    
    script:
    """
    quast.py ${assembly} -o ${sample_id}_quast
    """
}

process PROKKA {
    tag "$sample_id"
    container 'biocontainers/prokka:v1.14.6-1'
    publishDir "${params.outdir}/annotations", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "${sample_id}_annotation"
    
    script:
    """
    prokka ${assembly} --outdir ${sample_id}_annotation --prefix ${sample_id} \
        --kingdom Bacteria --cpus ${task.cpus}
    """
}

// Workflow definition
workflow {
    // Input data
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    
    // QC and trimming
    FASTQC(reads_ch)
    TRIMMOMATIC(reads_ch)
    
    // Assembly
    SPADES(TRIMMOMATIC.out)
    
    // Quality assessment
    QUAST(SPADES.out)
    
    // Annotation
    PROKKA(SPADES.out)
}
```

### Exercise 3: Pipeline Configuration (30 minutes)
Create configuration files for different execution environments.

```groovy
// nextflow.config
params {
    // Input parameters
    reads = "data/*_R{1,2}_001.fastq.gz"
    outdir = "results"
    
    // Resource defaults
    max_cpus = 16
    max_memory = '64.GB'
    max_time = '24.h'
}

// Profile configurations
profiles {
    local {
        process.executor = 'local'
        process.cpus = 4
        process.memory = '8.GB'
        docker.enabled = true
    }
    
    slurm {
        process.executor = 'slurm'
        process.queue = 'batch'
        process.cpus = 16
        process.memory = '32.GB'
        singularity.enabled = true
    }
    
    conda {
        conda.enabled = true
        process.conda = 'envs/microbial_genomics.yml'
    }
}

// Process-specific configurations
process {
    withName: SPADES {
        cpus = 8
        memory = '16.GB'
        time = '4.h'
    }
    
    withName: PROKKA {
        cpus = 4
        memory = '8.GB'
        time = '2.h'
    }
}
```

## Key Concepts

### Nextflow Language Features
- **Processes**: Encapsulated tasks with inputs, outputs, and script
- **Channels**: Asynchronous data streams connecting processes
- **Operators**: Transform and manipulate channel data
- **Workflows**: Orchestrate process execution and data flow

### Best Practices
- **Modularity**: Break complex workflows into reusable components
- **Parameterization**: Make pipelines configurable and flexible
- **Error handling**: Implement robust error recovery mechanisms
- **Testing**: Validate pipeline components and outputs
- **Documentation**: Clear usage instructions and examples

### Configuration Management
```groovy
// Environment-specific settings
profiles {
    standard {
        process.container = 'ubuntu:20.04'
    }
    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
    }
}
```

## Pipeline Development Workflow

### 1. Planning Phase
- Define pipeline objectives and scope
- Identify required tools and dependencies
- Design modular architecture
- Plan testing and validation strategy

### 2. Development Phase
- Create process definitions for each step
- Implement workflow logic and data flow
- Add configuration and parameterization
- Include error handling and logging

### 3. Testing Phase
- Test individual processes with sample data
- Validate complete workflow execution
- Performance testing and optimization
- Edge case and error condition testing

### 4. Deployment Phase
- Containerize tools and dependencies
- Create user documentation
- Set up continuous integration
- Community sharing and feedback

## Assessment Activities

### Individual Tasks
- Create simple Nextflow process for bioinformatics tool
- Develop multi-step pipeline with proper data flow
- Configure pipeline for different execution environments
- Test pipeline with provided datasets

### Group Exercise
- Collaborate on pipeline feature development
- Review and optimize pipeline code
- Share experiences with different configuration profiles
- Troubleshoot common pipeline issues

## Common Challenges

### Resource Management
```groovy
// Dynamic resource allocation
process ASSEMBLY {
    memory { 4.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    
    script:
    """
    spades.py --memory ${task.memory.toGiga()} -t ${task.cpus} ...
    """
}
```

### Data Staging Issues
```groovy
// Proper input handling
process ANALYSIS {
    stageInMode 'copy'  // Ensure file availability
    scratch true        // Use local scratch space
    
    input:
    path(large_file)
    
    script:
    """
    # Process large file locally
    analysis_tool ${large_file}
    """
}
```

### Container Compatibility
```groovy
// Container-specific configurations
process {
    withName: SPECIAL_TOOL {
        container = 'custom/special-tool:v1.0'
        containerOptions = '--user root'
    }
}
```

## Resources

### Documentation
- [Nextflow Documentation](https://nextflow.io/docs/latest/)
- [nf-core Website](https://nf-co.re/)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)

### Training Materials
- [Nextflow Training](https://training.nextflow.io/)
- [nf-core Bytesize Talks](https://nf-co.re/events/bytesize)
- [Seqera Community Forum](https://community.seqera.io/)

### Development Tools
- [nf-core tools](https://nf-co.re/tools)
- [Nextflow VS Code Extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow)
- [nf-test Framework](https://nf-test.com/)

## Looking Ahead

**Day 8 Preview**: Pipeline Development continues with:
- Advanced pipeline optimization techniques
- Comprehensive testing and validation
- Performance monitoring and profiling
- Production deployment strategies

---

**Key Learning Outcome**: Nextflow enables development of reproducible, scalable bioinformatics workflows that can run consistently across different computing environments, from laptops to high-performance clusters.