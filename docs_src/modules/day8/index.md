# Day 8: Comparative Genomics

**Date**: September 10, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Pan-genome analysis, phylogenetic inference, tree construction and visualization

## Overview

Day 8 focuses on comparative genomics approaches for understanding microbial diversity and evolutionary relationships. We'll explore pangenome analysis to understand core and accessory gene content, and phylogenomic methods to infer evolutionary relationships from genomic data, including SNP-based phylogeny and tree visualization techniques.

## Learning Objectives

By the end of Day 8, you will be able to:

- Understand pangenome concepts and perform core/accessory genome analysis
- Identify conserved and variable genomic regions across strains
- Construct phylogenetic trees from core genome SNPs
- Visualize and interpret phylogenomic relationships
- Apply comparative genomics to understand pathogen evolution
- Use tools like Roary, Panaroo, and IQ-TREE for comparative analysis

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Pangenomics* | | Arash Iranzadeh |
| **10:30** | *Phylogenomics: Inferring evolutionary relationships from core SNPs* | | Arash Iranzadeh |
| **11:30** | **Break** | | |
| **12:00** | *Phylogenomics: Tree construction and visualisation* | | Arash Iranzadeh |

## Key Topics

### 1. Advanced Pipeline Architecture
- Multi-sample processing strategies
- Conditional execution and branching
- Pipeline modularity and reusability
- Configuration management across environments

### 2. Testing and Validation
- Unit testing for individual processes
- Integration testing for complete workflows
- Continuous integration setup
- Regression testing strategies

### 3. Performance Optimization
- Resource allocation strategies
- Parallelization patterns
- Caching and resume functionality
- Profile-based optimization

### 4. Production Deployment
- Environment-specific configurations
- Error handling and retry strategies
- Logging and monitoring
- Version control and release management

## Advanced Tools

### Testing Frameworks
- **nf-test** - Modern testing framework for Nextflow
- **pytest-workflow** - Python-based workflow testing
- **Nextflow Tower** - Pipeline monitoring and management

### Development Tools
- **nf-core lint** - Code quality checking
- **pre-commit hooks** - Automated code validation
- **GitHub Actions** - Continuous integration
- **Nextflow plugins** - Extended functionality

## Hands-on Exercises

### Exercise 1: Multi-Sample Pipeline (75 minutes)
Develop a comprehensive pipeline that processes multiple samples in parallel.

```nextflow
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameter definitions with validation
params.input_dir = null
params.outdir = "results"
params.reference = null
params.min_quality = 20
params.threads = 4

// Input validation
if (!params.input_dir) {
    error "Please provide input directory with --input_dir"
}
if (!params.reference) {
    error "Please provide reference genome with --reference"
}

// Include modules from separate files
include { FASTQC } from './modules/fastqc.nf'
include { TRIMMOMATIC } from './modules/trimmomatic.nf'
include { BWA_MEM } from './modules/bwa.nf'
include { VARIANT_CALLING } from './modules/variants.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow VARIANT_ANALYSIS {
    take:
    reads_ch
    reference
    
    main:
    // Quality control
    FASTQC(reads_ch)
    
    // Read trimming
    TRIMMOMATIC(reads_ch)
    
    // Alignment
    BWA_MEM(TRIMMOMATIC.out.trimmed, reference)
    
    // Variant calling
    VARIANT_CALLING(BWA_MEM.out.bam, reference)
    
    // Aggregate QC
    qc_files = FASTQC.out.html.mix(TRIMMOMATIC.out.log).collect()
    MULTIQC(qc_files)
    
    emit:
    variants = VARIANT_CALLING.out.vcf
    reports = MULTIQC.out.html
}

// Main workflow
workflow {
    // Create channel from input directory
    reads_ch = Channel
        .fromFilePairs("${params.input_dir}/*_R{1,2}_001.fastq.gz")
        .ifEmpty { error "No input files found in ${params.input_dir}" }
    
    // Load reference
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    
    // Run analysis
    VARIANT_ANALYSIS(reads_ch, reference_ch)
    
    // Output summary
    VARIANT_ANALYSIS.out.variants.view { "Variants called for: ${it[0]}" }
}

workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    Results: ${params.outdir}
    """.stripIndent()
}
```

### Exercise 2: Pipeline Testing (60 minutes)
Implement comprehensive testing for the variant calling pipeline.

```yaml
# .github/workflows/test.yml
name: Pipeline Testing

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        nextflow: ['23.04.0', 'latest']
    
    steps:
    - uses: actions/checkout@v4
    
    - name: Setup Nextflow
      uses: nf-core/setup-nextflow@v1
      with:
        version: ${{ matrix.nextflow }}
    
    - name: Setup Test Data
      run: |
        mkdir test_data
        wget -O test_data/sample_R1.fastq.gz https://example.com/test_R1.fastq.gz
        wget -O test_data/sample_R2.fastq.gz https://example.com/test_R2.fastq.gz
        wget -O test_data/reference.fasta https://example.com/reference.fasta
    
    - name: Run Pipeline Tests
      run: |
        nextflow run main.nf \
          --input_dir test_data \
          --reference test_data/reference.fasta \
          --outdir test_results \
          -profile test,docker
    
    - name: Validate Outputs
      run: |
        # Check if expected output files exist
        test -f test_results/variants/*.vcf
        test -f test_results/multiqc/multiqc_report.html
        
        # Validate variant file format
        bcftools view test_results/variants/*.vcf | head -20
```

### Exercise 3: Performance Optimization (45 minutes)
Optimize pipeline performance through profiling and resource tuning.

```nextflow
// Performance optimized configuration
process {
    // Default resources
    cpus = 2
    memory = 4.GB
    time = '1.hour'
    
    // Process-specific optimization
    withName: BWA_MEM {
        cpus = { 8 * task.attempt }
        memory = { 16.GB * task.attempt }
        time = { 4.hour * task.attempt }
        errorStrategy = 'retry'
        maxRetries = 3
    }
    
    withName: VARIANT_CALLING {
        cpus = 4
        memory = { 8.GB + (2.GB * task.attempt) }
        time = { 2.hour * task.attempt }
        
        // Use faster local storage when available
        scratch = '/tmp'
    }
    
    withName: FASTQC {
        // Lightweight process can use minimal resources
        cpus = 1
        memory = 2.GB
        time = 30.min
    }
}

// Profile-specific optimizations
profiles {
    standard {
        process.executor = 'local'
        process.cpus = 2
        process.memory = '4 GB'
    }
    
    hpc {
        process.executor = 'slurm'
        process.queue = 'compute'
        process.clusterOptions = '--account=genomics --qos=normal'
        
        // Optimize for cluster environment
        process {
            withName: BWA_MEM {
                cpus = 16
                memory = 32.GB
                time = 2.hour
            }
        }
    }
    
    cloud {
        process.executor = 'awsbatch'
        aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'
        aws.region = 'us-west-2'
        
        // Cloud-specific optimizations
        process {
            withName: '.*' {
                container = 'your-ecr-repo/pipeline:latest'
            }
        }
    }
}
```

## Advanced Concepts

### Error Handling Strategies

```nextflow
process ROBUST_ASSEMBLY {
    errorStrategy 'retry'
    maxRetries 3
    
    // Dynamic resource allocation
    memory { 8.GB * task.attempt }
    cpus { 4 * task.attempt }
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    tuple val(sample_id), path("${sample_id}_assembly.log"), emit: log
    
    script:
    """
    # Log system information for debugging
    echo "Hostname: \$(hostname)" > ${sample_id}_assembly.log
    echo "Memory: ${task.memory}" >> ${sample_id}_assembly.log
    echo "CPUs: ${task.cpus}" >> ${sample_id}_assembly.log
    echo "Attempt: ${task.attempt}" >> ${sample_id}_assembly.log
    
    # Run assembly with error checking
    spades.py --careful -1 ${reads[0]} -2 ${reads[1]} \
        -o spades_out --threads ${task.cpus} --memory ${task.memory.toGiga()} \
        2>&1 | tee -a ${sample_id}_assembly.log
    
    if [ ! -f spades_out/scaffolds.fasta ]; then
        echo "ERROR: Assembly failed" >> ${sample_id}_assembly.log
        exit 1
    fi
    
    cp spades_out/scaffolds.fasta ${sample_id}_assembly.fasta
    """
}
```

### Conditional Workflows

```nextflow
workflow ADAPTIVE_ANALYSIS {
    take:
    samples
    
    main:
    // Initial QC
    FASTQC(samples)
    
    // Conditional trimming based on quality
    samples
        .join(FASTQC.out.stats)
        .branch { sample_id, reads, qc_stats ->
            high_quality: qc_stats.mean_quality > 30
            needs_trimming: qc_stats.mean_quality <= 30
        }
        .set { qc_branched }
    
    // Process high quality samples directly
    high_qual_samples = qc_branched.high_quality.map { sample_id, reads, stats -> [sample_id, reads] }
    
    // Trim lower quality samples
    TRIMMOMATIC(qc_branched.needs_trimming.map { sample_id, reads, stats -> [sample_id, reads] })
    
    // Combine processed samples
    all_samples = high_qual_samples.mix(TRIMMOMATIC.out.trimmed)
    
    // Continue with assembly
    SPADES(all_samples)
    
    emit:
    assemblies = SPADES.out.assembly
}
```

## Best Practices

### 1. Code Organization
- Use modules for reusable processes
- Implement clear naming conventions
- Document complex logic
- Version control all configurations

### 2. Resource Management
```groovy
// Dynamic resource allocation
process {
    withLabel: 'high_memory' {
        memory = { 32.GB * task.attempt }
        errorStrategy = 'retry'
        maxRetries = 2
    }
    
    withLabel: 'cpu_intensive' {
        cpus = { Math.min(16, task.attempt * 4) }
        time = { 4.hour * task.attempt }
    }
}
```

### 3. Monitoring and Debugging
```nextflow
// Enable comprehensive reporting
trace {
    enabled = true
    file = "${params.outdir}/trace.txt"
    fields = 'task_id,hash,native_id,process,tag,name,status,exit,module,container,cpus,time,disk,memory,attempt,submit,start,complete,duration,realtime,queue,%cpu,%mem,rss,vmem,peak_rss,peak_vmem,rchar,wchar,syscr,syscw,read_bytes,write_bytes,vol_ctxt,inv_ctxt'
}

report {
    enabled = true
    file = "${params.outdir}/report.html"
}

timeline {
    enabled = true
    file = "${params.outdir}/timeline.html"
}

dag {
    enabled = true
    file = "${params.outdir}/dag.html"
}
```

## Assessment Activities

### Individual Projects
- Optimize existing pipeline for specific use case
- Implement comprehensive error handling
- Create test suite for pipeline validation
- Deploy pipeline to different computing environment

### Group Collaboration
- Code review session for pipeline improvements
- Troubleshooting complex pipeline failures
- Sharing optimization strategies
- Planning production deployment

## Common Challenges

### Memory Management
```nextflow
// Handle large datasets efficiently
process LARGE_DATA_PROCESSING {
    memory { task.attempt < 3 ? 16.GB : 32.GB }
    time { 2.hour * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    
    // Use streaming where possible
    script:
    """
    # Process data in chunks to manage memory
    split -l 1000000 ${large_input} chunk_
    
    for chunk in chunk_*; do
        process_chunk.py \$chunk >> results.txt
        rm \$chunk  # Clean up as we go
    done
    """
}
```

### Workflow Resume Issues
```bash
# Best practices for resumable workflows
nextflow run pipeline.nf -resume -with-report report.html

# Clean resume when needed
nextflow clean -f
rm -rf work/
```

### Container Compatibility
```nextflow
process {
    withName: PROBLEMATIC_TOOL {
        container = 'custom/fixed-tool:v2.0'
        containerOptions = '--user root --privileged'
    }
}
```

## Production Deployment

### Environment Setup
```yaml
# production.config
process {
    executor = 'slurm'
    queue = 'production'
    
    // Production-level resource allocation
    cpus = 16
    memory = '64 GB'
    time = '12 hours'
    
    // Enhanced error handling
    errorStrategy = 'terminate'  // Fail fast in production
    maxRetries = 1
}

// Enable comprehensive logging
trace.enabled = true
report.enabled = true
timeline.enabled = true
```

### Monitoring Setup
```bash
# Set up pipeline monitoring
nextflow run pipeline.nf -with-tower -profile production

# Custom monitoring hooks
nextflow run pipeline.nf \
  --hook-url https://monitoring.example.com/webhook \
  --notify-on-completion \
  --notify-on-failure
```

## Resources

### Documentation
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)
- [nf-core Developer Guide](https://nf-co.re/developers/)
- [Nextflow Tower Documentation](https://help.tower.nf/)

### Testing Tools
- [nf-test](https://nf-test.com/)
- [pytest-workflow](https://pytest-workflow.readthedocs.io/)
- [Nextflow CI/CD Examples](https://github.com/nextflow-io/nextflow/tree/master/.github/workflows)

### Performance Optimization
- [Nextflow Performance Tips](https://nextflow.io/docs/latest/tracing.html)
- [Resource Requirements Guide](https://nf-co.re/docs/usage/configuration)

## Looking Ahead

**Day 9 Preview**: Bring Your Own Data session including:
- Applying learned skills to participant datasets
- Troubleshooting real-world challenges
- Customizing pipelines for specific needs
- Preparing for independent analysis

---

**Key Learning Outcome**: Advanced Nextflow development enables creation of production-ready, scalable bioinformatics pipelines that can handle complex datasets across diverse computing environments while maintaining reproducibility and reliability.