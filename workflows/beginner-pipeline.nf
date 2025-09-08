#!/usr/bin/env nextflow

/*
 * Beginner's Microbial Genomics Pipeline
 * =====================================
 * 
 * This pipeline demonstrates basic genomics analysis steps:
 * 1. Quality control with FastQC
 * 2. Read trimming with Trimmomatic  
 * 3. Genome assembly with SPAdes
 * 4. Assembly quality assessment with QUAST
 * 
 * Perfect for learning Nextflow basics!
 */

// Pipeline parameters - you can change these
params.reads = "data/*_{R1,R2}.fastq.gz"
params.outdir = "results"
params.help = false

// Show help message
if (params.help) {
    log.info """
    Beginner's Microbial Genomics Pipeline
    =====================================
    
    Usage:
    nextflow run beginner-pipeline.nf --reads 'data/*_{R1,R2}.fastq.gz'
    
    Parameters:
    --reads     Path to paired-end FASTQ files (default: data/*_{R1,R2}.fastq.gz)
    --outdir    Output directory (default: results)
    --help      Show this help message
    
    Example:
    nextflow run beginner-pipeline.nf --reads 'samples/*_{R1,R2}.fastq.gz' --outdir my_results
    """
    exit 0
}

// Print pipeline info
log.info """
Starting Beginner's Microbial Genomics Pipeline
==============================================
Reads: ${params.reads}
Output: ${params.outdir}
"""

// Create input channel
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "No read files found matching: ${params.reads}" }
    .set { read_pairs_ch }

/*
 * STEP 1: Quality Control with FastQC
 * Checks the quality of your sequencing data
 */
process fastqc {
    tag "$sample_id"
    container 'biocontainers/fastqc:v0.11.9'
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{zip,html}"
    tuple val(sample_id), path(reads) into trimming_ch
    
    script:
    """
    echo "Running FastQC on sample: ${sample_id}"
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

/*
 * STEP 2: Read Trimming with Trimmomatic
 * Removes low-quality parts of sequences
 */
process trimming {
    tag "$sample_id"
    container 'staphb/trimmomatic:latest'
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    cpus 2
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz") into assembly_ch
    path "${sample_id}_trimming.log"
    
    script:
    """
    echo "Trimming reads for sample: ${sample_id}"
    trimmomatic PE -threads ${task.cpus} \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_trimmed.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \\
        ${sample_id}_R2_trimmed.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \\
        SLIDINGWINDOW:4:20 MINLEN:50 \\
        2> ${sample_id}_trimming.log
    """
}

/*
 * STEP 3: Genome Assembly with SPAdes
 * Puts the sequencing reads back together to reconstruct the genome
 */
process assembly {
    tag "$sample_id"
    container 'staphb/spades:latest'
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    cpus 4
    memory '8 GB'
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta") into quast_ch
    path "${sample_id}_assembly"
    
    script:
    """
    echo "Assembling genome for sample: ${sample_id}"
    spades.py \\
        -1 ${trimmed_reads[0]} \\
        -2 ${trimmed_reads[1]} \\
        -o ${sample_id}_assembly \\
        --threads ${task.cpus} \\
        --memory \$(echo ${task.memory} | sed 's/ GB//')
    
    # Copy the final contigs with a clear name
    cp ${sample_id}_assembly/contigs.fasta ${sample_id}_contigs.fasta
    """
}

/*
 * STEP 4: Assembly Quality Assessment with QUAST
 * Checks how good your genome assembly is
 */
process quast {
    tag "$sample_id"
    container 'staphb/quast:latest'
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    
    input:
    tuple val(sample_id), path(contigs)
    
    output:
    path "${sample_id}_quast"
    
    script:
    """
    echo "Assessing assembly quality for sample: ${sample_id}"
    quast.py \\
        --output-dir ${sample_id}_quast \\
        --threads ${task.cpus} \\
        ${contigs}
    """
}

/*
 * Pipeline completion message
 */
workflow.onComplete {
    log.info """
    Pipeline completed!
    =================
    Results are in: ${params.outdir}
    
    Check these folders:
    - 01_fastqc: Quality control reports
    - 02_trimmed: Cleaned sequencing reads  
    - 03_assembly: Assembled genomes
    - 04_quast: Assembly quality reports
    
    Success: ${workflow.success}
    Duration: ${workflow.duration}
    """
}
