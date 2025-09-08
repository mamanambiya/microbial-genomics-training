#!/usr/bin/env nextflow

// Parameters
params.input = "samplesheet.csv"
params.outdir = "results"

// Read sample sheet and create channel
Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row ->
        def sample = row.sample
        def fastq1 = file(row.fastq_1)
        def fastq2 = file(row.fastq_2)
        return [sample, [fastq1, fastq2]]
    }
    .set { read_pairs_ch }

// No need to duplicate in DSL2 - we can use the same channel multiple times

// FastQC on raw reads
process fastqc_raw {
    module 'fastqc/0.12.1'
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    echo "Running FastQC on raw reads: ${sample_id}"
    fastqc ${reads}
    """
}

// Trimmomatic for quality trimming
process trimmomatic {
    module 'trimmomatic/0.39'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_*_paired.fastq.gz")
    path "${sample_id}_*_unpaired.fastq.gz"

    script:
    """
    echo "Running Trimmomatic on ${sample_id}"
    
    trimmomatic PE -threads 2 \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_paired.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \\
        ${sample_id}_R2_paired.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \\
        LEADING:3 TRAILING:3 \\
        SLIDINGWINDOW:4:15 \\
        MINLEN:36
    
    echo "Trimming completed for ${sample_id}"
    """
}

// FastQC on trimmed reads
process fastqc_trimmed {
    module 'fastqc/0.12.1'
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    echo "Running FastQC on trimmed reads: ${sample_id}"
    fastqc ${reads}
    """
}

// SPAdes genome assembly
process spades_assembly {
    module 'spades/4.2.0'
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    tag "$sample_id"

    // Assembly can be memory intensive - good for cluster submission
    memory '8 GB'
    cpus 4
    time '2h'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_assembly")
    path "${sample_id}_assembly/contigs.fasta"

    script:
    """
    echo "Running SPAdes assembly for ${sample_id}"
    echo "Input files: ${reads}"

    # Run SPAdes with paired-end reads
    spades.py \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --threads ${task.cpus} \\
        --memory \$(echo ${task.memory} | sed 's/ GB//') \\
        --only-assembler \\
        -o ${sample_id}_assembly

    echo "Assembly completed for ${sample_id}"
    echo "Contigs written to: ${sample_id}_assembly/contigs.fasta"

    # Show basic assembly stats
    if [ -f "${sample_id}_assembly/contigs.fasta" ]; then
        echo "=== Assembly Statistics for ${sample_id} ==="
        echo "Number of contigs: \$(grep -c '>' ${sample_id}_assembly/contigs.fasta)"
        echo "Total assembly size: \$(grep -v '>' ${sample_id}_assembly/contigs.fasta | wc -c) bp"
    fi
    """
}

// Prokka genome annotation
process prokka_annotation {
    module 'prokka/1.14.6'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    tag "$sample_id"

    // Annotation requires moderate resources
    memory '4 GB'
    cpus 2
    time '1h'

    input:
    tuple val(sample_id), path(assembly_dir)
    path contigs_file

    output:
    tuple val(sample_id), path("${sample_id}_annotation")
    path "${sample_id}_annotation/${sample_id}.gff"

    script:
    """
    echo "Running Prokka annotation for ${sample_id}"
    echo "Input assembly: ${contigs_file}"

    # Run Prokka annotation for Mycobacterium tuberculosis
    prokka --outdir ${sample_id}_annotation \\
           --prefix ${sample_id} \\
           --cpus ${task.cpus} \\
           --genus Mycobacterium \\
           --species tuberculosis \\
           --kingdom Bacteria \\
           ${contigs_file}

    echo "Annotation completed for ${sample_id}"
    echo "Results written to: ${sample_id}_annotation/"

    # Show basic annotation stats
    if [ -f "${sample_id}_annotation/${sample_id}.gff" ]; then
        echo "=== Annotation Statistics for ${sample_id} ==="
        echo "Total features: \$(grep -v '^#' ${sample_id}_annotation/${sample_id}.gff | wc -l)"
        echo "CDS features: \$(grep -v '^#' ${sample_id}_annotation/${sample_id}.gff | grep 'CDS' | wc -l)"
        echo "Gene features: \$(grep -v '^#' ${sample_id}_annotation/${sample_id}.gff | grep 'gene' | wc -l)"
    fi
    """
}

workflow {
    // Run FastQC on raw reads
    fastqc_raw_results = fastqc_raw(read_pairs_ch)

    // Run Trimmomatic
    trimmomatic_results = trimmomatic(read_pairs_ch)

    // Run FastQC on trimmed reads
    fastqc_trimmed_results = fastqc_trimmed(trimmomatic_results[0])

    // Run SPAdes assembly on trimmed reads
    assembly_results = spades_assembly(trimmomatic_results[0])

    // Run Prokka annotation on assembled contigs
    annotation_results = prokka_annotation(assembly_results[0], assembly_results[1])

    // Show results
    fastqc_raw_results.view { "Raw FastQC: $it" }
    fastqc_trimmed_results.view { "Trimmed FastQC: $it" }
    assembly_results[0].view { "Assembly completed: $it" }
    assembly_results[1].view { "Contigs file: $it" }
    annotation_results[0].view { "Annotation completed: $it" }
    annotation_results[1].view { "GFF file: $it" }
}
