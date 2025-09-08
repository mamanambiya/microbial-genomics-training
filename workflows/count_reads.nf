#!/usr/bin/env nextflow

// Parameters you can change
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
        return [sample, fastq1, fastq2]
    }
    .set { samples_ch }

// Process to count reads in paired FASTQ files
process countReads {
    // Where to save results
    publishDir params.outdir, mode: 'copy'
    
    // Use sample name for process identification
    tag "$sample"

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output:
    path "${sample}.count"

    script:
    """
    echo "Counting reads in sample: ${sample}"
    echo "Forward reads (${fastq1}):"
    
    # Count reads in both files (compressed FASTQ)
    reads1=\$(zcat ${fastq1} | wc -l | awk '{print \$1/4}')
    reads2=\$(zcat ${fastq2} | wc -l | awk '{print \$1/4}')
    
    echo "Sample: ${sample}" > ${sample}.count
    echo "Forward reads: \$reads1" >> ${sample}.count
    echo "Reverse reads: \$reads2" >> ${sample}.count
    echo "Total read pairs: \$reads1" >> ${sample}.count
    
    echo "Finished counting ${sample}: \$reads1 read pairs"
    """
}

workflow {
    countReads(samples_ch)
}
