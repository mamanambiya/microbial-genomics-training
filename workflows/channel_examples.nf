#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

workflow {
    // Channel from file pairs
    reads_ch = Channel.fromFilePairs("/data/Dataset_Mt_Vc/tb/raw_data/*_{1,2}.fastq.gz")
    reads_ch.view { sample, files -> "Sample: $sample, Files: $files" }

    // Channel from list
    samples_ch = Channel.from(['sample1', 'sample2', 'sample3'])
    samples_ch.view { "Processing: $it" }

    // Channel from path pattern
    ref_ch = Channel.fromPath("*.fasta")
    ref_ch.view { "Reference file: $it" }

    // Channel from CSV
    csv_ch = Channel
        .fromPath("samplesheet.csv")
        .splitCsv(header: true)
        .map { row -> [row.sample, row.fastq_1, row.fastq_2] }
    csv_ch.view { sample, r1, r2 -> "Sample: $sample, R1: $r1, R2: $r2" }

    // Channel operations
    numbers_ch = Channel.from(1, 2, 3, 4, 5)
    
    // Filter even numbers
    even_ch = numbers_ch.filter { it % 2 == 0 }
    even_ch.view { "Even number: $it" }
    
    // Map operation
    squared_ch = Channel.from(1, 2, 3, 4, 5).map { it * it }
    squared_ch.view { "Squared: $it" }
    
    // Combine channels
    letters_ch = Channel.from('A', 'B', 'C')
    numbers2_ch = Channel.from(1, 2, 3)
    combined_ch = numbers2_ch.combine(letters_ch)
    combined_ch.view { num, letter -> "Combined: $num-$letter" }
}
