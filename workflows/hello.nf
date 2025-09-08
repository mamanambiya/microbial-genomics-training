#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Parameters
params.samples = ['sample1', 'sample2', 'sample3']

// Define a process (a step in your pipeline)
process sayHello {
    // What this process does
    input:
    val sample_name

    // What it produces
    output:
    stdout

    // The actual command
    script:
    """
    echo "Hello from ${sample_name}!"
    """
}

// Main workflow (DSL2 style)
workflow {
    // Create a channel (think of it as a conveyor belt for data)
    samples_ch = Channel.from(params.samples)

    // Run the process
    sayHello(samples_ch)

    // Show the results
    sayHello.out.view()
}
