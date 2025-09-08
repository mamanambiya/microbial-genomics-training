#!/usr/bin/env nextflow

process sayHello {
    output:
    stdout

    script:
    """
    echo 'Hello, Nextflow!'
    """
}

workflow {
    sayHello()
}
