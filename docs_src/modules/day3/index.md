# Day 3: Genomic Characterization

**Date**: September 3, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Quality control, genome assembly, assessment 

## Overview

Day 3 introduces participants to the Ilifu high-performance computing infrastructure and the essential data pre-processing steps. These include quality control, species identification, followed by genome assembly and quality assessment. These foundational skills are critical for all downstream genomic analyses.

## Learning Objectives

By the end of Day 3, you will be able to:

- Submit and manage jobs using the SLURM scheduler
- Understand resource allocation and job queue management on HPC systems
- Perform quality checking and control on sequencing data using FastQC
- Identify species from genomic data using Kraken2 and other tools
- Execute de novo genome assembly using SPAdes or other assemblers
- Assess assembly quality using QUAST and other metrics

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Introduction High Performance Computing (HPC) – Ilifu* | [Notes](../../day2/hpc-ilifu-training.md) • [Practical 1](../../day2/slurm-practical-tutorial.md) • [Practical 2](../../day2/unix-commands-pathogen-examples.md) | Mamana Mbiyavanga |
| **10:00** | *Quality checking and control, as well as species identification* | | Arash Iranzadeh |
| **11:30** | **Break** | | |
| **10:00** | *Genome assembly, quality assessment* | [Notes](../../day3/genome_assembly_notes.md)  • [Practical](../../day3/practical_genome_assmbly.md) | Ephifania Geza |

## Key Topics

### 1. Genome Assembly and Quality Assessment
- Quality checking, and adapter and low read quality filtering
- Contamination detection and removal
- Species identification or confirmation
- De novo assembly algorithms and approaches
- Short-read vs long-read assembly strategies
- Assembly quality metrics and interpretation
- Assembly polishing and gap filling


## Tools and Software

### Assembly Tools
- **FASTQC** - Quality checking
- **MULTIQC** - Quality checking and amalgation of reports
- **SPAdes** - De novo genome assembler
- **Unicycler** - Hybrid assembly pipeline
- **Flye** - Long-read assembly
- **QUAST** - Assembly quality assessment

### Annotation Tools
- **Prokka** - Automated prokaryotic annotation
- **RAST** - Rapid Annotation using Subsystem Technology
- **NCBI PGAP** - Prokaryotic Genome Annotation Pipeline
- **Bakta** - Rapid bacterial genome annotation

## Hands-on Exercises

### Exercise 1: Genome Assembly and Quality Assessment (60 minutes)
Assemble bacterial genomes and evaluate assembly quality.

```bash
# De novo assembly with SPAdes
spades.py --careful -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o assembly_output/

# Assembly quality assessment
quast.py assembly_output/scaffolds.fasta -o quast_results/

# Check assembly statistics
assembly-stats assembly_output/scaffolds.fasta

# Contamination screening
kraken2 --db minikraken2_v2 assembly_output/scaffolds.fasta --report contamination_check.txt
```



## Key Concepts

### Assembly Quality Metrics
| Metric | Good Assembly | Poor Assembly | Action |
|--------|---------------|---------------|--------|
| N50 | >50 kb | <10 kb | Optimize parameters |
| Contigs | <100 | >500 | Check contamination |
| Genome size | Expected ±10% | >20% difference | Review input data |
| Coverage | >50x | <20x | Sequence more |


## Assessment Activities

### Individual Analysis
- Complete genome assembly workflow
- Perform quality assessment and interpretation

### Group Discussion
- Compare assembly strategies and results

## Common Challenges

### Assembly Issues
```bash
# Low coverage assemblies
spades.py --careful --cov-cutoff 5 -1 R1.fastq -2 R2.fastq -o low_cov_assembly/

# Contamination removal
# Remove contaminant contigs based on taxonomy
seqtk subseq scaffolds.fasta clean_contigs.txt > clean_assembly.fasta
```


## Clinical Applications

### Routine Surveillance
- Mantaining data quality control standards
- Rapid species identification and typing

## Resources

### Assembly Resources
- [SPAdes Manual](http://cab.spbu.ru/software/spades/)
- [QUAST Documentation](http://quast.sourceforge.net/)
- [Assembly Best Practices](https://github.com/rrwick/Perfect-bacterial-genome-tutorial)

## Looking Ahead

**Day 4 Preview**: Comparative Genomics including:
- Genome annotation and intepretation
- Antimicrobial resistance gene prediction
- Sequence typing

---

**Key Learning Outcome**: Genome assembly, MLST, serotyping, AMR detection, and mobile genetic element analysis provide comprehensive characterization capabilities essential for clinical genomics and public health surveillance.
