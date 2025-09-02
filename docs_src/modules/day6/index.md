# Day 6: Nextflow Pipeline Development & GitHub

**Date**: September 8, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Pipeline completion, version control with GitHub, and collaborative development

## Overview

Day 6 continues Nextflow pipeline development with a focus on completing the genomic analysis pipeline, implementing version control using GitHub, and establishing collaborative development practices. Participants will learn to manage their pipelines using Git, create repositories, handle pull requests, and implement continuous integration for their bioinformatics workflows.

## Learning Objectives

By the end of Day 6, you will be able to:

- Complete a functional Nextflow pipeline for genomic analysis
- Use Git for version control of bioinformatics pipelines
- Create and manage GitHub repositories for pipeline projects
- Implement collaborative development workflows using branches and pull requests
- Set up continuous integration for pipeline testing
- Document pipelines effectively using README files and wikis
- Share and publish pipelines for community use

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Pipeline development continued and GitHub integration* | | Mamana Mbiyavanga |
| **11:30** | **Break** | | |
| **12:00** | *Version control, collaboration, and pipeline sharing* | | Mamana Mbiyavanga |

## Key Topics

### 1. Completing the Nextflow Pipeline
- Finalizing QC, assembly, and annotation processes
- Adding error handling and recovery mechanisms
- Implementing pipeline configuration profiles
- Testing with real genomic data

### 2. Version Control with Git
- Git fundamentals for bioinformatics projects
- Tracking pipeline changes and versions
- Managing configuration files and parameters
- Best practices for commit messages

### 3. GitHub for Pipeline Development
- Creating and organizing GitHub repositories
- Writing effective README documentation
- Using GitHub Issues for project management
- Releases and semantic versioning

### 4. Collaborative Development
- Branching strategies for pipeline development
- Creating and reviewing pull requests
- Resolving merge conflicts
- Team workflows and code review

### 5. Continuous Integration
- Setting up GitHub Actions for pipeline testing
- Automated testing strategies
- Container building and registry management
- Documentation generation and deployment

## Tools and Software

### Pipeline Development
- **Nextflow** - Workflow orchestration
- **nf-core tools** - Pipeline templates and tools
- **Docker/Singularity** - Container platforms

### Version Control
- **Git** - Version control system
- **GitHub** - Code hosting and collaboration
- **GitHub CLI** - Command-line interface for GitHub
- **GitHub Desktop** - GUI for Git operations

### Continuous Integration
- **GitHub Actions** - CI/CD platform
- **nf-test** - Nextflow testing framework
- **pytest-workflow** - Pipeline testing tool
- **DESMAN** - De novo extraction of strains from metagenomes

## Hands-on Exercises

### Exercise 1: Rare Pathogen Detection (60 minutes)
Use metagenomic data to identify rare pathogens in clinical samples.

```bash
# Run Kraken2 with comprehensive database
kraken2 --db kraken2_plus --paired sample_R1.fastq sample_R2.fastq \
    --output sample.kraken --report sample_report.txt

# Filter for potential pathogens
python3 scripts/filter_pathogens.py sample_report.txt \
    --pathogen-list pathogen_database.txt \
    --min-reads 10 --output rare_pathogens.txt

# Validate findings with targeted BLAST
blastn -query extracted_sequences.fasta -db nt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 5 -out validation_blast.txt
```

### Exercise 2: Co-infection Analysis (75 minutes)
Analyze samples with multiple pathogen species.

```bash
# Quantify pathogen abundances
metaphlan sample_metagenome.fastq.gz --input_type fastq \
    --bowtie2out sample_bowtie2.bz2 --nproc 4 \
    --output_file sample_profile.txt

# Identify co-infections
python3 scripts/coinfection_analysis.py sample_profile.txt \
    --abundance-threshold 0.1 --output coinfection_report.txt

# Strain-level analysis for dominant pathogens
strainest sample_bowtie2.bz2 reference_genomes/ \
    --output strain_analysis/
```

### Exercise 3: Variant Calling (90 minutes)
Call variants and correlate with phenotypic data.

```bash
# Map reads to reference genome
bwa mem reference.fasta sample_R1.fastq sample_R2.fastq | \
    samtools sort -o sample_sorted.bam
samtools index sample_sorted.bam

# Call variants with GATK
gatk HaplotypeCaller -R reference.fasta -I sample_sorted.bam \
    -O sample_variants.vcf

# Annotate variants
snpEff ann reference_genome sample_variants.vcf > sample_annotated.vcf

# Correlate with phenotype data
python3 scripts/genotype_phenotype_correlation.py \
    sample_annotated.vcf phenotype_data.csv \
    --output correlation_analysis.txt
```

## Key Concepts

### Detection Sensitivity
- **Limit of detection**: Minimum number of reads for reliable identification
- **False positive control**: Contamination and database artifacts
- **Validation requirements**: Confirmatory testing approaches
- **Clinical significance**: Distinguishing colonization from infection

### Co-infection Dynamics
- **Competitive exclusion**: Pathogens competing for resources
- **Synergistic interactions**: Enhanced pathogenicity in co-infections
- **Treatment complications**: Multi-drug considerations
- **Epidemiological implications**: Transmission patterns

### Variant Classification
| Variant Type | Impact | Example |
|--------------|--------|---------|
| Synonymous | Silent (no amino acid change) | c.123C>T (p.Leu41=) |
| Missense | Amino acid change | c.123G>A (p.Gly41Ser) |
| Nonsense | Premature stop codon | c.123C>T (p.Gln41*) |
| Frameshift | Reading frame alteration | c.123delA (p.Leu41fs) |

## Clinical Applications

### Diagnostic Applications
- **Syndromic testing**: Multiple pathogen panels
- **Outbreak investigation**: Source identification
- **Treatment monitoring**: Resistance development
- **Infection control**: Transmission prevention

### Antimicrobial Resistance
- **Rapid resistance testing**: Genomic predictions
- **Resistance mechanism identification**: Gene and mutation analysis
- **Treatment guidance**: Targeted therapy selection
- **Surveillance**: Population-level resistance monitoring

## Assessment Activities

### Individual Analysis
- Identify rare pathogens in provided metagenomic dataset
- Quantify co-infection burden in clinical samples
- Perform variant calling and annotation on genomic data
- Correlate genomic findings with resistance phenotypes

### Group Discussion
- Compare sensitivity of different pathogen detection methods
- Interpret co-infection patterns in clinical context
- Evaluate clinical significance of identified variants
- Discuss challenges in genotype-phenotype correlation

## Common Challenges

### False Positive Control
```bash
# Include negative controls in analysis
kraken2 --db comprehensive_db --paired \
    negative_control_R1.fastq negative_control_R2.fastq \
    --output neg_control.kraken --report neg_control_report.txt

# Remove contaminants found in negative controls
python3 scripts/remove_contaminants.py sample_report.txt \
    neg_control_report.txt --output cleaned_sample_report.txt
```

### Low Coverage Variants
```bash
# Adjust variant calling parameters for low coverage
freebayes -f reference.fasta --min-coverage 5 --min-alternate-count 2 \
    --min-alternate-fraction 0.1 sample_sorted.bam > low_cov_variants.vcf

# Validate with manual inspection
samtools tview sample_sorted.bam reference.fasta -p chr1:12345
```

### Complex Structural Variants
```bash
# Use specialized tools for structural variants
delly call -g reference.fasta sample_sorted.bam -o structural_variants.bcf

# Convert and filter
bcftools view structural_variants.bcf | \
    bcftools filter -i 'QUAL>20' > filtered_sv.vcf
```

## Resources

### Essential Reading
- Simner et al. (2018). Understanding the promises and hurdles of metagenomic next-generation sequencing
- Chiu & Miller (2019). Clinical metagenomics
- Ellington et al. (2017). The role of whole genome sequencing in antimicrobial susceptibility testing

### Software Documentation
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)
- [PathSeq Documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037224932)
- [SnpEff Manual](https://pcingola.github.io/SnpEff/)

### Databases
- [CARD](https://card.mcmaster.ca/) - Comprehensive Antibiotic Resistance Database
- [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) - Resistance gene identification
- [NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) - Pathogen surveillance

## Looking Ahead

**Day 7 Preview**: Nextflow workflows including:
- Introduction to reproducible workflow systems
- Nextflow basics and nf-core community
- Building automated analysis pipelines
- Pipeline development best practices

## Homework (Optional)

1. Analyze co-infection patterns in additional clinical datasets
2. Practice variant calling with different parameter settings
3. Explore resistance gene databases for target pathogens
4. Read case studies of metagenomic pathogen detection

---

**Key Learning Outcome**: Advanced metagenomic analysis enables detection of complex infection scenarios and genomic variants, providing critical information for clinical decision-making and public health surveillance.