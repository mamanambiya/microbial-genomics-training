# Day 9: Bring Your Own Data

**Date**: September 11, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Independent analysis of participant datasets

## Overview

Day 9 is dedicated to hands-on application of all techniques learned throughout the course. Participants will analyze their own datasets with guidance from trainers, troubleshoot real-world challenges, and develop customized analysis approaches for their specific research questions.

## Learning Objectives

By the end of Day 9, you will be able to:

- Apply learned bioinformatics techniques to your own research data
- Troubleshoot common issues encountered in real-world analyses
- Adapt standard protocols to meet specific research requirements
- Develop analysis strategies for novel research questions
- Integrate multiple analysis approaches into comprehensive workflows
- Plan sustainable bioinformatics practices for ongoing research

## Schedule

| Time (CAT) | Topic | Trainer |
|------------|-------|------------|
| **09:00** | *Participants to analyse their own data* | All trainers |
| **11:30** | **Break** | |
| **12:00** | *Participants to analyse their own data* | All trainers |

## Session Structure

### Opening Session (09:00-09:15)
- Brief recap of key course concepts
- Overview of available support and resources
- Formation of analysis groups based on data types
- Technical setup verification

### Individual Analysis Time (09:15-11:30)
- Independent work on participant datasets
- One-on-one consultation with trainers
- Peer collaboration and knowledge sharing
- Documentation of analysis approaches

### Progress Sharing (11:45-12:00)
- Brief presentations of initial findings
- Discussion of challenges encountered
- Sharing of successful analysis strategies

### Advanced Analysis (12:00-13:00)
- Continued independent work
- Focus on complex analyses or troubleshooting
- Preparation for Day 10 presentations
- Final consultations with trainers

## Data Types and Analysis Approaches

### Genomic Data
**Common analyses for participants with genomic datasets:**

- **Quality control and preprocessing**
  - FastQC assessment and interpretation
  - Adapter trimming and quality filtering
  - Species identification and contamination detection

- **Genome assembly and annotation**
  - De novo assembly optimization
  - Assembly quality assessment
  - Functional annotation and gene prediction

- **Comparative genomics**
  - MLST and serotyping
  - Antimicrobial resistance gene detection
  - Phylogenetic analysis and clustering

### Metagenomic Data
**For participants with microbiome or metagenomic samples:**

- **Community profiling**
  - Taxonomic classification
  - Abundance estimation and normalization
  - Diversity analysis (alpha and beta)

- **Functional analysis**
  - Pathway reconstruction
  - Antimicrobial resistance profiling
  - Metabolic potential assessment

- **Clinical applications**
  - Pathogen detection in complex samples
  - Co-infection analysis
  - Treatment response monitoring

### Outbreak Investigation Data
**For epidemiological and outbreak datasets:**

- **Transmission analysis**
  - SNP-based clustering
  - Phylogenetic reconstruction
  - Temporal and geographic analysis

- **Resistance surveillance**
  - Multi-drug resistance patterns
  - Resistance gene distribution
  - Treatment outcome correlations

## Technical Support Available

### Computational Resources
- Access to high-performance computing cluster
- Pre-installed bioinformatics software environments
- Container images for reproducible analysis
- Shared storage for large datasets

### Analysis Pipelines
- Nextflow workflows developed during the course
- Customizable analysis templates
- Pre-configured environment profiles
- Automated reporting tools

### Expert Guidance
**Trainer specializations available:**

| Trainer | Expertise Areas |
|---------|----------------|
| **Ephifania Geza** | Genomic surveillance, AMR analysis, metagenomics, clinical applications |
| **Arash Iranzadeh** | Phylogenomics, comparative genomics, outbreak investigation |
| **Sindiswa Lukhele** | Sequencing technologies, quality control, species identification |
| **Mamana Mbiyavanga** | Workflow development, HPC systems, pipeline optimization |

## Common Analysis Workflows

### Genomic Surveillance Workflow
```bash
# 1. Initial data assessment
fastqc raw_data/*.fastq.gz
multiqc fastqc_results/

# 2. Species identification
kraken2 --db minikraken2_v2 --paired sample_R1.fastq sample_R2.fastq

# 3. Quality trimming
trimmomatic PE sample_R1.fastq sample_R2.fastq \
    sample_R1_trimmed.fastq sample_R1_unpaired.fastq \
    sample_R2_trimmed.fastq sample_R2_unpaired.fastq \
    SLIDINGWINDOW:4:20 MINLEN:50

# 4. Assembly
spades.py -1 sample_R1_trimmed.fastq -2 sample_R2_trimmed.fastq -o assembly/

# 5. Assembly quality assessment
quast.py assembly/scaffolds.fasta -o quast_results/

# 6. Annotation
prokka assembly/scaffolds.fasta --outdir annotation/ --prefix sample
```

### Metagenomic Analysis Workflow
```bash
# 1. Host DNA removal (if applicable)
kneaddata --input sample_R1.fastq --input sample_R2.fastq \
    --reference-db human_genome --output cleaned_data/

# 2. Taxonomic profiling
metaphlan cleaned_data/sample_paired_1.fastq,cleaned_data/sample_paired_2.fastq \
    --bowtie2out sample.bowtie2.bz2 --nproc 4 --input_type fastq \
    --output_file sample_profile.txt

# 3. Functional profiling
humann --input cleaned_data/sample_paired.fastq \
    --output functional_analysis/ --nucleotide-database chocophlan \
    --protein-database uniref90

# 4. Diversity analysis in R
Rscript diversity_analysis.R sample_profile.txt metadata.csv
```

## Troubleshooting Guide

### Common Issues and Solutions

#### Low-Quality Data
```bash
# Check read quality distribution
fastqc *.fastq.gz

# Aggressive quality trimming if needed
trimmomatic PE input_R1.fastq input_R2.fastq \
    output_R1.fastq unpaired_R1.fastq \
    output_R2.fastq unpaired_R2.fastq \
    SLIDINGWINDOW:4:25 LEADING:20 TRAILING:20 MINLEN:75

# Consider different assembly strategies
# For poor quality data, try more conservative parameters
spades.py --careful --cov-cutoff auto -1 R1.fastq -2 R2.fastq -o assembly/
```

#### Contamination Issues
```bash
# Check for contamination
kraken2 --db standard --paired sample_R1.fastq sample_R2.fastq \
    --report contamination_report.txt

# Remove contaminant sequences
extract_kraken_reads.py -k sample.kraken -s1 sample_R1.fastq \
    -s2 sample_R2.fastq -o clean_R1.fastq -o2 clean_R2.fastq \
    --exclude --taxid 9606  # Exclude human reads
```

#### Assembly Problems
```bash
# If SPAdes fails, try different assemblers
# Unicycler for hybrid assembly
unicycler -1 short_R1.fastq -2 short_R2.fastq -l long_reads.fastq -o assembly/

# Or SKESA for quick assembly
skesa --reads short_R1.fastq,short_R2.fastq --cores 8 > assembly.fasta
```

#### Memory/Resource Issues
```bash
# Monitor resource usage
htop

# Reduce memory usage for large datasets
spades.py --memory 16 -1 R1.fastq -2 R2.fastq -o assembly/

# Use subsampling for initial testing
seqtk sample -s100 input_R1.fastq 100000 > subset_R1.fastq
seqtk sample -s100 input_R2.fastq 100000 > subset_R2.fastq
```

## Analysis Documentation

### Laboratory Notebook Template
```markdown
# Analysis Log: [Your Dataset Name]
**Date**: [Current Date]
**Analyst**: [Your Name]
**Data Source**: [Description of samples]

## Objectives
- Primary research question
- Specific analyses planned
- Expected outcomes

## Data Description
- Sample type and collection method
- Sequencing platform and parameters
- Data quality metrics

## Analysis Steps
### Step 1: Quality Control
- Command used: `fastqc *.fastq.gz`
- Results: [Summary of quality metrics]
- Decision: [Any quality filtering applied]

### Step 2: [Next Analysis]
- Command: [Exact command used]
- Parameters chosen: [Rationale for parameter selection]
- Results: [Key findings]

## Challenges Encountered
- Issue: [Description of problem]
- Solution attempted: [What was tried]
- Outcome: [Whether resolved]

## Key Findings
- [Major results from analysis]
- [Statistical summaries]
- [Biological interpretations]

## Next Steps
- Additional analyses needed
- Questions raised
- Follow-up experiments
```

## Resource Management

### Data Organization
```bash
# Recommended directory structure
project_name/
├── raw_data/          # Original sequencing files
├── quality_control/   # QC reports and cleaned data
├── analysis/          # Main analysis outputs
├── scripts/           # Custom scripts and commands
├── results/           # Final results and figures
└── documentation/     # Analysis logs and notes
```

### Backup Strategies
```bash
# Regular backup of important results
rsync -av results/ backup_drive/project_results/
tar -czf analysis_$(date +%Y%m%d).tar.gz analysis/

# Version control for scripts
git init
git add scripts/
git commit -m "Initial analysis scripts"
```

## Collaboration Guidelines

### Peer Support
- Form analysis groups based on similar data types
- Share successful parameter combinations
- Collaborate on troubleshooting challenging datasets
- Review each other's analysis approaches

### Trainer Consultation
- Prepare specific questions about your data
- Document issues with exact error messages
- Have your analysis objectives clearly defined
- Be ready to explain your research context

## Assessment and Preparation for Day 10

### Presentation Preparation
Participants should prepare a 5-minute presentation covering:

1. **Research Question**: What you aimed to investigate
2. **Data Overview**: Type and source of your dataset
3. **Methods Applied**: Which course techniques you used
4. **Key Results**: Main findings from your analysis
5. **Challenges**: Obstacles encountered and solutions found
6. **Future Directions**: Next steps for your research

### Technical Documentation
- Save all commands used in a script file
- Document parameter choices and rationale
- Prepare summary statistics and key figures
- Note any analysis limitations or assumptions

## Success Metrics

By the end of Day 9, participants should have:

- [ ] Successfully processed their own dataset
- [ ] Applied at least 3 different analysis techniques from the course
- [ ] Documented their analysis workflow
- [ ] Identified key findings relevant to their research
- [ ] Prepared materials for Day 10 presentation
- [ ] Established ongoing analysis plan

## Resources for Continued Learning

### Online Communities
- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [BioStars Forum](https://www.biostars.org/)
- [Galaxy Community](https://galaxyproject.org/community/)

### Software Documentation
- Tool-specific manuals and tutorials
- GitHub repositories for pipeline development
- Container registries for reproducible environments

### Professional Development
- Local bioinformatics user groups
- International conferences and workshops
- Online course platforms for advanced topics

## Looking Ahead

**Day 10 Preview**: Wrap-up session including:
- Participant presentations of analysis results
- Discussion of lessons learned and best practices
- Information about ongoing support resources
- Course completion and next steps planning

---

**Key Learning Outcome**: Independent application of bioinformatics skills to real research data builds confidence and reveals the practical challenges and rewards of computational biology in actual research contexts.