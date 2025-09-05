# Day 4: Genome Assembly Essentials: QC, Identification, and assembly

**Date**: September 4, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Genome assembly and assessment
## Overview

Day 4 builds on Day 3 by starting with a recap on quality control, followed by genome assembly and assessment of assemblies to ensure.

## Learning Objectives

By the end of Day 4, you will be able to:

- Execute de novo genome assembly using SPAdes or other assemblers
- Assess assembly quality using QUAST and other metrics

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Recap: Quality checking and control, and species identification* | [Practical](https://github.com/Arash-Iranzadeh/Microbial-Genomics/blob/main/scripts/qc_cleaning_detection.sh) | Arash Iranzadeh |
| **10:30** | *Genome assembly, quality assessment* |[Notes](../day4/genome_assembly.md). | Ephifania Geza |
| **11:00** | **Break** | | |
| **11:30** | *Genome assembly, quality assessment: Continuation* | [Practical](../day4/practical_assembly.md) | Ephifania Geza |

## Key Topics

### 1. Genome assembly and assessing contigs
- De novo assembly algorithms and approaches
- Short-read vs long-read assembly strategies
- Assembly quality metrics and interpretation
- Assembly polishing and gap filling


## Datasets Used

### *V. cholerae* Outbreak Collection
- **Source**: Multi-country cholera outbreak (2019)
- **Samples**: 25 clinical isolates + environmental samples
- **Timespan**: 8-month outbreak period
- **Geographic**: Three countries, coastal regions
- **Epidemiological data**: Case demographics, travel history


## Tools Introduced

### Pangenome Analysis
- **SPAdes** - De novo genome assembler
- **Unicycler** - Hybrid assembly pipeline
- **Flye** - Long-read assembly
- **QUAST** - Assembly quality assessment

### Annotation Tools
- **Prokka** - Automated prokaryotic annotation
- **RAST** - Rapid Annotation using Subsystem Technology
- **NCBI PGAP** - Prokaryotic Genome Annotation Pipeline
- **Bakta** - Rapid bacterial genome annotation

- **Panaroo** - Pangenome pipeline
- **Roary** - Rapid large-scale prokaryote pangenome analysis
- **PPanGGOLiN** - Depicting microbial diversity via pangenomes

### Phylogenetic Analysis
- **IQ-TREE** - Maximum likelihood phylogenetic inference
- **RAxML** - Randomized Axelerated Maximum Likelihood
- **FastTree** - Approximately maximum-likelihood trees

### SNP Analysis
- **Snippy** - Rapid haploid variant calling
- **ParSNP** - Rapid core genome SNP typing
- **Gubbins** - Recombination detection in bacterial genomes

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


### Exercise 1: Pangenome Analysis (60 minutes)
Analyze the core and accessory genome of *V. cholerae* outbreak strains.

```bash
# Run Panaroo pangenome analysis
panaroo -i *.gff -o panaroo_output --clean-mode strict

# Visualize results
python3 scripts/visualize_pangenome.py panaroo_output/
```

### Exercise 2: Phylogenetic Tree Construction (60 minutes)
Build maximum likelihood trees from core genome alignments.

```bash
# Generate core genome alignment
snippy-core --ref reference.gbk snippy_output/*

# Build phylogenetic tree
iqtree -s core_alignment.aln -m GTR+G -bb 1000 -nt 4

# Visualize tree
figtree core_alignment.aln.treefile
```

### Exercise 3: Outbreak Investigation (45 minutes)
Integrate genomic and epidemiological data to investigate transmission.

```bash
# Calculate pairwise SNP distances
snp-dists core_alignment.aln > snp_distances.tsv

# Identify transmission clusters
cluster_analysis.py --snp-threshold 10 --epi-data metadata.csv
```

## Key Concepts

### Assembly Quality Metrics
| Metric | Good Assembly | Poor Assembly | Action |
|--------|---------------|---------------|--------|
| N50 | >50 kb | <10 kb | Optimize parameters |
| Contigs | <100 | >500 | Check contamination |
| Genome size | Expected ±10% | >20% difference | Review input data |
| Coverage | >50x | <20x | Sequence more |

### Pangenome Structure
- **Core genome**: Genes present in all strains (housekeeping functions)
- **Accessory genome**: Variable genes (adaptation, virulence, resistance)
- **Singleton genes**: Present in single strain only
- **Shell genes**: Present in several but not all strains

### Phylogenetic Inference
- **Substitution models**: GTR, HKY, JC69 for nucleotide evolution
- **Rate heterogeneity**: Gamma distribution for variable sites
- **Bootstrap support**: Statistical confidence in tree topology
- **Branch lengths**: Evolutionary distance (substitutions per site)

### SNP Thresholds for Outbreak Investigation
| Pathogen | SNP Threshold | Timeframe | Context |
|----------|---------------|-----------|---------|
| *M. tuberculosis* | 0-5 SNPs | Recent transmission | Same household |
| *S. Typhi* | 0-20 SNPs | Outbreak cluster | Weeks to months |
| *V. cholerae* | 0-10 SNPs | Epidemic spread | Days to weeks |
| *E. coli* O157 | 0-15 SNPs | Foodborne outbreak | Days |

## Advanced Topics

### Molecular Dating
- **Tip dating**: Using collection dates for molecular clock
- **Bayesian methods**: BEAST, MrBayes for time-resolved phylogenies
- **Substitution rates**: Pathogen-specific evolutionary rates

### Network Analysis
- **Minimum spanning trees**: Alternative to bifurcating phylogenies
- **Median-joining networks**: Visualization of reticulate evolution
- **Transmission networks**: Direct transmission inference

## Assessment Activities

### Individual Analysis
- Generate pangenome summary statistics
- Construct phylogenetic tree with bootstrap support
- Calculate SNP distances between outbreak isolates
- Identify transmission clusters based on genomic data

### Group Discussion
- Compare assembly strategies and results

## Common Challenges

- Interpret pangenome diversity in outbreak context
- Evaluate phylogenetic support for transmission hypotheses
- Discuss integration of genomic and epidemiological evidence
- Consider limitations of genomic outbreak investigation

## Common Challenges

### Low Phylogenetic Resolution
```bash
# Try different substitution models
iqtree -s alignment.aln -m MFP -bb 1000  # Model selection

# Use more informative sites only
iqtree -s alignment.aln -m GTR+G -bb 1000 --rate
```

### Recombination Issues
```bash
# Detect and remove recombinant regions
gubbins alignment.aln

# Use recombination-free alignment
iqtree -s alignment.filtered_polymorphic_sites.fasta
```

### Missing Epidemiological Links
```bash
# Lower SNP threshold for exploration
cluster_analysis.py --snp-threshold 20 --epi-data metadata.csv

# Consider longer transmission chains
transmission_chains.py --max-generations 3
```

## Resources

## Resources

### Assembly Resources
- [SPAdes Manual](http://cab.spbu.ru/software/spades/)
- [QUAST Documentation](http://quast.sourceforge.net/)
- [Assembly Best Practices](https://github.com/rrwick/Perfect-bacterial-genome-tutorial)

### Key Publications
- Page et al. (2015). Roary: rapid large-scale prokaryote pangenome analysis
- Tonkin-Hill et al. (2020). Producing polished prokaryotic pangenomes
- Croucher et al. (2015). Rapid phylogenetic analysis of bacterial genomes

### Software Documentation
- [IQ-TREE Manual](http://www.iqtree.org/doc/)
- [Panaroo Documentation](https://gtonkinhill.github.io/panaroo/)
- [Snippy Manual](https://github.com/tseemann/snippy)

### Online Resources
- [Microreact](https://microreact.org/) - Visualization platform
- [iTOL](https://itol.embl.de/) - Interactive tree visualization
- [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) - Tree viewing software

## Looking Ahead

**Day 5 Preview**: Metagenomics analysis including:
- Microbiome profiling from clinical samples
- Pathogen detection in complex communities
- Functional analysis of metagenomes
- Association with host health outcomes

## Homework (Optional)

1. Analyze pangenome diversity in additional pathogen datasets
2. Experiment with different phylogenetic methods and compare results
3. Read case studies of genomic outbreak investigations
4. Practice tree interpretation with published examples

---

**Key Takeaway**: Genomic analysis provides unprecedented resolution for understanding pathogen evolution and transmission, but interpretation requires careful integration with epidemiological context and understanding of method limitations.
