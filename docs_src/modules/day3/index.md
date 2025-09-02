# Day 3: Genomic Characterization

**Date**: September 3, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Quality control, genome assembly, assessment and annotation

## Overview

Day 3 introduces essential genomic characterization techniques starting with quality control and species identification, followed by genome assembly and quality assessment, and concluding with genome annotation. These foundational skills are critical for all downstream genomic analyses.

## Learning Objectives

By the end of Day 3, you will be able to:

- Perform quality checking and control on sequencing data using FastQC
- Identify species from genomic data using Kraken2 and other tools
- Execute de novo genome assembly using SPAdes or other assemblers
- Assess assembly quality using QUAST and other metrics
- Annotate bacterial genomes using Prokka
- Interpret and utilize genome annotation outputs

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Quality checking and control, as well as species identification* | | Arash Iranzadeh |
| **10:00** | *Genome assembly, quality assessment* | | Ephifania Geza |
| **11:30** | **Break** | | |
| **12:00** | *Genome Annotation* | | Arash Iranzadeh |

## Key Topics

### 1. Genome Assembly and Quality Assessment
- De novo assembly algorithms and approaches
- Short-read vs long-read assembly strategies
- Assembly quality metrics and interpretation
- Contamination detection and removal
- Assembly polishing and gap filling

### 2. Genome Annotation
- Gene prediction and functional annotation
- Prokka automated annotation pipeline
- Manual curation and quality control
- Annotation databases and resources
- Comparative annotation approaches

### 3. Multi-locus Sequence Typing (MLST)
- MLST principles and methodology
- Housekeeping gene selection and analysis
- Allelic profiling and sequence type assignment
- Population structure analysis
- Database resources and submission

### 4. Serotyping
- Serological classification principles
- In silico serotyping methods
- Species-specific typing schemes
- Clinical relevance of serotypes
- Concordance with traditional methods

### 5. Antimicrobial Resistance Detection
- Resistance gene databases and resources
- Sequence-based resistance prediction
- Point mutations and resistance mechanisms
- Phenotype-genotype correlation
- Clinical interpretation guidelines

### 6. Mobile Genetic Elements
- Plasmid detection and characterization
- Integron structure and analysis
- Transposon identification
- Horizontal gene transfer mechanisms
- Evolution of resistance dissemination

## Tools and Software

### Assembly Tools
- **SPAdes** - De novo genome assembler
- **Unicycler** - Hybrid assembly pipeline
- **Flye** - Long-read assembly
- **QUAST** - Assembly quality assessment

### Annotation Tools
- **Prokka** - Automated prokaryotic annotation
- **RAST** - Rapid Annotation using Subsystem Technology
- **NCBI PGAP** - Prokaryotic Genome Annotation Pipeline
- **Bakta** - Rapid bacterial genome annotation

### Typing Tools
- **mlst** - Multi-locus sequence typing
- **SeqSero2** - Salmonella serotype prediction
- **SeroFinder** - Serotype identification
- **cgMLSTFinder** - Core genome MLST

### Resistance Analysis Tools
- **ABRicate** - Mass screening of contigs for antimicrobial resistance
- **ResFinder** - Identification of acquired antimicrobial resistance genes
- **PointFinder** - Detection of chromosomal point mutations
- **RGI** - Resistance Gene Identifier

### Mobile Element Tools
- **PlasmidFinder** - Plasmid identification
- **IntFinder** - Integron detection
- **ISfinder** - Insertion sequence identification
- **MOB-suite** - Plasmid reconstruction and typing

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

### Exercise 2: Genome Annotation (45 minutes)
Annotate assembled genomes and explore gene content.

```bash
# Automated annotation with Prokka
prokka assembly_output/scaffolds.fasta --outdir annotation/ --genus Mycobacterium --species tuberculosis

# Extract protein sequences
grep -c "^>" annotation/*.faa  # Count proteins

# Functional classification
python3 scripts/analyze_annotations.py annotation/*.gff
```

### Exercise 3: MLST and Serotyping (45 minutes)
Perform molecular typing and serological classification.

```bash
# MLST typing
mlst assembly_output/scaffolds.fasta > mlst_results.txt

# Serotyping (example for E. coli)
python3 serotypefinder/serotypefinder.py -i assembly_output/scaffolds.fasta -o serotype_results/

# Visualize results
cat mlst_results.txt
cat serotype_results/results_tab.txt
```

### Exercise 4: AMR Gene Detection (60 minutes)
Identify antimicrobial resistance genes and predict phenotypes.

```bash
# Screen for resistance genes with ABRicate
abricate --db resfinder assembly_output/scaffolds.fasta > resistance_genes.txt
abricate --db card assembly_output/scaffolds.fasta > card_results.txt

# Point mutation analysis
abricate --db pointfinder assembly_output/scaffolds.fasta > point_mutations.txt

# Summarize resistance profile
python3 scripts/summarize_resistance.py resistance_genes.txt card_results.txt
```

### Exercise 5: Mobile Genetic Element Analysis (45 minutes)
Identify and characterize mobile genetic elements.

```bash
# Plasmid detection
abricate --db plasmidfinder assembly_output/scaffolds.fasta > plasmids.txt

# Integron analysis
python3 integron_finder/integron_finder.py assembly_output/scaffolds.fasta

# Comprehensive mobile element screen
mob_recon --infile assembly_output/scaffolds.fasta --outdir mob_results/
```

## Key Concepts

### Assembly Quality Metrics
| Metric | Good Assembly | Poor Assembly | Action |
|--------|---------------|---------------|--------|
| N50 | >50 kb | <10 kb | Optimize parameters |
| Contigs | <100 | >500 | Check contamination |
| Genome size | Expected ±10% | >20% difference | Review input data |
| Coverage | >50x | <20x | Sequence more |

### MLST Interpretation
- **Sequence Types (STs)**: Unique allelic combinations
- **Clonal Complexes**: Related STs sharing alleles
- **Population structure**: Understanding strain relationships
- **Epidemiological significance**: Linking to clinical outcomes

### Resistance Mechanisms
| Mechanism | Examples | Detection Method |
|-----------|----------|------------------|
| Enzymatic | β-lactamases | Gene presence |
| Efflux pumps | AcrAB-TolC | Gene/regulation |
| Target modification | rRNA methylation | Point mutations |
| Membrane impermeability | Porin loss | Sequence analysis |

### Mobile Element Types
- **Plasmids**: Autonomous replicating elements
- **Integrons**: Gene capture and expression platforms
- **Transposons**: Mobile sequences within genomes
- **Insertion sequences**: Simple transposable elements

## Assessment Activities

### Individual Analysis
- Complete genome assembly workflow
- Perform quality assessment and interpretation
- Execute MLST typing and serotype prediction
- Identify resistance genes and mobile elements
- Generate comprehensive characterization report

### Group Discussion
- Compare assembly strategies and results
- Evaluate MLST typing reliability
- Discuss resistance gene interpretation
- Analyze mobile element distribution patterns

## Common Challenges

### Assembly Issues
```bash
# Low coverage assemblies
spades.py --careful --cov-cutoff 5 -1 R1.fastq -2 R2.fastq -o low_cov_assembly/

# Contamination removal
# Remove contaminant contigs based on taxonomy
seqtk subseq scaffolds.fasta clean_contigs.txt > clean_assembly.fasta
```

### MLST Problems
```bash
# Novel alleles
mlst --novel alleles.fasta assembly.fasta

# Missing loci
# Check for fragmented assemblies or contamination
mlst --debug assembly.fasta
```

### Resistance Prediction Issues
- **Low identity matches**: Consider distant homologs
- **Truncated genes**: Check assembly quality
- **Novel variants**: Manual curation needed
- **False positives**: Validate with phenotypic data

### Mobile Element Challenges
```bash
# Complex rearrangements
# Use long-read sequencing when available
flye --nano-raw long_reads.fastq --out-dir flye_assembly/

# Plasmid circularization
# Check for overlapping ends
python3 scripts/check_circular.py plasmid_contigs.fasta
```

## Clinical Applications

### Routine Surveillance
- Rapid species identification and typing
- Resistance profile determination
- Outbreak detection and investigation
- Infection control decision support

### Public Health Applications
- Antimicrobial resistance monitoring
- Vaccine development support
- Epidemiological investigations
- Policy development and guidelines

## Resources

### Assembly Resources
- [SPAdes Manual](http://cab.spbu.ru/software/spades/)
- [QUAST Documentation](http://quast.sourceforge.net/)
- [Assembly Best Practices](https://github.com/rrwick/Perfect-bacterial-genome-tutorial)

### Typing Resources
- [PubMLST Database](https://pubmlst.org/)
- [MLST.net](https://www.mlst.net/)
- [SerotypeFinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/)

### AMR Resources
- [CARD Database](https://card.mcmaster.ca/)
- [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/)
- [NCBI AMRFinderPlus](https://github.com/ncbi/amr)

### Mobile Element Resources
- [PlasmidFinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/)
- [ISfinder Database](https://isfinder.biotoul.fr/)
- [IntFinder Documentation](https://integronfinder.readthedocs.io/)

## Looking Ahead

**Day 4 Preview**: Comparative Genomics including:
- Pangenome analysis methods
- Phylogenomic reconstruction from core genome SNPs
- Tree construction and visualization
- Evolutionary relationship inference

---

**Key Learning Outcome**: Genome assembly, MLST, serotyping, AMR detection, and mobile genetic element analysis provide comprehensive characterization capabilities essential for clinical genomics and public health surveillance.