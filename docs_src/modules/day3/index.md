# Day 3: Pathogen Genomics

**Duration**: Full day (09:00-13:00)  
**Focus**: Pathogen-specific analysis, AMR detection, and comparative genomics

## Overview

Day 3 dives deep into pathogen-specific genomic analysis using real clinical isolates. We'll explore the unique characteristics of major pathogens and learn to detect antimicrobial resistance genes and mobile genetic elements.

## Learning Objectives

By the end of Day 3, you will be able to:

- Perform pathogen-specific genomic characterization
- Detect and interpret antimicrobial resistance profiles
- Identify mobile genetic elements and their role in resistance spread
- Apply comparative genomics approaches to outbreak investigation
- Understand pathogen typing methods (MLST, cgMLST, SNP typing)

## Schedule

| Time | Topic | Instructor |
|------|-------|------------|
| 09:00-09:45 | [*M. tuberculosis* Analysis](mtb-analysis.md) | Ephifania Geza |
| 09:45-10:30 | [*V. cholerae* Analysis](vibrio-analysis.md) | Arash Iranzadeh |
| 10:30-10:45 | *Coffee Break* | |
| 10:45-11:30 | [Antimicrobial Resistance Detection](amr-detection.md) | Sindiswa Lukhele |
| 11:30-12:15 | Mobile Genetic Elements | Ephifania Geza |
| 12:15-13:00 | Comparative Analysis Workshop | All Instructors |

## Key Pathogens

### *Mycobacterium tuberculosis*
- **Genome Size**: ~4.4 Mb
- **Key Features**: Slow growth, drug resistance, lineage diversity
- **Analysis Focus**: Drug resistance mutations, lineage typing, transmission clusters

### *Vibrio cholerae*
- **Genome Size**: ~4.0 Mb (2 chromosomes)
- **Key Features**: Pandemic potential, toxin genes, environmental survival
- **Analysis Focus**: Pathogenicity islands, outbreak tracking, environmental adaptation

### *Escherichia coli* & *Shigella sonnei*
- **Genome Size**: ~5.0 Mb
- **Key Features**: Diverse pathotypes, plasmid-mediated resistance
- **Analysis Focus**: Virulence factors, resistance plasmids, outbreak investigation

## Tools and Databases

### Resistance Detection
- **ABRicate** - Screening for resistance genes
- **ResFinder** - Comprehensive resistance gene database
- **CARD** - Comprehensive Antibiotic Resistance Database
- **PointFinder** - Point mutations in resistance genes

### Typing Methods
- **mlst** - Multi-locus sequence typing
- **chewBBACA** - cgMLST analysis
- **SNP-sites** - Core SNP extraction
- **ParSNP** - Rapid core genome SNP typing

### Specialized Tools
- **TB-Profiler** - *M. tuberculosis* drug resistance prediction
- **Mykrobe** - Species identification and resistance prediction
- **SISTR** - *Salmonella* in silico typing resource

## Hands-on Exercises

### Exercise 1: TB Drug Resistance Analysis (60 minutes)
Analyze *M. tuberculosis* genomes for drug resistance mutations.

```bash
# TB-Profiler analysis
tb-profiler profile -1 TB_sample_R1.fastq.gz -2 TB_sample_R2.fastq.gz \
    -p TB_sample --txt --csv

# Generate summary report
tb-profiler collate -d results/ -p batch_analysis
```

### Exercise 2: *V. cholerae* Outbreak Investigation (45 minutes)
Investigate a cholera outbreak using comparative genomics.

```bash
# Core genome alignment
parsnp -g reference.gbk -d genomes/ -p 4 -o parsnp_output/

# Generate phylogenetic tree
FastTree -nt -gtr parsnp_output/parsnp.xmfa > outbreak_tree.newick
```

### Exercise 3: AMR Gene Detection (45 minutes)
Screen multiple pathogens for antimicrobial resistance genes.

```bash
# Screen for resistance genes
abricate --db resfinder genome_assembly.fasta > resistance_profile.txt

# Screen multiple databases
for db in resfinder card argannot; do
    abricate --db $db assembly.fasta > ${db}_results.txt
done
```

## Dataset Details

### TB Clinical Isolates
- **n=20** diverse clinical strains
- **Drug resistance**: Mixed MDR/XDR and susceptible strains
- **Lineages**: Representatives from major global lineages
- **Metadata**: Treatment history, geographic origin

### Cholera Outbreak Collection
- **n=15** outbreak isolates + environmental samples
- **Source**: Recent epidemic investigation
- **Timespan**: 6-month outbreak period
- **Geography**: Urban coastal setting

### *E. coli* Surveillance Data
- **n=25** clinical isolates
- **Source**: Hospital surveillance program  
- **Resistance**: ESBL and carbapenem resistance focus
- **Pathotypes**: Mixed UPEC, STEC, and commensal strains

## Key Concepts

### Drug Resistance Mechanisms
- **Target modification**: Mutations in drug targets
- **Drug inactivation**: Enzymatic resistance
- **Efflux pumps**: Active drug removal
- **Reduced permeability**: Altered uptake

### Mobile Genetic Elements
- **Plasmids**: Circular DNA elements
- **Transposons**: Mobile DNA sequences
- **Integrons**: Gene capture systems
- **Prophages**: Integrated viral sequences

### Typing Methods Comparison
| Method | Resolution | Speed | Coverage | Best For |
|--------|------------|-------|----------|----------|
| MLST | Low | Fast | Universal | Species identification |
| cgMLST | High | Medium | Species-specific | Outbreak investigation |
| SNP typing | Very High | Slow | Universal | Fine-scale epidemiology |

## Assessment Activities

### Individual Analysis
- Complete resistance profiling for assigned isolates
- Generate phylogenetic analysis of outbreak data
- Interpret AMR gene presence/absence patterns

### Group Discussion
- Compare resistance mechanisms across pathogens
- Evaluate outbreak investigation strategies
- Discuss public health implications of findings

## Common Challenges

### False Positive Resistance Calls
```bash
# Verify with multiple tools
abricate --db resfinder assembly.fasta
rgi main -i assembly.fasta -o rgi_output -t contig -a BLAST
```

### Incomplete Assemblies
```bash
# Check assembly quality first
quast.py assembly.fasta -o quality_check/
# Consider using long-read data or hybrid assembly
```

### Phylogenetic Artifacts  
```bash
# Mask recombinogenic regions
gubbins aligned_sequences.aln
# Use recombination-free SNPs for tree building
```

## Advanced Topics

### Machine Learning in AMR Prediction
- Feature selection for resistance prediction
- Training datasets and validation
- Integration with genomic epidemiology

### Pan-genome Analysis
- Core vs accessory genome concepts
- Virulence factor distribution
- Population structure analysis

## Resources

### Databases
- [CARD Database](https://card.mcmaster.ca/)
- [ResFinder Database](https://cge.cbs.dtu.dk/services/ResFinder/)
- [PubMLST](https://pubmlst.org/)
- [TB-Profiler Database](https://github.com/jodyphelan/tbdb)

### Key Publications
- Zankari et al. (2012). Identification of acquired antimicrobial resistance genes
- Hunt et al. (2019). ARIBA: rapid antimicrobial resistance genotyping
- Phelan et al. (2019). Integrating informatics tools and portable sequencing

### Online Tools
- [CGE Server](https://cge.cbs.dtu.dk/) - Multiple typing tools
- [Galaxy Project](https://usegalaxy.org/) - Web-based analysis platform
- [BacWGSTdb](https://bacdb.org/) - Bacterial genomics database

## Looking Ahead

**Day 4 Preview**: Advanced workflows and metagenomics including:
- Metagenomics analysis principles
- Workflow management with Nextflow
- Data interpretation and visualization

---

**Remember**: Resistance gene presence doesn't always equal phenotypic resistance. Consider gene expression, mutations, and clinical context in your interpretations!