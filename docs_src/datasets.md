# Training Datasets

## Overview

This course uses carefully curated datasets representing real-world scenarios in microbial genomics and metagenomics. All data has been quality-controlled and prepared for educational use.

## Dataset Categories

### 1. Primary Training Datasets

#### *Mycobacterium tuberculosis* Collection
- **Sample Size**: 20 clinical isolates
- **Geographic Origin**: Global collection (South Africa, India, UK, Peru)
- **Drug Resistance**: Mixed MDR, XDR, and drug-susceptible strains
- **Lineages**: Representatives from Lineages 1-4
- **Sequencing**: Illumina paired-end (2×150bp, 80-120x coverage)
- **Size**: ~2.5 GB
- **Use Cases**: Drug resistance analysis, lineage typing, phylogenetic analysis

#### *Vibrio cholerae* Outbreak Investigation
- **Sample Size**: 15 outbreak isolates + 3 environmental samples
- **Source**: Simulated cholera outbreak (based on real data)
- **Timeframe**: 6-month epidemic period
- **Geographic**: Coastal urban setting
- **Sequencing**: Illumina paired-end (2×150bp, 60-100x coverage)
- **Size**: ~1.8 GB
- **Use Cases**: Outbreak tracking, source attribution, transmission analysis


### 2. Reference Materials

#### Reference Genomes
- High-quality reference assemblies for all target pathogens
- Annotation files (GFF, GenBank formats)
- Resistance gene databases
- **Size**: ~500 MB

#### Validation Datasets
- Known outbreak collections with confirmed epidemiological links
- Quality control standards
- Benchmark datasets for method comparison
- **Size**: ~1.5 GB

## Dataset Access

### File Organization
```
datasets/
├── genomics/
│   ├── mtb/                    # M. tuberculosis isolates
│   └── vibrio/                 # V. cholerae outbreak
├── references/
│   ├── genomes/                # Reference assemblies
│   ├── databases/              # Resistance/virulence databases
│   └── annotations/            # Gene annotations
└── validation/
    ├── benchmarks/             # Method comparison datasets
    └── qc_standards/           # Quality control references
```

### Download Instructions

#### During Course
Data is pre-loaded on course HPC systems:
```bash
# Access course data directory
cd /data/course/datasets/

# Copy to your workspace
cp -r /data/course/datasets/ ~/workspace/
```

#### Post-Course Access
Datasets remain available through:
```bash
# Clone dataset repository
git clone https://github.com/CIDRI-Africa/microbial-genomics-datasets.git

# Download specific collections
wget https://datasets.microbial-genomics.org/mtb_collection.tar.gz
```

## Data Formats

### Raw Sequencing Data
- **Format**: FASTQ (compressed with gzip)
- **Quality**: Phred+33 encoding
- **Naming**: `SampleID_R1.fastq.gz`, `SampleID_R2.fastq.gz`

### Processed Data
- **Assemblies**: FASTA format
- **Annotations**: GFF3, GenBank
- **Alignments**: SAM/BAM format
- **Variants**: VCF format

### Metadata
- **Sample Information**: CSV/TSV format
- **Study Design**: Detailed README files
- **Quality Metrics**: MultiQC reports included

## Metadata Schema

### Genomic Samples
| Field | Description | Example |
|-------|-------------|---------|
| sample_id | Unique identifier | MTB_001 |
| species | Organism name | Mycobacterium tuberculosis |
| collection_date | Sample date | 2023-01-15 |
| location | Geographic origin | Cape Town, South Africa |
| resistance_profile | Known resistance | INH-R, RIF-R |
| sequencing_platform | Technology | Illumina MiSeq |
| coverage_depth | Average coverage | 85x |


## Quality Control

### Pre-processing Standards
- **Quality Score**: Minimum Q30 for 80% of bases
- **Contamination**: <2% non-target DNA
- **Coverage**: Minimum 30x for genomic samples
- **Assembly Quality**: N50 >50kb, <200 contigs

### Validation Procedures
- Species confirmation by 16S rRNA or genome similarity
- Contamination screening with multiple tools
- Assembly quality assessment with standard metrics
- Metadata validation and consistency checking

## Ethical Considerations

### Data Privacy
- All clinical data de-identified according to HIPAA standards
- Geographic information limited to city/region level
- No patient identifiers or medical record linkage possible

### Usage Rights
- Educational use permitted under Creative Commons License
- Commercial use requires separate permission
- Attribution required for publications using these datasets
- Redistribution allowed with proper citation

### Responsible Use
- Data should not be used to identify individuals
- Results should not be used for clinical decision-making
- Sharing outside course requires instructor approval

## Dataset-Specific Notes

### *M. tuberculosis* Collection
- Lineage assignments based on SNP typing
- Drug resistance confirmed by phenotypic testing
- Geographic sampling represents global diversity
- Suitable for phylogeographic analysis

### *V. cholerae* Outbreak
- Temporal sampling allows transmission inference
- Environmental samples included for source attribution
- Metadata includes case demographics and exposure history
- Excellent dataset for outbreak investigation training


## Troubleshooting

### Common Issues

#### File Access Problems
```bash
# Check file permissions
ls -la datasets/
# Fix permissions if needed
chmod -R 755 datasets/
```

#### Corrupted Downloads
```bash
# Verify file integrity
md5sum -c checksums.md5
# Re-download corrupted files
wget -c https://datasets.url/file.tar.gz
```

#### Storage Space Issues
```bash
# Check available space
df -h
# Compress unused files
gzip *.fastq
# Remove temporary files
rm -rf temp/
```

## Citation Information

When using these datasets in publications, please cite:

> CIDRI-Africa Microbial Genomics Training Consortium. (2024). 
> Comprehensive training datasets for microbial genomics and metagenomics education. 
> *Microbial Genomics Education*, Dataset Repository.

Individual dataset citations available in respective README files.

## Support

For dataset-related questions:
- **Technical Issues**: Submit issue on GitHub repository
- **Scientific Questions**: Contact course instructors
- **Access Problems**: Email dataset-admin@cidri-africa.org

## Updates and Versioning

- **Current Version**: v2.1 (September 2025)
- **Update Frequency**: Annually or as needed
- **Change Log**: Available in repository documentation
- **Notification**: Users notified of major updates via email

---

**Remember**: These datasets represent real scientific data and should be treated with appropriate care and respect for the original sample donors and research contexts.