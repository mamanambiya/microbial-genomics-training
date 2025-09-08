# Beginner's Microbial Genomics Pipeline

A simple, well-documented Nextflow pipeline perfect for learning bioinformatics workflow management.

## What This Pipeline Does

This pipeline takes raw sequencing data from bacteria and:

1. **Checks data quality** (FastQC)
2. **Cleans the data** (Trimmomatic)
3. **Assembles the genome** (SPAdes)
4. **Evaluates assembly quality** (QUAST)

## Quick Start

### Prerequisites

- Nextflow installed
- Docker or Singularity available
- Some bacterial sequencing data (paired-end FASTQ files)

### Running the Pipeline

```bash
# Basic run
nextflow run beginner-pipeline.nf

# With your own data
nextflow run beginner-pipeline.nf --reads 'my_data/*_{R1,R2}.fastq.gz'

# Custom output directory
nextflow run beginner-pipeline.nf --reads 'data/*_{R1,R2}.fastq.gz' --outdir my_results
```

### Expected File Structure

Your data should be organized like this:

```
data/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
├── sample2_R2.fastq.gz
└── sample3_R1.fastq.gz
└── sample3_R2.fastq.gz
```

## Understanding the Results

After the pipeline completes, you'll find:

```
results/
├── 01_fastqc/          # Quality control reports (open HTML files in browser)
├── 02_trimmed/         # Cleaned sequencing reads
├── 03_assembly/        # Assembled genomes (FASTA files)
└── 04_quast/           # Assembly quality reports
```

### Key Files to Check

1. **FastQC reports** (`01_fastqc/*.html`): Open in web browser to see data quality
2. **Assembly files** (`03_assembly/*_contigs.fasta`): Your assembled genomes
3. **QUAST reports** (`04_quast/*/report.html`): Assembly quality metrics

## Troubleshooting

### Common Issues

**Pipeline fails immediately:**
- Check that Nextflow is installed: `nextflow -version`
- Verify Docker is running: `docker ps`
- Make sure your file paths are correct

**No input files found:**
- Check your file naming pattern matches `*_{R1,R2}.fastq*`
- Verify files exist: `ls data/*_{R1,R2}.fastq*`
- Try absolute paths if relative paths don't work

**Out of memory errors:**
- Reduce the number of samples you're processing
- Check available system memory
- Consider using a high-performance computing cluster

### Getting Help

1. Check the Nextflow log: `cat .nextflow.log`
2. Look at process-specific logs in the `work/` directory
3. Ask for help on the course Slack channel
4. Review the Nextflow documentation: https://www.nextflow.io/docs/

## Learning Exercises

### Beginner Level
1. Run the pipeline with the provided test data
2. Examine the FastQC reports - what do they tell you?
3. Compare assembly statistics between different samples

### Intermediate Level
1. Modify the trimming parameters in the pipeline
2. Add a new process to count genes in the assemblies
3. Create a summary report combining all sample results

### Advanced Level
1. Add conditional logic to handle single-end reads
2. Implement error handling and retry mechanisms
3. Create a configuration file for different computing environments

## Pipeline Details

### Software Used
- **FastQC v0.11.9**: Quality control for sequencing data
- **Trimmomatic**: Read trimming and quality filtering
- **SPAdes**: Genome assembly for bacterial genomes
- **QUAST**: Assembly quality assessment

### Resource Requirements
- **CPU**: 2-4 cores per sample recommended
- **Memory**: 8GB minimum for assembly step
- **Storage**: ~5GB per sample for intermediate files
- **Time**: 30-60 minutes per sample (depends on data size)

### Container Images
All software runs in Docker containers for reproducibility:
- `biocontainers/fastqc:v0.11.9`
- `staphb/trimmomatic:latest`
- `staphb/spades:latest`
- `staphb/quast:latest`

## Next Steps

Once you're comfortable with this pipeline:

1. **Explore nf-core**: Try community-built pipelines like `nf-core/bacass`
2. **Add annotation**: Include gene prediction with Prokka
3. **Comparative analysis**: Add phylogenetic analysis steps
4. **Customize**: Modify the pipeline for your specific research needs

## Support

This pipeline is part of the CIDRI-Africa Microbial Genomics Training Course.

- **Course website**: https://cidri-africa.github.io/microbial-genomics-training/
- **Issues**: Create an issue in the course GitHub repository
- **Questions**: Ask on the course Slack channel

Happy learning! 🧬
