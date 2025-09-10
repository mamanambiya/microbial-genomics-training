# Metagenomics in Clinical & Public Health

## 1. Overview

This section introduces **core principles of metagenomics**. Metagenomics is the study of genetic material recovered directly from environmental or clinical samples, allowing 
analysis of **entire microbial communities** without the need for cultivation. This section includes conceptual notes, rationale for each step, practical commands, 
and a **mini toy dataset exercise**. Unlike genomics (WGS) which focuses on analyzing individual genomes (such as bacterial isolate), metagenomics studies the collective genomes or markers from microbial communities.

### Genomics vs Metagenomics


| Description | Whole Genome Sequencing (WGS) | Metagenomics |
|-------------|-------------------------------|--------------|
| **Umbrella term** | Microbial Genomics | Microbial Profiling |
| **Scope** | Sequencing of the entire genome of a **single isolate** (pure culture) | Sequences **all genetic material** in mixed community (without culturing)|
| **Applications** | Comparative genomics, AMR/virulence genes/MGE, outbreak tracking | Community profiling, functional potential, pathogen detection in mixed samples, ecology studies (diversity studies) |
| **Features** | Organism is known \& isolated before sequencing | Captures both known and unknown microbes from environment or host |

### Metagenomic Strategies: Shotgun vs 16S rRNA

| Description | Shotgun Metagenomics | 16S rRNA Amplicon Gene Sequencing |
|-------------|-------------------------------|--------------|
| **Definition** | Random sequencing of **All DNA** in a sample | **Targeted amplicon sequencing** of the 16 S rRNA gene |
| **Resolution** | Species- and strain-level, **functional genes** (AMR, metabolism, plasmids, MGEs) | Genus-level (sometimes species-level), **no direct functional** information |
| **Merits** | Comprehensive (taxonomy + function), **detects viruses and fungi** |  Limited to taxonomic resolution, can't detect fungi/viruses, lacks functional insights, **good for bacterial** surveys |
| **Demerits** | More expensive, higher computational load | Cost-effective, standardized |

---

## 2. Workflow: From Sequencing to Interpretation

### Step 1. Study Design & Sequencing

- **Why**: Sequencing depth, read length, and platform choice directly influence resolution of taxa, detection of low-abundance organisms, and assembly quality. For example, if you aim to do strain-level tracking of low abundance organisms, increase depth or perform targeted enrichment.
- **Example**:
  - Illumina (short reads): accurate, cost-effective, widely used for clinical metagenomics.
  - ONT/PacBio (long reads): useful for resolving repeats, plasmids, MGEs.

---

### Step 2. Quality Control

- **Why**: To reduce the effects of sequencing errors, remove adapters, and low-quality 
reads which may reduce false positives.
- **Tools**: `fastqc`, `multiqc`

**What to think about?**
- Which parts of the data are flagged as potentially problematic? GC% content of the dataset.
NB: We are dealing with mixed community and organisms, so its difficult to  have a perfect normal distribution around the average.
As such we consider this as normal given the biology of the system.
- Does the sequencing length matches the libraries used? If sequences are shorter than expected, are adapters
 a concern?
- Are adapters and/or barcodes removed?
  - Look at the Per base sequence content to diagnose this.
- Is there unexpected sequence duplication?
  - This can occur when low-input library preparations are used.
- Are over-represented k-mers present?
  - This can be a sign of adapter and barcode contamination.


---

### Step 3. Trimming & Filtering

- **Why**: Removes adapters (which can cause false alignments \& false taxonomic assignments), low-quality bases (increase error rates), and very short 
reads (to standardize read lengths) that bias downstream analysis.
- **Tool**: `fastp`, `trimmomatic`, `BBMap`, `sickle`, `cutadapt`, etc.


```bash
#!/bin/bash

# Load required modules/tools
module load trimmomatic

# Define input and output dir
indata="/data/users/user24/metagenomes/"
wkdir="/data/users/${USER}/metagenomics/shotgun/"
outtrimmomatic=${wkdir}"/data_analysis/02_trimmomatic"

mkdir -p ${outtrimmomatic}

# Trimming with Trimmomatic
for file in `ls ${data}*1.fastq.gz`
do
  sample=$(basename ${file} _1.fastq.gz)
  trimmomatic PE -threads 8 \
      ${indata}${sample}_R1.fastq.gz {indata}${sample}_R2.fastq.gz \
      ${outtrimmomatic}${sample}_R1.fastq.gz ${outtrimmomatic}${sample}_R1_unpaired.fastq.gz \
      ${outtrimmomatic}${sample}_R2.fastq.gz ${outtrimmomatic}${sample}_R2_unpaired.fastq.gz \
      ILLUMINACLIP:adapters.fa:2:30:10 \
      LEADING:3 TRAILING:3 \
      SLIDINGWINDOW:4:20 \
      MINLEN:50
done
echo "Trimming with Trimmomatic completed"
```
Parameter |	Type |	Description
----------|------|-----------------------------
PE |	positional |	Specifies whether we are analysing single- or paired-end reads
-threads 2 |	keyword |	Specifies the number of threads to use when processing
-phred33 |	keyword |	Specifies the fastq encoding used
${indata}${sample}_R1.fastq.gz {indata}${sample}_R2.fastq.gz |	positional 	| The paired forward and reverse reads to trim
ILLUMINACLIP:adapters.fa:2:30:10 |	positional |	Adapter trimming allowing for 2 seed mismatch, palindrome clip score threshold of 30, and simple clip score threshold of 10
SLIDINGWINDOW:4:20 	| positional 	|Quality filtering command. Analyse each sequence in a 4 base pair sliding window and then truncate if the average quality drops below Q20
MINLEN:50 	| positional |	Length filtering command. Discard sequences that are shorter than 80 base pairs after trimming

> **Note:** The trimming parameters in `trimmomatic` are processed in the order they are specified. For instance, 

---

```bash
for file in `ls ${data}*1.fastq.gz`
do
  sample=$(basename ${file} _1.fastq.gz)
  trimmomatic PE -threads 8 \
        ${indata}${sample}_R1.fastq.gz {indata}${sample}_R2.fastq.gz \
        ${outtrimmomatic}${sample}_R1.fastq.gz ${outtrimmomatic}${sample}_R1_unpaired.fastq.gz \
        ${outtrimmomatic}${sample}_R2.fastq.gz ${outtrimmomatic}${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:adapters.fa:2:30:10 \
        LEADING:3 TRAILING:3 \
        MINLEN:50 \
        SLIDINGWINDOW:4:20 
done
```
means we remove sequences shorter than 50 bps and then qualiyty trim, thus if a sequence is trimmed to a length shorter than 50bps after trimming, the `MINLEN` filtering does not execute a second time.

```bash
# Delete unnecessary files
rm -rf ${outtrimmomatic}*${sample}*_unpaired.fastq.gz

# Quality checking
module load fastqc multiqc
fastqc ${outtrimmomatic}* -o ${outtrimmomatic}
multiqc  ${outtrimmomatic} -o  ${outtrimmomatic}
```


---


### Step 4. Deduplication & Host DNA Removal

- **Why**:
- Often metagenomes are obtained from host-associated microbial communities. As a result, they contain significant amount of host DNA which may interfere with microbial analysis
  and create privacy concerns.
- Specifically any studies invloving human subjects or samples derived from Taonga species.
- Although several approaches are used for this, the most popular is to map reads to a reference genome (includin human genome). That is remove all reads that map to the reference of the dataset.
  
- **Tools**: `clumpify` (dedup), `BBMap` (bbmap.sh),  `bowtie2` or `bwa mem` (host removal).

```bash
#!/bin/bash
set -e  # Exit on error

# Define dirs
genome_dir="/data/users/user24/refs/human_reference/"

## Remove human contamination using BWA and SAMtools
## Create dir
mkdir -p ${genome_dir}
# Download the Human refence genome
cd ${genome_dir}

# Configuration
GENOME_VERSION="GRCh38"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38"

echo "=== Downloading Human Reference Genome (${GENOME_VERSION}) ==="

# Option 1: Download complete genome assembly (recommended for most applications)
echo "Downloading complete genome assembly..."
wget -c "${BASE_URL}/GCA_000001405.15_GRCh38_genomic.fna.gz" \
     -O "GRCh38_genomic.fna.gz"
# Decompress
echo "Decompressing genome file..."
gunzip -f GRCh38_genomic.fna.gz

# Option 2: Download chromosome-only version (excludes contigs/scaffolds)
echo "Downloading chromosome-only version..."
wget -c "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" \
     -O "hg38_chromosomes_only.fa.gz"

gunzip -f hg38_chromosomes_only.fa.gz
```

```bash
#!/bin/bash

# Load modules
module load bwa
module load bowtie2

echo "=== Building BWA Index ==="
# Build BWA index for alignment (required for host removal, one time set-up)
bwa index GRCh38_genomic.fna

echo "=== Verifying Download ==="
# Check file integrity
echo "Genome file size:"
ls -lh *.fna

echo "Number of sequences:"
grep -c ">" GRCh38_genomic.fna

echo "First few sequence headers:"
grep ">" GRCh38_genomic.fna | head -10

echo "=== Download Complete ==="
echo "Reference genome files are ready in: $(pwd)"
echo "Main genome file: GRCh38_genomic.fna"
echo "Chromosome-only file: hg38_chromosomes_only.fa"
echo ""
echo "Files generated:"
echo "- GRCh38_genomic.fna (main reference)"
echo "- GRCh38_genomic.fna.amb, .ann, .bwt, .pac, .sa (BWA index)"
echo "- GRCh38_genomic.fna.fai (samtools index)"

# Align reads to human genome
for file in `ls ${indata}*1.fastq.gz`
do
  sample=$(basename ${file} _1.fastq.gz)
  bowtie2 -x ${refs}human_index -1 {sample}${data}_trimmed_R1.fastq.gz -2 sample_trimmed_R2.fastq.gz \
    --un-conc sample_host_removed.fastq.gz -S /dev/null


  bwa mem -t 8 human_reference.fasta \
      trimmed_R1_paired.fastq.gz trimmed_R2_paired.fastq.gz \
      | samtools view -b -f 4 - > unmapped_reads.bam
  
  # Convert unmapped reads back to FASTQ
  samtools fastq -1 dehosted_R1.fastq.gz -2 dehosted_R2.fastq.gz unmapped_reads.bam

echo "Host removal completed"
```

---

### Step 5. Taxonomic Profiling (Read-Based)

- **Why**: Directly assigns taxonomy without assembly, faster and less computationally intensive.
- Research question: Who is in this sample?
- **Tools**: Kraken2 + Bracken (abundance estimation), MetaPhlAn/ MetaSpades + Mapping + MetaBAT2 (or CONCOT/MaxBin2, DAS Tool) + checkM + GTDB-tk + ANI (species-level contamination)
- Reconstruct rRNA using `PhyloFlash - EMIRGE` (by leveraging the expansive catalogue of 16S rRNA genes available in databases such as SILVA in order to subset reads and then reconstruct the full-length gene).

```bash
kraken2 --db kraken2_db --paired sample_host_removed.1.fq.gz sample_host_removed.2.fq.gz \
  --output kraken2_output.txt --report kraken2_report.txt
bracken -d kraken2_db -i kraken2_report.txt -o bracken_species.txt -r 150 -l S
```

---

### Step 6. Functional Profiling (Read-Based)

- **Why**: Identifies pathways/genes present without needing assembly.
- **Tool**: HUMAnN

```bash
humann --input sample_host_removed.fastq.gz --output humann_out/
```

---

### Step 7. Assembly-Based Profiling

- **Why**: Reconstructs contigs/MAGs b enables strain-level analysis, discovery of new genes, plasmids, MGEs.
- **Tools**: **MEGAHIT** or **metaSPAdes**

#### MEGAHIT

```bash
megahit -1 sample_host_removed.1.fq.gz -2 sample_host_removed.2.fq.gz -o megahit_out/
```


```bash
spades.py --meta -1 sample_host_removed.1.fq.gz -2 sample_host_removed.2.fq.gz -o metaspades_out/
```

**MEGAHIT vs metaSPAdes:**

- **MEGAHIT Advantages:**
  - Extremely fast, lower memory usage.
  - Scales better for very large datasets (e.g., population metagenomes).
- **MEGAHIT Disadvantages:**
  - May produce slightly shorter contigs than metaSPAdes.
- **metaSPAdes Advantages:**
  - Produces higher-quality, longer assemblies (useful for MAG recovery).
- **metaSPAdes Disadvantages:**
  - Requires more RAM and CPU time.

---

### Step 8. Mapping & Binning

- **Why**: Group contigs into MAGs, quantify abundances.
- **Tools**: `bowtie2`, `samtools`, `MetaBAT2` / `CONCOT` / `MaxBin2`, `DAS Tool`

```bash
bowtie2 -x megahit_out/final.contigs.fa -1 sample_host_removed.1.fq.gz -2 sample_host_removed.2.fq.gz | samtools sort -o aln.bam
metabat2 -i megahit_out/final.contigs.fa -a depth.txt -o bins_dir/bin
```

---

### Step 9. MAG Quality Control

- **Why**: Ensures completeness & contamination are acceptable.
- **Tool**: CheckM

```bash
checkm lineage_wf bins_dir/ checkm_out/
```

---

### Step 10. Annotation & Specialized Analyses

- **Why**: Identify AMR, virulence, plasmids, metabolic capacity.
- **Tools**: Prokka, Bakta, AMRFinderPlus, ABRicate, DRAM
- Use **DRAM** in place of Prokka/Bakta for **functional annotation of MAGs**.
- DRAM produces detailed metabolic profiles, pathway reconstruction, and microbial ecology insights.

```bash
DRAM.py annotate -i bins_dir/ -o dram_out/ --threads 16
```

---

### Step 11. Abundance Estimation & Visualization

- **Why**: Quantifies taxa/genes  links to clinical or epidemiological metadata.
- **Tools**: CoverM, Krona, R for plots.

```bash
coverm contig --bam-files aln.bam --reference megahit_out/final.contigs.fa --methods tpm > coverm_tpm.tsv
```

---

### Step 12. Reporting & Reproducibility

- **Why**: Essential for public health applications.
- Use **Nextflow + Singularity/Conda** for reproducible pipelines.
- Summarize results in Excel or RMarkdown reports.

---


## Shotgun metagenomics

For this training we will use the data in `/data/users/user29/metagenomes/shotgun/` for shotgun metagenomics. It is important to note that, one can download shotgun metagenome sequences from NCBI-SRA using `ncbi-tools`. Install `ncbi-tools` and run

```bash
## Fetch the data from NCBI-SRA
# fasterq-dump SRR13827118 --progress --threads 8

## Compress the files
gzip SRR13827118*.fastq
```

---

```bash
# Load modules
module load fastqc
module load multiqc

# Data directory
indata="/data/users/user29/metagenomes/shotgun/"

## Create Dir
mkdir -p /data/users/$USER/metagenomes/shotgun/scripts /data/users/$USER/metagenomes/shotgun/logs

# create a samplesheet.csv with three columns samplename,fastq_1,fastq_2
python /data/users/user24/metagenomes/shotgun/scripts/samplesheet_generator.py ${indata} \
  /data/users/${USER}/metagenomes/shotgun/samplesheet.csv

# Run raw QC
nextflow run /data/users/$USER/metagenomes/shotgun/scripts/qc_pipeline_v1.nf \
  --input /data/users/${USER}/metagenomes/shotgun/samplesheet.csv \
  --outdir /data/users/${USER}/metagenomes/shotgun/results/rawfastqc
```

### Quality assessment
```bash
module load bbmap
# Generate statistics for filtered SPAdes assembly
stats.sh in=spades_assembly/spades.fna
```

```bash
# Assuming your contigs have .fa as suffix
CONTIG_DIR=""
OUT_DIR=""
mkdir ${OUT_DIR}
for contig in ${CONTIG_DIR}*.fa
module load bowtie2
do
    SAMPLE=$(basename ${contig} .fa)
    R1=${rawdata}/${SAMPLE}_R1_001.fastq.gz
    R2=${rawdata}/${SAMPLE}_R2_001.fastq.gz

    echo "=== Processing ${SAMPLE} ==="

    CONTIGS=${CONTIG_DIR}${SAMPLE}".fa"
    mkdir -p ${OUT_DIR}/${SAMPLE}/"align"

    ## Step 2: Bowtie2 index
    #echo "Building Bowtie2 index..."
    bowtie2-build ${CONTIGS} ${OUT_DIR}/${SAMPLE}/align/${SAMPLE}

    ## Step 3: Align reads
    #echo "Aligning reads..."
    bowtie2 -x ${OUT_DIR}/${SAMPLE}/align/${SAMPLE} -1 "$R1" -2 "$R2" \
    #  -S ${OUT_DIR}/${SAMPLE}/align/${SAMPLE}.sam -p $THREADS

    #samtools view -Sb ${OUT_DIR}/${SAMPLE}/align/${SAMPLE}.sam | \
    #  samtools sort -o ${OUT_DIR}/${SAMPLE}/align/${SAMPLE}.bam

    #samtools index ${OUT_DIR}/${SAMPLE}/align/${SAMPLE}.bam

    ## Step 4: Depth file
    #echo "Computing depth..."
    #jgi_summarize_bam_contig_depths --outputDepth ${OUT_DIR}/${SAMPLE}/${SAMPLE}_depth.txt \
    #  ${OUT_DIR}/${SAMPLE}/align/${SAMPLE}.bam

    ## Step 5: Binning
    #mkdir -p ${OUT_DIR}/${SAMPLE}/bins

    #echo "Binning with MetaBAT2..."
    #metabat2 -i ${CONTIGS} -a ${OUT_DIR}/${SAMPLE}/${SAMPLE}_depth.txt \
    #  -o ${OUT_DIR}/${SAMPLE}/bins/${SAMPLE}_bin -t $THREADS
done
```

### Bacteriophage analysis

- Useful [protol](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=3) that requires `virsorter2`, `checkV` and `DRAMv`
  
#### nfcore/mag Pipeline <https://nf-co.re/mag/4.0.0>

- Use an [nf-core/mag](https://nf-co.re/mag/4.0.0/) pipeline for assembly, binning and annotation of metagenomes, [github repository](https://github.com/nf-core/mag/tree/4.0.0).
- This pipeline works for short- and/or long-reads.
- ✅**Key Features:**
  - **Preprocessing:**
     - Short reads: fastp, Bowtie2, FastQC
     - Long reads: Porechop, NanoLyse, Filtlong, NanoPlot
  - **Assembly:**
     - Short reads: MEGAHIT, SPAdes
     - Hybrid: hybridSPAdes
  - **Binning:**
     - Tools: MetaBAT2, MaxBin2, CONCOCT, DAS Tool
     - Quality checks: BUSCO, CheckM, GUNC
 - **Taxonomic Classification:**
    - Tools: GTDB-Tk, CAT/BAT
    - Co-assembly and co-abundance:
    - Supports sample grouping for co-assembly and binning
- 📦 **Reproducibility:**
- Uses Nextflow DSL2, Docker/Singularity containers
- Fully portable across HPC, cloud, and local systems
- Includes test datasets and CI testing
- It requires a sample sheet in `csv`  format with five columns: sample,group,short_reads_1,short_reads_2,long_reads
- Assuming all your raw reads (short- and long-reads) are in the same folder, run the python script:

```bash
proj="/data/users/${USER}/metagenomes/shotgun/"
# Create required directories
mkdir -p ${proj}scripts ${proj}logs
# Get the script to create samplesheet
cp /data/users/user24/metagenomes/shotgun/scripts/generate_mag_samplesheet.py /data/users/${USER}/metagenomes/shotgun/scripts/

# Create samplesheet for running mag
python3 /data/users/${USER}/metagenomes/shotgun/scripts/generate_mag_samplesheet.py /data/users/user29/metagenomes/shotgun/ mag-samplesheet.csv
# Create submission script
nano data/users/${USER}/metagenomes/shotgun/scripts/mag-nf_submit.sh

#!/bin/bash
#SBATCH --job-name='mag'
#SBATCH --time=24:00:00
#SBATCH --mem=128g
#SBATCH --ntasks=16
#SBATCH --output=/data/users/user24/metagenomes/shotgun/logs/nfcore-mag-stdout.log
#SBATCH --error=/data/users/user24/metagenomes/shotgun/logs/nfcore-mag-stderr.log
#SBATCH --mail-user=ephie.geza@uct.ac.za

proj="/data/users/user24/metagenomes/shotgun/"

module load nextflow/25.04.6
#### Unload JAVA 18 as it doesn't work and load JAVA 17
module unload java/openjdk-18.0.2
module load java/openjdk-17.0.2
####Unset conflicting environment variables (optional but recommended)
unset JAVA_CMD
unset JAVA_HOME

#### Run pipeline
nextflow run ${proj}mag \
      --input ${proj}mag-samplesheet.csv \
      --outdir ${proj}nfcore-mag \
      -w ${work}work/nfcore-mag \
      -profile singularity \
      -resume --skip_gtdbtk

## Save and submit
```

### [nf-core/funcscan](https://nf-co.re/funcscan/2.1.0/)

Using contigs to screen for functional and natural gene sequences

### Viralgenie [nf-core/viralmetagenome](https://nf-co.re/viralmetagenome/0.1.2/)

### [nf-core/taxprofiler](https://nf-co.re/taxprofiler/1.2.4/)
## Targeted Metagenomics 

### [nf-core/ampliseq](https://nf-co.re/ampliseq/2.14.0/)

