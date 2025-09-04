# De novo Genome Assembly Practical

Before setting up you need to know the current workig directory 
```bash
# Check current working directory
pwd
# List the contents of the working diresctory
ls
# Define output dir
outdir="/data/users/${USER}/data_analysis/assembly/"
outspades=${outdir}"01_spades/"
outquast=${outdir}"02_quast/"
# Create relevant output directories. -p so that it creates parent dir if it doesn't exist
mkdir -p ${outspades}"tb"
mkdir -p ${outspades}"vc"
mkdir -p ${outquast}"tb"
mkdir -p ${outquast}"vc"
```

Delete ${outspades} and ${outquast} and create them using one line of code and not four as above

```bash
# Request for an interactive node (Resources)
## Use this command to retrieve the previously used "srun" command
history | grep "srun"
srun --cpus-per-tasks=16 --mem=128g --time 3:00:00 --job-name "${USER}-assembly" --pty /bin/bash

# # Clear all modules
module purge
# Load the required modules
module load spades/4.2.0
## How can I get the help documentation?
spades.py --help | less
# Lets run spades with a test dataset
/software/bio/spades/4.2.0/bin/spades.py  --test  --careful
```
We  will use trimmed reads and file with IDs that we created during the initial data cleaning process.
Check where the ID files are  and CHANGE DIR to where these are. Use "cd /path/to/tb_IDs"

```bash
## Subset the data
head -2 /data/users/${USER}/data_analysis/tb_IDs > /data/users/${USER}/02_tb_IDs
head -2 /data/users/${USER}/data_analysis/vc_IDs > /data/users/${USER}/02_vc_IDs
# Define dirs
indir="/data/users/${USER}/data_analysis/trimmed_trimmomatic/"
outdir="/data/users/${USER}/data_analysis/assembly-test/"
outspades=${outdir}"01_spades/"
outquast=${outdir}"02_quast/"

mkdir -p ${outspades}tb ${outspades}vc ${outquast}tb ${outquast}vc

cd /data/users/${USER}/
# Run SPAdes for M. tuberculosis
echo "Starting M. tuberculosis assembly at $(date)"

for SAMPLE in $(cat 02_tb_IDs)
do
  echo "[TB|SPADES] $SAMPLE"
  spades.py -1 ${indir}tb/${SAMPLE}_1.fastq.gz \
    -2 ${indir}tb/${SAMPLE}_2.fastq.gz \
    --careful \
    --cov-cutoff auto \
    -t 16 \
    -o ${outspades}tb/${SAMPLE}
done
```
## Notes
- -1 <filename>               file with forward paired-end reads
- -2 <filename>               file with reverse paired-end reads
- --careful                   error and mismatch correction
- -o                          path to output file (directory shouldn't be existing)
- --cov-cutoff auto           Automatic coverage cutoff (usiful mostly for high GC genomes)


```bash
# Perform denovo assembly for all TB ands VC samples
mkdir -p /users/${USER}/scripts /users/${USER}/logs

## Now let's write our submission script
nano /users/${USER}/scripts/assembly_01.sh

#!/bin/bash
#SBATCH --job-name='assembly-${USER}'
#SBATCH --nodes=1 --ntasks=16
#SBATCH --partition=Main
#SBATCH --mem=128GB
#SBATCH --output=/users/${USER}/logs/assembly_01-stdout.txt
#SBATCH --error=/users/${USER}/logs/assembly_01-stdout.txt
#SBATCH --time=12:00:00

# Define dir
indir="/data/users/${USER}/data_analysis/trimmed_trimmomatic/"
outdir="/data/users/${USER}/data_analysis/assembly/"
outspades=${outdir}"01_spades/"
outquast=${outdir}"02_quast/"

# Change into the DIR with "tb_IDs"
cd /data/users/${USER}/data_analysis/
# Run SPAdes for M. tuberculosis
echo "Starting M. tuberculosis assembly at $(date)"
for SAMPLE in $(cat tb_IDs); do
  echo "[TB|SPADES] $SAMPLE"
  spades.py -1 ${indir}tb/${SAMPLE}_1.fastq.gz \
    -2 ${indir}tb/${SAMPLE}_2.fastq.gz \
    --careful --cov-cutoff auto -t 16 \
    -o ${outspades}tb/${SAMPLE}
done
echo "M. tuberculosis assembly completed at $(date)"

for SAMPLE in $(cat vc_IDs); do
  echo "[VC|SPADES] $SAMPLE"
  spades.py -1 ${indir}vc/${SAMPLE}_1.fastq.gz \
    -2 ${indir}vc/${SAMPLE}_2.fastq.gz \
    --careful -t 16 \
    -o ${outspades}vc/${SAMPLE}
done
## SAVE THIS IN NANO BY
# Ctrl + O

# CLOSE THE FILE
# ctrl + x
```
**Before submission of your script, verify if all the input directory and files exist.**

```bash
## SUBMIT THE ASSEMBLY SCRIPT
sbatch /users/${USER}/scripts/assembly_01.sh
```

# Genome Assembly Quality 
Assess the quality of your assemblies to be confident of your downstream analysis. 

## Common Assembly Quality Metrics
- Expected genome size
- Expected GC content %
- Number of contigs
- N50 and L50

## Common Assembly Issues and Solutions
===================================

| 🚨 **Issue** | 🔍 **Possible Cause** | 🛠️ **Solution** |
|--------------|----------------------|-----------------|
| **High number of contigs (>100)** | • Low coverage (<30x) <br> • Poor quality reads <br> • Highly repetitive genome <br> • Contamination | • **Increase sequencing depth** <br> • Better **read trimming/filtering** <br> • Use **hybrid assembly** with long reads <br> • Check contamination with **Kraken2** |
| **Low N50 (<50kb for bacteria)** | • Repetitive sequences <br> • Low coverage <br> • Assembly parameter issues | • Adjust **k-mer sizes** in *SPAdes* <br> • Use `--careful` flag <br> • Try different assemblers (*Unicycler*, *SKESA*) |
| **Assembly size much larger than expected** | • Contamination <br> • Duplication in assembly <br> • Presence of plasmids | • Perform **contamination screening** <br> • Check duplication ratio in **QUAST** <br> • Separate **plasmid sequences** |
| **High number of hypothetical proteins (>50%)** | • Novel organism/strain <br> • Poor annotation database match <br> • Assembly quality issues | • Use specialized databases (**RefSeq**, **UniProt**) <br> • **Manual curation** of key genes <br> • Functional annotation with **KEGG/COG** |


## Running QUAST Tool

```bash
# Run QUAST for a subset of M. tuberculosis dataset
echo "Running QUAST analysis for M. tuberculosis..."

indir="/data/users/${USER}/data_analysis/trimmed_trimmomatic/"
outdir="/data/users/${USER}/data_analysis/assembly-test/"
outspades=${outdir}"01_spades/"
outquast=${outdir}"02_quast/"

mkdir -p ${outdir}"02_quast/"
# Remove modules in your env
modules purge
# Load module
module load quast
# Use the subset of our data
cd /data/users/${USER}/
for SAMPLE in $(cat 02_tb_IDs)
do
  echo "[TB|SPADES] $SAMPLE"
  quast.py ${outspades}tb/${SAMPLE}/contigs.fasta \
    #-r /data/TB_H37Rv.fasta -g /data/TB_H37rv.gff \
    -o ${outquast}tb/${SAMPLE} \
    --threads 4 --min-contig 200 \
    --labels "MTB_Assembly"
done
```

### Quast Submission script for all samples

```bash
# Assess all assembly results for all TB ands VC samples
# In case you didn't create these above
mkdir -p /users/${USER}/scripts /users/${USER}/logs

## Now let's write our submission script
nano /users/${USER}/scripts/assembly_02.sh

#!/bin/bash
#SBATCH --job-name='quast-${USER}'
#SBATCH --nodes=1 --ntasks=16
#SBATCH --partition=Main
#SBATCH --mem=128GB
#SBATCH --output=/users/${USER}/logs/assembly_02-stdout.txt
#SBATCH --error=/users/${USER}/logs/assembly_02-stdout.txt
#SBATCH --time=12:00:00

# Define dir
outdir="/data/users/${USER}/data_analysis/assembly/"
outspades=${outdir}"01_spades/"
outquast=${outdir}"02_quast/"

mkdir -p ${outdir}"01_spades" ${outdir}"02_quast"
# Run QUAST for a subset of M. tuberculosis dataset
echo "Running QUAST analysis for M. tuberculosis..."
# Remove modules in your env
modules purge
# Load module
module load quast
# Use the sample IDs for our data
cd /data/users/${USER}/data_analysis
for SAMPLE in $(cat tb_IDs)
do
  echo "[TB|QUAST] $SAMPLE"
  quast.py ${outspades}tb/${SAMPLE}/contigs.fasta \
    #-r /data/TB_H37Rv.fasta -g /data/TB_H37rv.gff \
    -o ${outquast}tb/${SAMPLE} \
    --threads 4 --min-contig 200 \
    --labels "MTB_Assembly"
done

for SAMPLE in $(cat vc_IDs); do
  echo "[VC|QUAST] $SAMPLE"
  quast.py ${outspades}vc/${SAMPLE}/contigs.fasta \
    #-r /data/TB_H37Rv.fasta -g /data/TB_H37rv.gff \
    -o ${outquast}tb/${SAMPLE} \
    --threads 4 --min-contig 200 \
    --labels "MTB_Assembly"
done
## SAVE THIS IN NANO BY
# Ctrl + O

# CLOSE THE FILE
# ctrl + x
```

```bash
# Assess all assembly results for all TB ands VC samples
# In case you didn't create these above
mkdir -p /users/${USER}/scripts /users/${USER}/logs

## Now let's write our submission script
nano /users/${USER}/scripts/assembly_02.sh

#!/bin/bash
#SBATCH --job-name='quast-${USER}'
#SBATCH --nodes=1 --ntasks=16
#SBATCH --partition=Main
#SBATCH --mem=128GB
#SBATCH --output=/users/${USER}/logs/assembly_02-stdout.txt
#SBATCH --error=/users/${USER}/logs/assembly_02-stdout.txt
#SBATCH --time=12:00:00

# Define dir
outdir="/data/users/${USER}/data_analysis/assembly/"
outspades=${outdir}"01_spades/"
outquast=${outdir}"02_quast/"

mkdir -p ${outdir}"02_quast/tb" ${outdir}"02_quast/vc"
# Run Prokka for M. tuberculosis
echo "Starting M. tuberculosis annotation at $(date)"
for SAMPLE in $(cat tb_IDs); do
  echo "[TB|SPADES] $SAMPLE"
  prokka ${outspades}tb/${SAMPLE}/contigs.fasta \
       --outdir ${outquast}/${SAMPLE} \
       --cpus 8 --genus Mycobacterium
done
echo "M. tuberculosis annotation completed at $(date)"

for SAMPLE in $(cat vc_IDs); do
  echo "[VC|SPADES] $SAMPLE"
  spades.py -1 ${indir}vc/${SAMPLE}_1.fastq.gz \
    -2 ${indir}vc/${SAMPLE}_2.fastq.gz \
    --careful -t 8 \
    -o ${outspades}vc/${SAMPLE}
done
## SAVE THIS IN NANO BY
# Ctrl + O

# CLOSE THE FILE
# ctrl + x
```
