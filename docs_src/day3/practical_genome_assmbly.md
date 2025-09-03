# De novo Genome Assembly Practical

Before setting up you need to know the current workig directory 
```bash
# Check current working directory
pwd

# List the contents of the working diresctory
ls

# Create relevant output directories. -p so that it creates parent dir if it doesn't exist
mkdir -p /users/${USER}/results/02_assembly/tb
mkdir -p /users/${USER}/results/02_assembly/vc
mkdir -p /users/${USER}/results/03_quality/tb
mkdir -p /users/${USER}/results/03_quality/vc
```

Delete /users/${USER}/results/02_assembly and /users/user24/results/03_quality and create thenm using one line code and not four as above

```bash
mkdir /users/${USER}/results/02_assembly/tb /users/user24/results/02_assembly/vc \
  /users/${USER}/results/03_quality/tb /users/user24/results/03_quality/vc

# Request for an interactive node
srun --cpus-per-tasks=32 --mem=128g --time 5:00:00 --job-name "ephie" --pty /bin/bash

# # Clear all modules
module purge
 
# Load the required modules
module load spades/4.2.0

## How can I get the help documentation?
spades.py --help

# Lets run spades with a test dataset
/software/bio/spades/4.2.0/bin/spades.py  --test  --careful

# Run spades with test data but with more details on the reads
spades.py  -1 /software/bio/spades/4.2.0/share/spades/test_dataset/ecoli_1K_1.fq.gz \
	-2 /software/bio/spades/4.2.0/share/spades/test_dataset/ecoli_1K_2.fq.gz \
  --careful -o results/02_assembly-test

## Subset the data
head -5 tb_IDs > 02_tb_IDs
head -5 vc_IDs > 02_vc_IDs

# Perform denovo assembly for all TB samples
mkdir -p /users/${USER}/scripts
nano /users/${USER}/scripts/02_assembly.sh

#!/bin/bash
#SBATCH --job-name='02_assembly'
#SBATCH --nodes=1 --ntasks=16
#SBATCH --partition=Main
#SBATCH --mem=120GB
#SBATCH --output=/users/${USER}/logs/02_assembly-stdout.txt
#SBATCH --error=/users/${USER}/logs/02_assembly-stdout.txt
#SBATCH --time=12:00:00

for SAMPLE in $(cat tb_IDs); do
  echo "[TB|SPADES] $SAMPLE"
  spades.py -1 /users/${USER}/results/trimmed_trimmomatic/tb/${SAMPLE}_1.fastq.gz \
    -2 /users/${USER}/results/trimmed_trimmomatic/tb/${SAMPLE}_2.fastq.gz \
    --careful \
    -o /users/${USER}/results/02_assembly/tb/${SAMPLE}
done

## Save
# Ctrl + O

# CLOSE
# ctrl + x

sbatch /users/${USER}/scripts/02_assembly.sh
```

## Notes
- -1 <filename>               file with forward paired-end reads
- -2 <filename>               file with reverse paired-end reads
- --careful                   will reduce number of mismatches and short indels
- -o                          path to output file (directory shouldn't be existing)


# Quality Assessment
We need assess the quality of the contigs so that we are confident that they are of good quality.
