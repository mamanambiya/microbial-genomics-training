# De novo Genome Assembly Practical

Before setting up you need to know the current workig directory 
```bash
# Check current working directory
pwd

# List the contents of the working diresctory
ls

# Create relevant output directories. -p so that it creates parent dir if it doesn't exist
mkdir -p /users/user24/results/02_assembly/tb
mkdir -p /users/user24/results/02_assembly/vc
mkdir -p /users/user24/results/03_quality/tb
mkdir -p /users/user24/results/03_quality/vc
```

Delete /users/user24/results/02_assembly and /users/user24/results/03_quality and create thenm using one line code and not four as above

```bash
mkdir /users/user24/results/02_assembly/tb /users/user24/results/02_assembly/vc \
  /users/user24/results/03_quality/tb /users/user24/results/03_quality/vc

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

# Perform denovo assembly
for SAMPLE in $(cat 02_tb_IDs); do
  echo "[TB|SPADES] $SAMPLE"
  spades.py -1 /users/${USER}/results/trimmed_trimmomatic/tb/${SAMPLE}_1.fastq.gz \
    -2 /users/${USER}/results/trimmed_trimmomatic/tb/${SAMPLE}_2.fastq.gz \
    --careful \
    -o /users/${USER}/results/02_assembly/tb/${SAMPLE}
done
```

## Notes
- -1 <filename>               file with forward paired-end reads
- -2 <filename>               file with reverse paired-end reads
- --careful                   will reduce number of mismatches and short indels
- -o                          path to output file (directory shouldn't be existing)
