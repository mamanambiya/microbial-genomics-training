# Unix Commands for Pathogen Genomics - Practical Tutorial

*Adapted from [Microbial-Genomics practice scripts](https://github.com/Arash-Iranzadeh/Microbial-Genomics/blob/main/scripts/practice_unix_commands.sh)*

## Tutorial Overview

### What You'll Learn
This hands-on tutorial will teach you essential Unix commands for pathogen genomics analysis. By the end, you'll be able to:
- Navigate and organize genomics project directories
- Process FASTQ and FASTA files efficiently
- Extract meaningful information from sequencing data
- Build simple analysis pipelines
- Prepare data for HPC analysis

### Prerequisites
- Basic terminal/command line access
- No prior Unix experience required
- Access to the training server

### Time Required
- **Complete tutorial**: 2-3 hours
- **Quick essentials**: 45 minutes

### Learning Path
1. **Setup** → 2. **Basic Navigation** → 3. **File Operations** → 4. **Data Processing** → 5. **Pipeline Building**

## Setup Instructions

Before starting the Unix command exercises, you need to prepare your workspace with sample genomics data. Here's how:

### Step 1: Create Your Practice Directory

```bash
# Create a new directory for your practice exercises
# mkdir = "make directory"
# -p = create parent directories if needed (won't error if directory exists)
# ~ = shortcut for your home directory (/home/username)
mkdir -p ~/hpc_practice

# Navigate into your new directory
# cd = "change directory"
cd ~/hpc_practice
```

### Step 2: Copy Sample Genomics Data

```bash
# Copy all sample files from the shared course folder to your current location
# cp = "copy" command
# -r = "recursive" - copy directories and all their contents
# * = wildcard that matches all files
# . = current directory (destination)
cp -r /cbio/training/courses/2025/micmet-genomics/sample-data/* .
```

### Step 3: Verify Your Files

```bash
# List all files with details to confirm the copy was successful
# ls = "list" command
# -l = long format (shows permissions, size, dates)
# -a = show all files (including hidden files starting with .)
ls -la
```

### Understanding Your Sample Files

The following files are now in your directory for practice:

| File | Description | Used For |
|------|-------------|----------|
| `sample.fastq.gz` | Compressed DNA sequencing reads | Learning `zcat`, `gunzip`, file compression |
| `sample1.fastq` | Uncompressed sequencing data (3 reads) | Practicing `grep`, `wc`, sequence counting |
| `sample2.fastq` | Another FASTQ file (2 reads) | Array processing, file comparisons |
| `reference.fasta` | Reference genome sequences | Learning `grep` with FASTA headers, sequence extraction |
| `data.txt` | Tab-delimited sample metadata | Practicing `awk`, `sed`, `cut`, `sort` commands |

**File Formats Explained:**
- **FASTQ**: Contains sequences + quality scores (4 lines per read)
- **FASTA**: Contains sequences only (2 lines per sequence: header + sequence)
- **GZ**: Gzip compressed file (saves space, common in genomics)

---

## Module 1: Directory Organization for Genomics Projects

### Learning Objectives
✓ Create organized project directories
✓ Navigate between directories efficiently
✓ Understand genomics project structure

### Tutorial 1.1: Creating Your First Project Structure

#### Step 1: Start with a Simple Structure

```bash
# Create your main project directory
mkdir my_first_project

# Enter the directory
cd my_first_project

# Check where you are
pwd
# Output: /home/username/hpc_practice/my_first_project
```

#### Step 2: Add Subdirectories

```bash
# Create data directories
mkdir data
mkdir results
mkdir scripts

# List what you created
ls
# Output: data  results  scripts
```

#### Step 3: Create a Complex Structure

```bash
# Use -p to create nested directories
mkdir -p data/{raw_reads,reference_genomes,metadata}
mkdir -p results/{qc,alignment,variants,phylogeny}
mkdir -p scripts logs tmp

# View the structure
ls -la
# The -la flags show: l=long format, a=all files
```

**Try It Yourself:**
```bash
# Exercise: Create this structure
# project/
#   ├── input/
#   │   ├── sequences/
#   │   └── references/
#   └── output/
#       ├── aligned/
#       └── reports/

# Solution:
mkdir -p project/{input/{sequences,references},output/{aligned,reports}}
```

**Real-world application:**
```bash
# Set up M. tuberculosis outbreak analysis
mkdir -p mtb_outbreak_2025/{data,results,scripts,logs}
cd mtb_outbreak_2025
mkdir -p data/{fastq,references,clinical_metadata}
mkdir -p results/{fastqc,trimming,bwa_alignment,vcf_files,phylogenetic_tree}
```

---

## Module 2: File Management for Sequencing Data

### Learning Objectives
✓ Create and edit text files
✓ Copy and rename sequencing files
✓ Organize data systematically

### Tutorial 2.1: Working with Sample Lists

#### Step 1: Create a Sample List

```bash
# First, ensure you're in the right place
pwd
cd ~/hpc_practice

# Create an empty file
touch sample_list.txt

# Check it was created
ls -la sample_list.txt
# Output: -rw-r--r-- 1 user group 0 Sep 2 10:00 sample_list.txt
```

#### Step 2: Add Content to the File

```bash
# Method 1: Using echo (for single lines)
echo "Sample_001" > sample_list.txt
echo "Sample_002" >> sample_list.txt  # >> appends, > overwrites!

# Method 2: Using cat with heredoc (for multiple lines)
cat > sample_list.txt << EOF
MTB_sample_001
MTB_sample_002
MTB_sample_003
MTB_sample_004
EOF

# View what you created
cat sample_list.txt
```

#### Step 3: Count and Verify

```bash
# Count lines in file
wc -l sample_list.txt
# Output: 4 sample_list.txt

# Count words
wc -w sample_list.txt
# Output: 4 sample_list.txt

# Get full statistics
wc sample_list.txt
# Output: 4  4  58 sample_list.txt
#        (lines, words, characters)
```

### Tutorial 2.2: Organizing Sequencing Files

#### Step 1: Copy Files Safely

```bash
# Copy sample files to practice with
cp sample*.fastq .

# List files before renaming
ls sample*.fastq
# Output: sample1.fastq  sample2.fastq

# Create a backup first (always!)
mkdir backups
cp sample*.fastq backups/
```

#### Step 2: Rename Files Systematically

```bash
# Rename a single file
mv sample1.fastq patient001_reads.fastq

# Batch rename using a loop
for file in sample*.fastq; do
    # Extract the number from filename
    num=$(echo $file | grep -o '[0-9]\+')
    # Create new name
    newname="patient_${num}_sequences.fastq"
    echo "Renaming $file to $newname"
    mv "$file" "$newname"
done

# Verify the renaming
ls *.fastq
```

#### Practice Exercise:
```bash
# Exercise: Create copies with dates
# Copy sample.fastq.gz to sample_20250902.fastq.gz

# Solution:
date_stamp=$(date +%Y%m%d)
cp sample.fastq.gz sample_${date_stamp}.fastq.gz
```

---

## Module 3: Viewing and Inspecting Genomics Files

### Learning Objectives
✓ View compressed and uncompressed files
✓ Count sequences in FASTQ/FASTA files
✓ Extract specific parts of files

### Tutorial 3.1: Working with FASTQ Files

#### Step 1: View Compressed Files

```bash
# View first 4 lines (1 complete read) of compressed file
zcat sample.fastq.gz | head -4
# Output:
# @SEQ_001
# ACGTACGTACGTACGTACGTACGTACGTACGTACGT
# +
# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# View first 2 reads (8 lines)
zcat sample.fastq.gz | head -8
```

#### Step 2: Count Sequences

```bash
# Count total lines
zcat sample.fastq.gz | wc -l
# Output: 12  (for 3 reads)

# Count number of reads (FASTQ has 4 lines per read)
zcat sample.fastq.gz | wc -l | awk '{print $1/4}'
# Output: 3

# Alternative: Count sequence headers
zcat sample.fastq.gz | grep -c "^@"
# Output: 3
```

#### Step 3: View with Line Numbers

```bash
# Add line numbers to output
zcat sample.fastq.gz | head -8 | cat -n
# Output:
#      1	@SEQ_001
#      2	ACGTACGTACGTACGTACGTACGTACGTACGTACGT
#      3	+
#      4	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
#      5	@SEQ_002
#      6	TGCATGCATGCATGCATGCATGCATGCATGCATGCA
#      7	+
#      8	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

### Tutorial 3.2: Working with FASTA Files

#### Step 1: View FASTA Headers

```bash
# View first line (header) of FASTA
head -1 reference.fasta
# Output: >Sequence_1 gene=ABC123

# View all headers in file
grep "^>" reference.fasta
# Output:
# >Sequence_1 gene=ABC123
# >Sequence_2 gene=DEF456

# Count sequences
grep -c "^>" reference.fasta
# Output: 2
```

#### Step 2: Extract Sequences

```bash
# View first 2 lines (header + sequence start)
head -2 reference.fasta

# Skip header, view sequence only
tail -n +2 reference.fasta | head -1
# Output: ACGTACGTACGTACGTACGTACGTACGTACGTACGT

# Get sequence length
tail -n +2 reference.fasta | head -1 | wc -c
# Output: 37 (includes newline)
```

#### Practice Exercise:
```bash
# Exercise: Count total bases in all sequences
# Hint: Remove headers first

# Solution:
grep -v "^>" reference.fasta | tr -d '\n' | wc -c
```

---

## Module 4: Searching and Filtering Genomics Data

### Learning Objectives
✓ Search for patterns in files
✓ Filter genomics data
✓ Use regular expressions

### Tutorial 4.1: Basic Pattern Searching

#### Step 1: Simple Searches

```bash
# Search for a word in a file
grep "resistant" data.txt
# Output: Sample1 100 resistant
#         Sample3 150 resistant

# Case-insensitive search (-i flag)
grep -i "SAMPLE" data.txt
# Finds: Sample1, Sample2, etc.

# Count matches (-c flag)
grep -c "resistant" data.txt
# Output: 2

# Show line numbers (-n flag)
grep -n "resistant" data.txt
# Output: 1:Sample1 100 resistant
#         3:Sample3 150 resistant
```

#### Step 2: Search in Genomics Files

```bash
# Find sequence headers in FASTA
grep "^>" reference.fasta
# ^ means "start of line"

# Find adapter sequences in FASTQ
grep "AGATCGGAAGAG" sample1.fastq

# Search compressed files
zcat sample.fastq.gz | grep "ACGT"
```

### Tutorial 4.2: Advanced Pattern Matching

#### Using Regular Expressions

```bash
# Find lines with numbers
grep '[0-9]' data.txt
# [0-9] matches any digit

# Find lines ending with specific pattern
grep 'resistant$' data.txt
# $ means "end of line"

# Extract only the matching part (-o flag)
echo "Sample123" | grep -o '[0-9]\+'
# Output: 123
# \+ means "one or more"

# Multiple patterns (OR)
grep -E 'resistant|sensitive' data.txt
# -E enables extended regex
```

#### Practice Exercise:
```bash
# Exercise: Find all samples with values > 150
# Hint: Use awk instead of grep for numeric comparisons

# Solution:
awk '$2 > 150' data.txt
# Output: Sample2 200 sensitive
#         Sample4 300 sensitive
```

---

## Module 5: Text Processing for Genomics

### Learning Objectives
✓ Extract specific columns from files
✓ Perform calculations on data
✓ Transform text efficiently

### Tutorial 5.1: Using awk for Data Processing

#### Step 1: Extract Columns

```bash
# Print specific columns (1st and 2nd)
awk '{print $1, $2}' data.txt
# Output: Sample1 100
#         Sample2 200
#         Sample3 150
#         Sample4 300

# Print with custom separator
awk '{print $1 "," $2}' data.txt
# Output: Sample1,100

# Add text to output
awk '{print "Sample:" $1 " Value:" $2}' data.txt
```

#### Step 2: Perform Calculations

```bash
# Sum values in column 2
awk '{sum+=$2} END {print "Total:", sum}' data.txt
# Output: Total: 750

# Calculate average
awk '{sum+=$2; count++} END {print "Average:", sum/count}' data.txt
# Output: Average: 187.5

# Filter based on value
awk '$2 > 150 {print $0}' data.txt
# Output: Sample2 200 sensitive
#         Sample4 300 sensitive
```

### Tutorial 5.2: Using sed for Text Manipulation

#### Step 1: Basic Substitutions

```bash
# Replace text (s = substitute)
echo "Sample1" | sed 's/1/A/'
# Output: SampleA

# Global replacement (g = global)
echo "Sample111" | sed 's/1/A/g'
# Output: SampleAAA

# Replace in file and save
sed 's/Sample/Patient/g' data.txt > patients.txt

# Edit file in place (careful!)
sed -i.bak 's/Sample/Patient/g' data.txt
# Creates data.txt.bak as backup
```

#### Step 2: Advanced Manipulations

```bash
# Delete lines containing pattern
sed '/sensitive/d' data.txt
# Removes lines with "sensitive"

# Add text to beginning of lines
sed 's/^/PREFIX_/' data.txt
# Adds PREFIX_ to each line start

# Convert spaces to tabs
sed 's/ /\t/g' data.txt > data.tsv
```

---

## 6. File Permissions and Management

### Managing Permissions

```bash
# Make scripts executable
chmod +x scripts/analysis_pipeline.sh

# Protect raw data from accidental modification
chmod 444 data/raw_reads/*.fastq.gz

# Set directory permissions
chmod 755 results/

# Check permissions
ls -la data/raw_reads/
```

---

## 7. Sorting and Unique Operations

### Processing Sample Lists

```bash
# Sort sample names
sort data/sample_list.txt

# Sort numerically by coverage
sort -k2 -n coverage_stats.txt

# Get unique mutations
cut -f1,2 variants.txt | sort | uniq

# Count occurrences of each mutation
cut -f3 mutations.txt | sort | uniq -c | sort -rn
```

---

## 8. Pipelines and Redirection

### Creating Simple Analysis Pipelines

```bash
# Count reads per sample and save to file
for file in data/raw_reads/*.fastq.gz; do
    sample=$(basename $file .fastq.gz)
    count=$(zcat $file | wc -l | awk '{print $1/4}')
    echo -e "$sample\t$count"
done > results/qc/read_counts.txt

# Extract and count specific genes
grep "^>" reference.fasta | cut -d' ' -f1 | sed 's/>//' | sort | uniq -c > gene_counts.txt

# Process multiple VCF files
for vcf in results/variants/*.vcf; do
    sample=$(basename $vcf .vcf)
    pass_count=$(grep -c "PASS" $vcf)
    total_count=$(grep -v "^#" $vcf | wc -l)
    echo -e "$sample\t$total_count\t$pass_count"
done > results/variant_summary.txt
```

---

## 9. Practical SLURM Integration

### Preparing Files for HPC Analysis

```bash
#!/bin/bash
#SBATCH --job-name=prep_pathogen_data
#SBATCH --time=00:30:00
#SBATCH --mem=4GB

# Create directory structure
mkdir -p ${SLURM_JOB_ID}_analysis/{data,results,logs}

# Copy and organize files
cp /shared/data/*.fastq.gz ${SLURM_JOB_ID}_analysis/data/

# Generate file list for processing
ls ${SLURM_JOB_ID}_analysis/data/*.fastq.gz > file_list.txt

# Count and verify files
echo "Total files to process: $(wc -l < file_list.txt)"

# Create metadata file
for file in ${SLURM_JOB_ID}_analysis/data/*.fastq.gz; do
    size=$(du -h $file | cut -f1)
    reads=$(zcat $file | wc -l | awk '{print $1/4}')
    echo "$(basename $file)\t$size\t$reads"
done > ${SLURM_JOB_ID}_analysis/data/file_metadata.tsv

echo "Preparation complete. Ready for analysis."
```

---

## Hands-On Exercise: Complete Pathogen Analysis Workflow

### Exercise Overview
Build a complete analysis pipeline step-by-step, applying all the skills you've learned.

### Part 1: Setup and Data Preparation

```bash
# Step 1: Create project structure
mkdir -p pathogen_practice/{data,results,scripts}
cd pathogen_practice
pwd  # Verify you're in the right place

# Step 2: Create sample metadata
cat > data/samples.txt << EOF
Mtb_patient_001_resistant
Mtb_patient_002_susceptible  
Mtb_patient_003_resistant
Salmonella_outbreak_001
Salmonella_outbreak_002
EOF

# Verify the file was created
cat data/samples.txt

```

### Part 2: Data Analysis Tasks

#### Task 1: Find Resistant Samples
```bash
# Use grep to find resistant samples
grep "resistant" data/samples.txt

# Save results to file
grep "resistant" data/samples.txt > results/resistant_samples.txt

# Count how many
grep -c "resistant" data/samples.txt

```

#### Task 2: Count by Pathogen Type
```bash
# Count MTB samples
grep -c "Mtb" data/samples.txt

# Count Salmonella samples
grep -c "Salmonella" data/samples.txt

# Save counts
echo "MTB samples: $(grep -c 'Mtb' data/samples.txt)" > results/pathogen_counts.txt
echo "Salmonella samples: $(grep -c 'Salmonella' data/samples.txt)" >> results/pathogen_counts.txt

```

#### Task 3: Generate Summary Report
```bash
# Create a comprehensive summary
cat > results/summary_report.txt << EOF
=== Pathogen Analysis Summary ===
Date: $(date)
Total samples: $(wc -l < data/samples.txt)
Resistant samples: $(grep -c "resistant" data/samples.txt)
Susceptible samples: $(grep -c "susceptible" data/samples.txt)
MTB samples: $(grep -c "Mtb" data/samples.txt)
Salmonella samples: $(grep -c "Salmonella" data/samples.txt)
=================================
EOF

# Display the report
cat results/summary_report.txt
```

### Part 3: Challenge Exercises

#### Challenge 1: Extract Sample IDs
```bash
# Extract just the patient IDs (hint: use cut or awk)
# Try it yourself first!

# Solution:
cut -d'_' -f2,3 data/samples.txt
# Or using awk:
awk -F'_' '{print $2"_"$3}' data/samples.txt
```

#### Challenge 2: Sort and Count Unique Pathogens
```bash
# Extract pathogen names and count occurrences
# Try it yourself first!

# Solution:
cut -d'_' -f1 data/samples.txt | sort | uniq -c
```

#### Challenge 3: Create a Pipeline
```bash
# Find all resistant MTB samples in one command
# Try it yourself first!

# Solution:
grep "Mtb" data/samples.txt | grep "resistant"
# Or more elegantly:
grep "Mtb.*resistant" data/samples.txt
```

---

## Tips for Pathogen Genomics Unix Usage

1. **Always work with copies** of raw sequencing data
2. **Use meaningful file names** with sample IDs and dates
3. **Document your commands** in scripts for reproducibility
4. **Check file integrity** after transfers (md5sum)
5. **Compress large files** to save space (gzip/bgzip)
6. **Use screen or tmux** for long-running processes
7. **Regular backups** of analysis results
8. **Version control** for scripts (git)

---

## Common File Formats in Pathogen Genomics

| Extension | Format | View Command | Description |
|-----------|--------|--------------|-------------|
| `.fastq.gz` | Compressed FASTQ | `zcat file.fastq.gz \| head` | Raw sequencing reads |
| `.fasta` | FASTA | `cat file.fasta` | Reference genomes |
| `.sam/.bam` | SAM/BAM | `samtools view file.bam \| head` | Alignments |
| `.vcf` | VCF | `cat file.vcf` | Variant calls |
| `.gff/.gtf` | GFF/GTF | `cat file.gff` | Gene annotations |
| `.newick/.tree` | Newick | `cat file.tree` | Phylogenetic trees |

---

## Troubleshooting Guide

### Common Issues and Solutions

#### Issue 1: "Permission denied"
```bash
# Problem: Can't access or modify a file
# Solution: Check permissions
ls -la filename

# Fix: Change permissions if you own the file
chmod u+rw filename
```

#### Issue 2: "No such file or directory"
```bash
# Problem: File path is wrong
# Solution: Check your current location
pwd

# List files to verify
ls -la

# Use absolute paths to be sure
/full/path/to/file
```

#### Issue 3: "Command not found"
```bash
# Problem: Tool not installed or not in PATH
# Solution: Check if command exists
which command_name

# Load module if available
module avail
module load tool_name
```

#### Issue 4: File is empty or corrupted
```bash
# Check file size
ls -lh filename

# Check file type
file filename

# For compressed files, test integrity
gzip -t file.gz
```

#### Issue 5: Out of disk space
```bash
# Check available space
df -h

# Find large files
du -sh * | sort -h

# Clean up temporary files
rm -rf tmp/*
```

### Best Practices to Avoid Issues

1. **Always backup before modifying**
   ```bash
   cp important_file important_file.backup
   ```

2. **Use tab completion** to avoid typos
   ```bash
   cat sam[TAB]  # Completes to sample.fastq
   ```

3. **Preview commands with echo first**
   ```bash
   echo mv *.fastq backup/  # See what would happen
   mv *.fastq backup/       # Then run for real
   ```

4. **Check file contents before processing**
   ```bash
   head -5 file.txt  # Preview first 5 lines
   wc -l file.txt    # Count total lines
   ```

---

## Quick Reference Card

### Essential Commands Summary

| Task | Command | Example |
|------|---------|---------|
| List files | `ls -la` | `ls -la *.fastq` |
| Change directory | `cd` | `cd ~/hpc_practice` |
| Create directory | `mkdir -p` | `mkdir -p data/reads` |
| Copy files | `cp -r` | `cp sample.fastq backup/` |
| Move/rename | `mv` | `mv old.txt new.txt` |
| View compressed | `zcat` | `zcat file.gz \| head` |
| Count lines | `wc -l` | `wc -l sample.txt` |
| Search text | `grep` | `grep "pattern" file` |
| Extract columns | `awk` | `awk '{print $1}' file` |
| Replace text | `sed` | `sed 's/old/new/g' file` |
| Sort data | `sort` | `sort -n numbers.txt` |
| Get unique | `uniq` | `sort file \| uniq` |

---

## Next Steps

After mastering these Unix commands, you're ready to:
1. **Submit SLURM jobs** - See [High Performance Computing with SLURM: Practical Tutorial](./slurm-practical-tutorial.md)
2. **Learn HPC concepts** - See [HPC and ILIFU Training Materials](./hpc-ilifu-training.md)
3. **Build analysis pipelines** with Nextflow
4. **Perform quality control** with FastQC
5. **Align reads** with BWA
6. **Call variants** with SAMtools/BCFtools

Remember: Unix commands are the foundation of all bioinformatics pipelines!