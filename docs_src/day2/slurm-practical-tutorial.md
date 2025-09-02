# SLURM Job Examples and Templates

## Getting Started - Setup Instructions

Before starting the exercises, you need to set up your working environment and copy the sample data files to your home directory. Follow these steps:

### Step 1: Create Your Working Directory

```bash
# The mkdir command creates a new directory
# The -p flag creates parent directories if they don't exist
mkdir -p ~/hpc_practice

# Change to your new working directory
# The ~ symbol represents your home directory
cd ~/hpc_practice
```

### Step 2: Copy Sample Data Files

```bash
# Copy all sample data from the shared course directory to your current directory
# The -r flag means "recursive" - it copies directories and their contents
# The * wildcard matches all files in the source directory
# The . (dot) means "current directory" (where you are now)
cp -r /cbio/training/courses/2025/micmet-genomics/sample-data/* .
```

### Step 3: Verify Your Setup

```bash
# List all files in your directory to confirm they copied correctly
# The -l flag shows detailed information (permissions, size, date)
# The -a flag shows all files including hidden ones (starting with .)
ls -la

# You should see these files:
# - sample.fastq.gz    : Compressed DNA sequencing data (gzipped FASTQ format)
# - sample1.fastq      : Uncompressed sequencing reads for practice
# - sample2.fastq      : Another set of sequencing reads
# - reference.fasta    : Reference genome sequence for alignment exercises
# - data.txt          : Tab-delimited data for text processing examples
```

### What These Files Contain

- **FASTQ files**: Contain DNA sequences and quality scores from sequencing machines
- **FASTA files**: Contain reference sequences without quality scores
- **Text files**: Contain structured data for analysis practice

Now you're ready to start the exercises!

## Table of Contents

1. [Prerequisites - Unix Commands](#prerequisites)
2. [Getting Started - Your First Job](#getting-started)
3. [Basic Job Templates](#basic-templates)
4. [Python and Bash Examples](#python-and-bash-examples)
5. [Advanced Job Types](#advanced-jobs)
6. [Resource Optimization](#resource-optimization)
7. [Troubleshooting Examples](#troubleshooting)

---

## Prerequisites

### Essential Unix Commands for HPC

Before submitting SLURM jobs, master these Unix commands for pathogen genomics:

```bash
# Navigate and organize
mkdir -p project/{data,results,scripts}
cd project
pwd

# Inspect FASTQ files
zcat sample.fastq.gz | head -20
zcat sample.fastq.gz | wc -l | awk '{print $1/4}'  # Count reads

# Search and filter
grep "^>" reference.fasta  # Find FASTA headers
grep -c "PASS" variants.vcf  # Count PASS variants

# Process text
awk '{print $1, $2}' data.txt
sed 's/old/new/g' file.txt
```

📚 **Full Unix guide:** See `unix-commands-pathogen-examples.md` for comprehensive examples and exercises.

---

## Tutorial: Your First SLURM Jobs - Step by Step

### Tutorial Overview

In this hands-on tutorial, you'll learn to:
1. Write and submit your first SLURM job
2. Monitor job status and view outputs
3. Run Python scripts on HPC
4. Process genomics data with SLURM
5. Handle errors and optimize resources

**Time needed:** 30-45 minutes
**Prerequisites:** Basic Unix commands (covered above)

---

### Tutorial 1: Hello World on HPC

#### Step 1: Write Your First Job Script

Create a simple SLURM job that prints a greeting:

```bash
#!/bin/bash
#SBATCH --job-name=hello
#SBATCH --time=00:05:00

echo "Hello from HPC!"
echo "This job ran on node: $(hostname)"
echo "Current time: $(date)"
```

#### Step 2: Save the Script

```bash
# Use nano editor to create the file
nano hello.sh

# Paste the script above, then:
# Press Ctrl+X to exit
# Press Y to save
# Press Enter to confirm filename
```

#### Step 3: Submit Your Job

```bash
# Submit the job to SLURM
sbatch hello.sh
```

You'll see: `Submitted batch job 12345` (your job ID will differ)

#### Step 4: Monitor Your Job

```bash
# Check if your job is running
squeue -u $USER

# You'll see something like:
# JOBID PARTITION     NAME     USER ST       TIME  NODES
# 12345      Main    hello  yourname  R       0:01      1
# ST column: PD=Pending, R=Running, CG=Completing
```

#### Step 5: View the Output

```bash
# Once job completes (status disappears from squeue)
# View the output file (replace 12345 with your job ID)
cat slurm-12345.out
```

**Example run:**

```bash
$ sbatch hello.sh
Submitted batch job 10

$ squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                10  training    hello   mamana  R       0:01      1 compute-1-sep2025

$ cat slurm-10.out
Hello from HPC!
This job ran on node: compute-1-sep2025
Current time: Mon Sep 1 23:57:07 SAST 2025
```

---

### Tutorial 2: Running Python on HPC

#### Step 1: Create a Python Job Script

Let's run Python code on the cluster:

```bash
#!/bin/bash
#SBATCH --job-name=python_hello
#SBATCH --time=00:10:00
#SBATCH --mem=1GB

# Load Python (or use system python3 if modules not available)
module load python/3.12.3  # Or use system python3 || echo "Using system Python"

# Run your Python script
python3 << 'EOF'
print("Hello from Python on HPC!")
import os
print(f"Running on: {os.uname().nodename}")

# Simple calculation
result = sum(range(1000))
print(f"Sum of 0-999 = {result}")
EOF
```

#### Step 2: Submit and Monitor

```bash
# Save the script
nano python_job.sh
# (paste script, save with Ctrl+X, Y, Enter)

# Submit the job
sbatch python_job.sh

# Watch it run (updates every 2 seconds)
watch -n 2 squeue -u $USER
# Press Ctrl+C to stop watching
```

#### Step 3: Check the Output

```bash
# Find your output file
ls -lt slurm-*.out | head -5

# View the results
cat slurm-[JOBID].out
```

**Expected output:**
```text
Hello from Python on HPC!
Running on: compute-1-sep2025
Sum of 0-999 = 499500
```

#### Common Issues and Solutions:

| Problem | Solution |
|---------|----------|
| "Module not found" | Use `python3` instead of loading module |
| "Python: command not found" | Check with `which python3` |
| Job stays pending too long | Check resources with `sinfo` |

---

### Tutorial 3: Real Genomics Analysis

#### Objective
Process FASTQ files using SLURM, simulating a real bioinformatics pipeline.

#### Step 1: Create the Analysis Script

This script demonstrates a typical genomics workflow:

```bash
#!/bin/bash
#SBATCH --job-name=fastq_analysis
#SBATCH --time=00:05:00
#SBATCH --mem=2GB
#SBATCH --cpus-per-task=2

echo "=== FASTQ Analysis Pipeline Starting ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

# Create sample FASTQ files for analysis
echo "Creating sample FASTQ files..."
cat > sample1.fastq << 'EOF'
@SEQ_1
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
@SEQ_2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
EOF

cat > sample2.fastq << 'EOF'
@SEQ_3
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
EOF

echo "Step 1: Initial file validation..."
sleep 45  # Simulate file checking and validation

echo "Step 2: Sequence counting and basic stats..."
for file in sample*.fastq; do
    echo "Processing $file..."
    sequences=$(wc -l < "$file")
    sequences=$((sequences / 4))
    echo "  Found $sequences sequences"
    
    # Simulate per-file analysis time
    echo "  Analyzing sequence lengths..."
    sleep 25  # Processing time per file
    
    avg_length=60
    echo "  Average sequence length: ${avg_length}bp"
done

echo "Step 3: Quality score analysis..."
echo "Analyzing quality scores across all sequences..."
sleep 60  # Simulate quality analysis

echo "Step 4: Generating contamination check..."
echo "Checking for adapter sequences and contaminants..."
sleep 45  # Simulate contamination screening

echo "Step 5: Creating final summary report..."
total_sequences=0
for file in sample*.fastq; do
    seqs=$(wc -l < "$file")
    seqs=$((seqs / 4))
    total_sequences=$((total_sequences + seqs))
done

echo "Step 6: Finalizing results and cleanup..."
sleep 20  # Final processing and cleanup

echo "=== Analysis Complete ==="
echo "Total sequences analyzed: $total_sequences"
echo "Analysis completed at: $(date)"
echo "Total runtime: ~4 minutes"

# Create a summary file
cat > analysis_summary.txt << EOF
FASTQ Analysis Summary
=====================
Total files processed: 2
Total sequences: $total_sequences
Average sequence length: 60bp
Quality check: PASSED
Contamination check: CLEAN
Analysis date: $(date)
EOF

echo "Summary report saved to: analysis_summary.txt"
```

#### Step 2: Submit and Monitor the Job

```bash
# Save the script
nano fastq_analysis.sh

# Submit the job
sbatch fastq_analysis.sh
# Note your job ID (e.g., "Submitted batch job 12347")
```

#### Step 3: Monitor Job Progress in Real-Time

Open multiple terminal windows to watch different aspects:

**Terminal 1: Submit and monitor queue**
```bash
# Submit the job
sbatch fastq_analysis.sh
Submitted batch job 15

# Watch it in the queue (run multiple times)
# Watch the queue (repeat every 10 seconds)
squeue -u $USER
# Status codes: PD=Pending, R=Running, CG=Completing
```

**Terminal 2: Watch live output**
```bash
# Once job starts running (status = R), watch the output
tail -f slurm-15.out
# Press Ctrl+C to stop watching
```

**Terminal 3: Check job details**
```bash
# Get detailed job information
scontrol show job 15
```

#### Step 4: Understanding Job States

During the 4-minute runtime, you'll observe these states:

| Time | Status | What's Happening |
|------|--------|------------------|
| 0:00-0:05 | PD (Pending) | Job waiting for resources |
| 0:05-4:00 | R (Running) | Job executing on compute node |
| 4:00+ | - (Completed) | Job finished, no longer in queue |

**Timeline of analysis steps:**
- **0:00-0:45** - File validation
- **0:45-1:35** - Sequence counting (sample1.fastq)
- **1:35-2:25** - Sequence counting (sample2.fastq)  
- **2:25-3:25** - Quality score analysis
- **3:25-4:10** - Contamination screening
- **4:10-4:30** - Final report generation

**Learning opportunity:** This 4-minute window allows everyone to:
- Practice using `squeue` to monitor jobs multiple times
- See job state transitions and timing in real-time  
- Understand queue system behavior with sufficient time for discussion
- Watch live output with `tail -f` to see analysis progress
- Check intermediate results and final efficiency reports

> **💡 Training Tip:** Have participants submit this job, then use the 4-minute window to demonstrate:
> - Refreshing `squeue -u $USER` every 30 seconds to track progress
> - Using `scontrol show job JOBID` for detailed job information
> - Explaining what PENDING vs RUNNING states mean
> - Demonstrating `tail -f slurm-JOBID.out` to watch live step-by-step output
> - Discussing resource allocation while job runs
> - Explaining the difference between walltime and CPU time

**Expected final output files:**
- `slurm-JOBID.out` - Complete log of all analysis steps
- `analysis_summary.txt` - Final summary report
- `sample1.fastq` & `sample2.fastq` - Generated test data files

**Sample log output:**
```text
=== FASTQ Analysis Pipeline Starting ===
Job ID: 15
Node: compute-2-sep2025
Start time: Mon Sep 2 10:15:23 SAST 2025
Creating sample FASTQ files...
Step 1: Initial file validation...
Step 2: Sequence counting and basic stats...
Processing sample1.fastq...
  Found 2 sequences
  Analyzing sequence lengths...
  Average sequence length: 60bp
Processing sample2.fastq...
  Found 2 sequences
  Analyzing sequence lengths...
  Average sequence length: 60bp
Step 3: Quality score analysis...
Analyzing quality scores across all sequences...
Step 4: Generating contamination check...
Checking for adapter sequences and contaminants...
Step 5: Creating final summary report...
Step 6: Finalizing results and cleanup...
=== Analysis Complete ===
Total sequences analyzed: 4
Analysis completed at: Mon Sep 2 10:19:45 SAST 2025
Total runtime: ~4 minutes
Summary report saved to: analysis_summary.txt
```

---

## Practice Exercises

### Exercise 1: Modify and Submit a Job

**Task:** Modify the hello.sh script to include your name and the current date.

```bash
# Step 1: Edit the script
nano hello.sh

# Step 2: Add these lines:
echo "Submitted by: [YOUR NAME]"
echo "Analysis date: $(date +%Y-%m-%d)"

# Step 3: Submit and check
sbatch hello.sh
squeue -u $USER
```

### Exercise 2: Resource Monitoring

**Task:** Create a job that uses specific resources and monitor them.

```bash
#!/bin/bash
#SBATCH --job-name=resource_test
#SBATCH --time=00:02:00
#SBATCH --mem=500MB
#SBATCH --cpus-per-task=2

echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"
echo "Allocated Memory: $SLURM_MEM_PER_NODE MB"
echo "Running on node: $(hostname)"

# Use the allocated CPUs
stress --cpu $SLURM_CPUS_PER_TASK --timeout 30s
```

### Exercise 3: Array Jobs

**Task:** Process multiple files in parallel using array jobs.

```bash
#!/bin/bash
#SBATCH --job-name=array_demo
#SBATCH --array=1-3
#SBATCH --time=00:05:00

echo "Processing file number: $SLURM_ARRAY_TASK_ID"
# Your processing command here
```

---

## Basic Templates

### 1. Standard Job Template

```bash
#!/bin/bash
#SBATCH --job-name=my_job          # Give your job a name
#SBATCH --time=01:00:00            # Max runtime (1 hour)
#SBATCH --mem=4GB                  # Memory needed
#SBATCH --output=output_%j.log     # Output file (%j = job ID)

# Load software you need
module load python/3.12.3  # Or use system python3

# Run your command
echo "Job started on $(hostname) at $(date)"
python my_script.py
echo "Job completed at $(date)"
```

### 2. Multi-core Parallel Job

```bash
#!/bin/bash
#SBATCH --job-name=parallel_job
#SBATCH --partition=Main
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --output=parallel_%j.log

module load python/3.12.3  # Or use system python3

# Use all available cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Using $SLURM_CPUS_PER_TASK CPU cores"
python parallel_script.py
```

---

## Python and Bash Examples

### Python Jobs

#### Basic Python Analysis

```bash
#!/bin/bash
#SBATCH --job-name=python_analysis
#SBATCH --time=01:30:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4

# Load Python module
module load python/3.12.3  # Or use system python3

# Run your analysis
python genome_analysis.py sample_data.fasta
```

#### Python with Virtual Environment

```bash
#!/bin/bash
#SBATCH --job-name=python_venv
#SBATCH --time=02:00:00
#SBATCH --mem=16GB

module load python/3.12.3  # Or use system python3

# Create and activate virtual environment
python -m venv pathogen_env
source pathogen_env/bin/activate

# Install bioinformatics packages
pip install biopython pandas numpy matplotlib

# Run pathogen analysis
python pathogen_analysis.py
```

#### Pathogen Genomics - SNP Analysis

```bash
#!/bin/bash
#SBATCH --job-name=snp_analysis
#SBATCH --time=04:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

module load python/3.12.3  # Or use system python3

# Python script for SNP analysis
python << 'EOF'
import pandas as pd
from multiprocessing import Pool
import os

def analyze_sample(vcf_file):
    """Analyze SNPs in a VCF file"""
    print(f"Processing {vcf_file}")
    
    # Count SNPs (simplified example)
    with open(vcf_file, 'r') as f:
        snp_count = sum(1 for line in f if not line.startswith('#'))
    
    return vcf_file, snp_count

# Get all VCF files
vcf_files = [f for f in os.listdir('.') if f.endswith('.vcf')]

# Use all available CPU cores
with Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as pool:
    results = pool.map(analyze_sample, vcf_files)

# Save results
results_df = pd.DataFrame(results, columns=['Sample', 'SNP_Count'])
results_df.to_csv('snp_analysis_results.csv', index=False)
print(f"Analyzed {len(vcf_files)} samples")
EOF
```

### Bash/Shell Script Jobs

#### Basic FASTQ Processing

```bash
#!/bin/bash
#SBATCH --job-name=fastq_processing
#SBATCH --time=01:00:00
#SBATCH --mem=4GB

# Process multiple FASTQ files
for file in *.fastq; do
    echo "Processing $file..."
    
    # Count sequences (FASTQ has 4 lines per sequence)
    sequences=$(wc -l < "$file")
    sequences=$((sequences / 4))
    
    # Get basic stats
    echo "File: $file - Sequences: $sequences"
    
    # Count reads with quality scores above threshold
    good_reads=$(awk 'NR%4==0 && length($0)>20' "$file" | wc -l)
    echo "High quality reads: $good_reads"
done

echo "Processing complete!"
```

**Expected output:**

```text
Processing sample1.fastq...
File: sample1.fastq - Sequences: 3
High quality reads: 3
Processing sample2.fastq...
File: sample2.fastq - Sequences: 2
High quality reads: 2
Processing complete!
```

#### Pathogen Genomics Pipeline

```bash
#!/bin/bash
#SBATCH --job-name=pathogen_pipeline
#SBATCH --time=06:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16

# Load bioinformatics tools
module load fastqc/0.12.1  # Check available version with 'module avail'
# module load trimmomatic  # Install if needed
module load bwa/github  # Check available version
module load samtools/1.22.1
module load bcftools/1.22

# Sample information
SAMPLE="pathogen_sample"
REFERENCE="reference_genome.fasta"

echo "=== Pathogen Genomics Pipeline Starting ==="
echo "Sample: $SAMPLE"
echo "Reference: $REFERENCE"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# Step 1: Quality control
echo "Step 1: Running FastQC..."
mkdir -p qc_reports
fastqc "${SAMPLE}_R1.fastq" "${SAMPLE}_R2.fastq" -o qc_reports/

# Step 2: Trim low-quality reads and adapters
echo "Step 2: Trimming reads..."
trimmomatic PE -threads $SLURM_CPUS_PER_TASK \
    "${SAMPLE}_R1.fastq" "${SAMPLE}_R2.fastq" \
    "${SAMPLE}_R1_trimmed.fastq" "${SAMPLE}_R1_unpaired.fastq" \
    "${SAMPLE}_R2_trimmed.fastq" "${SAMPLE}_R2_unpaired.fastq" \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: Align reads to reference genome
echo "Step 3: Aligning to reference genome..."
bwa mem -t $SLURM_CPUS_PER_TASK "$REFERENCE" \
    "${SAMPLE}_R1_trimmed.fastq" "${SAMPLE}_R2_trimmed.fastq" | \
    samtools sort -@ $SLURM_CPUS_PER_TASK -o "${SAMPLE}_sorted.bam"

# Step 4: Index BAM file
echo "Step 4: Indexing BAM file..."
samtools index "${SAMPLE}_sorted.bam"

# Step 5: Variant calling
echo "Step 5: Calling variants..."
bcftools mpileup -f "$REFERENCE" "${SAMPLE}_sorted.bam" | \
    bcftools call -mv -Oz -o "${SAMPLE}_variants.vcf.gz"

# Step 6: Index VCF and get basic stats
echo "Step 6: Processing variants..."
bcftools index "${SAMPLE}_variants.vcf.gz"
bcftools stats "${SAMPLE}_variants.vcf.gz" > "${SAMPLE}_variant_stats.txt"

# Summary
echo "=== Pipeline Summary ==="
echo "Alignment stats:"
samtools flagstat "${SAMPLE}_sorted.bam"

echo "Variant counts:"
bcftools view -H "${SAMPLE}_variants.vcf.gz" | wc -l

echo "=== Pathogen Genomics Pipeline Complete ==="
```

#### Multi-Sample Outbreak Analysis

```bash
#!/bin/bash
#SBATCH --job-name=outbreak_analysis
#SBATCH --time=08:00:00
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=32

# Load required modules
module load python/3.12.3  # Or use system python3
module load iqtree/2.2.0
module load mafft/7.490

echo "=== Multi-Sample Outbreak Analysis ==="

# Step 1: Concatenate all consensus sequences
echo "Step 1: Preparing sequences for phylogenetic analysis..."
cat *.consensus.fasta > all_samples.fasta

# Step 2: Multiple sequence alignment
echo "Step 2: Performing multiple sequence alignment..."
mafft --auto --thread $SLURM_CPUS_PER_TASK all_samples.fasta > aligned_sequences.fasta

# Step 3: Build phylogenetic tree
echo "Step 3: Building phylogenetic tree..."
iqtree2 -s aligned_sequences.fasta -nt $SLURM_CPUS_PER_TASK -bb 1000

# Step 4: Calculate pairwise distances
echo "Step 4: Calculating genetic distances..."
python << 'EOF'
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import pandas as pd

# Read alignment
alignment = AlignIO.read("aligned_sequences.fasta", "fasta")

# Calculate distances
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Convert to DataFrame for easier handling
samples = [record.id for record in alignment]
dist_df = pd.DataFrame(distance_matrix.matrix, 
                      index=samples, 
                      columns=samples)

# Save distance matrix
dist_df.to_csv('genetic_distances.csv')

# Find closely related samples (distance < 0.001)
close_pairs = []
for i, sample1 in enumerate(samples):
    for j, sample2 in enumerate(samples[i+1:], i+1):
        distance = distance_matrix.matrix[i][j]
        if distance < 0.001:  # Very similar sequences
            close_pairs.append([sample1, sample2, distance])

if close_pairs:
    close_df = pd.DataFrame(close_pairs, 
                           columns=['Sample1', 'Sample2', 'Distance'])
    close_df.to_csv('potential_transmission_links.csv', index=False)
    print(f"Found {len(close_pairs)} potential transmission links")
else:
    print("No closely related samples found")
EOF

echo "=== Outbreak Analysis Complete ==="
echo "Results:"
echo "- Phylogenetic tree: aligned_sequences.fasta.treefile"
echo "- Genetic distances: genetic_distances.csv"
echo "- Potential links: potential_transmission_links.csv"
```

---

## Advanced Job Types

### 1. Array Jobs

```bash
#!/bin/bash
#SBATCH --job-name=array_processing
#SBATCH --partition=Main
#SBATCH --array=1-100%10        # 100 jobs, max 10 concurrent
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=00:30:00
#SBATCH --output=array_%A_%a.log

module load python/3.12.3  # Or use system python3

# Use array task ID to process different files
INPUT_FILE="input_${SLURM_ARRAY_TASK_ID}.txt"
OUTPUT_FILE="output_${SLURM_ARRAY_TASK_ID}.txt"

echo "Processing $INPUT_FILE on $(hostname)"
python process_file.py $INPUT_FILE $OUTPUT_FILE

echo "Task $SLURM_ARRAY_TASK_ID completed"
```

### 2. Job Dependencies

```bash
#!/bin/bash
# Submit first job
JOB1=$(sbatch --parsable preprocess.sh)

# Submit second job that depends on first
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 analysis.sh)

# Submit final job that depends on second
sbatch --dependency=afterok:$JOB2 postprocess.sh
```

### 3. Multi-node MPI Job

```bash
#!/bin/bash
#SBATCH --job-name=mpi_job
#SBATCH --partition=Main
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=04:00:00

# module load openmpi  # Check if MPI is available

# Total tasks = nodes * ntasks-per-node = 4 * 16 = 64
echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS total tasks"

mpirun ./my_mpi_program input.dat
```

### 4. Interactive Job

![Interactive vs Batch Jobs](diagrams/interactive_jobs.svg)
*Comparison between interactive and batch job workflows in SLURM*

```bash
# Request interactive session using sinteractive (ILIFU-specific)
sinteractive -c 1 --time 03:00                    # 1 CPU for 3 hours (default)
sinteractive -c 5 --time 5-00:00                 # 5 CPUs for 5 days (maximum)

# Alternative: Use srun for interactive session
srun --partition=Main --cpus-per-task=4 --mem=8GB --time=02:00:00 --pty bash

# Once in interactive session:
module load python/3.12.3  # Or use system python3
python  # Start interactive Python
```

**Note**: Resources on the Devel partition are shared (CPU and memory). For dedicated resources, use `srun` on the Main partition.

### 5. Jupyter Notebook on Compute Node

```bash
#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --output=jupyter_%j.log

module load python/3.12.3  # Or use system python3

# Install jupyter if needed
pip install --user jupyter

# Get node info
NODE=$(hostname -s)
PORT=8888

echo "Starting Jupyter notebook on node $NODE, port $PORT"
echo "SSH tunnel command:"
echo "ssh -N -L ${PORT}:${NODE}:${PORT} ${USER}@training.ilifu.ac.za"

# Start Jupyter
jupyter notebook --no-browser --port=$PORT --ip=0.0.0.0
```

---

## Resource Optimization

### 1. Memory Optimization Examples

#### Low Memory Job

```bash
#!/bin/bash
#SBATCH --job-name=low_mem
#SBATCH --partition=Main
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB           # Conservative memory request
#SBATCH --time=01:00:00

module load python/3.12.3  # Or use system python3

# Process data in chunks to save memory
python << 'EOF'
import pandas as pd

# Read in chunks instead of loading entire file
chunk_size = 10000
results = []

for chunk in pd.read_csv('large_file.csv', chunksize=chunk_size):
    # Process chunk
    processed = chunk.groupby('category').sum()
    results.append(processed)

# Combine results
final_result = pd.concat(results)
final_result.to_csv('output.csv')
EOF
```

#### Memory-intensive Job

```bash
#!/bin/bash
#SBATCH --job-name=high_mem
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB          # High memory for large datasets
#SBATCH --time=04:00:00

module load python/3.12.3  # Or use system python3

# Load large dataset into memory
python << 'EOF'
import pandas as pd
import numpy as np

# Load entire large dataset
df = pd.read_csv('very_large_file.csv')
print(f"Loaded dataset with shape: {df.shape}")

# Memory-intensive operations
correlation_matrix = df.corr()
correlation_matrix.to_csv('correlations.csv')
EOF
```

### 2. Time Optimization

#### Checkpointing Example

```bash
#!/bin/bash
#SBATCH --job-name=checkpointed_job
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --output=checkpoint_%j.log

module load python/3.12.3  # Or use system python3

python << 'EOF'
import pickle
import os
import time

checkpoint_file = 'checkpoint.pkl'

# Try to load previous state
if os.path.exists(checkpoint_file):
    with open(checkpoint_file, 'rb') as f:
        state = pickle.load(f)
    start_iteration = state['iteration']
    results = state['results']
    print(f"Resuming from iteration {start_iteration}")
else:
    start_iteration = 0
    results = []
    print("Starting from scratch")

# Main computation loop
for i in range(start_iteration, 1000):
    # Simulate some work
    time.sleep(1)
    result = i ** 2
    results.append(result)
    
    # Save checkpoint every 100 iterations
    if i % 100 == 0:
        state = {'iteration': i + 1, 'results': results}
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(state, f)
        print(f"Checkpoint saved at iteration {i}")

print("Computation completed")

# Clean up checkpoint file
if os.path.exists(checkpoint_file):
    os.remove(checkpoint_file)
EOF
```

---

## Troubleshooting Examples

### 1. Debug Job Failures

```bash
#!/bin/bash
#SBATCH --job-name=debug_job
#SBATCH --partition=Main
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=00:15:00
#SBATCH --output=debug_%j.log
#SBATCH --error=debug_%j.err

# Enable debugging
set -e  # Exit on any error
set -x  # Print commands as they execute

echo "=== Environment Information ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "User: $(whoami)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM CPUs: $SLURM_CPUS_PER_TASK"

echo "=== Module Information ==="
module list

echo "=== Python Information ==="
module load python/3.12.3  # Or use system python3
which python
python --version

echo "=== Running Script ==="
python my_script.py 2>&1 | tee python_output.log

echo "=== Job Completed ==="
echo "Exit code: $?"
```

### 2. Memory Usage Monitoring

```bash
#!/bin/bash
#SBATCH --job-name=memory_monitor
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=01:00:00

# Function to monitor memory usage
monitor_memory() {
    while true; do
        echo "$(date): Memory usage: $(free -h | grep '^Mem' | awk '{print $3}')"
        sleep 30
    done
}

# Start memory monitoring in background
monitor_memory &
MONITOR_PID=$!

# Load modules and run main task
module load python/3.12.3  # Or use system python3
python memory_intensive_script.py

# Stop monitoring
kill $MONITOR_PID
```

### 3. File Permission Issues

```bash
#!/bin/bash
#SBATCH --job-name=file_check
#SBATCH --partition=Main
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00

echo "=== File System Checks ==="

# Check input files exist and are readable
INPUT_FILES=("input1.txt" "input2.txt" "config.json")

for file in "${INPUT_FILES[@]}"; do
    if [[ -f "$file" ]]; then
        if [[ -r "$file" ]]; then
            echo "✓ $file exists and is readable"
        else
            echo "✗ $file exists but is not readable"
            ls -l "$file"
            exit 1
        fi
    else
        echo "✗ $file does not exist"
        exit 1
    fi
done

# Check output directory is writable
OUTPUT_DIR="results"
if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR" || {
        echo "✗ Cannot create output directory $OUTPUT_DIR"
        exit 1
    }
fi

if [[ -w "$OUTPUT_DIR" ]]; then
    echo "✓ Output directory $OUTPUT_DIR is writable"
else
    echo "✗ Output directory $OUTPUT_DIR is not writable"
    ls -ld "$OUTPUT_DIR"
    exit 1
fi

echo "All file checks passed!"

# Proceed with actual work
python main_script.py
```

---

## Job Submission Scripts

### Batch Submit Multiple Jobs

```bash
#!/bin/bash
# submit_multiple.sh - Submit multiple related jobs

# Array of input files
INPUT_FILES=(data1.txt data2.txt data3.txt data4.txt)

# Submit a job for each input file
for i in "${!INPUT_FILES[@]}"; do
    input_file="${INPUT_FILES[$i]}"
    job_name="process_$(basename $input_file .txt)"
    
    echo "Submitting job for $input_file"
    
    sbatch --job-name="$job_name" \
           --output="${job_name}_%j.log" \
           --export=INPUT_FILE="$input_file" \
           process_template.sh
    
    sleep 1  # Brief pause between submissions
done
```

### Template with Environment Variables

```bash
#!/bin/bash
#SBATCH --job-name=templated_job
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log  # %x = job name, %j = job id

# Use environment variables passed from submission script
echo "Processing file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Parameters: $PARAMS"

module load python/3.12.3  # Or use system python3

# Use the variables in your script
python analysis.py \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_DIR" \
    --params "$PARAMS"
```

---

## Performance Testing Template

```bash
#!/bin/bash
#SBATCH --job-name=performance_test
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --output=perf_%j.log

echo "=== Performance Test Started ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: ${SLURM_MEM_PER_NODE}MB"
echo "Start time: $(date)"

# Record resource usage
echo "=== Initial Resource Usage ==="
free -h
df -h $HOME
df -h /scratch/$USER

module load python/3.12.3  # Or use system python3

# Time the main computation
echo "=== Starting Main Computation ==="
start_time=$(date +%s)

python performance_test_script.py

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "=== Performance Summary ==="
echo "Runtime: ${runtime} seconds"
echo "End time: $(date)"

# Check final resource usage
echo "=== Final Resource Usage ==="
free -h

echo "=== Performance Test Completed ==="
```

This comprehensive set of SLURM examples covers most common use cases and provides templates that can be adapted for specific needs. Each example includes comments explaining the key parameters and concepts.
