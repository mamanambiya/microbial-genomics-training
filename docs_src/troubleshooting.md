# Troubleshooting Guide

## Common Issues and Solutions

### Setup and Installation Issues

#### Git Configuration Problems
**Problem**: Git not configured properly
```bash
# Check current configuration
git config --list

# Error: Please tell me who you are
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

**Problem**: SSH key authentication fails
```bash
# Generate new SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"

# Add to SSH agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

# Test connection
ssh -T git@github.com
```

#### Permission Denied Errors
**Problem**: Cannot execute scripts or access files
```bash
# Fix script permissions
chmod +x script.sh

# Fix directory permissions
chmod 755 directory/
chmod -R 755 directory/  # Recursive

# Fix SSH key permissions
chmod 600 ~/.ssh/id_ed25519
chmod 644 ~/.ssh/id_ed25519.pub
```

### Data Analysis Issues

#### FastQC Problems
**Problem**: FastQC fails to run
```bash
# Check Java installation
java -version

# Install Java if missing (Ubuntu/Debian)
sudo apt install default-jdk

# Run FastQC with memory limit
fastqc --memory 4096 *.fastq.gz
```

**Problem**: Out of memory errors
```bash
# Increase memory allocation
fastqc --memory 8192 file.fastq.gz

# Process files individually
for file in *.fastq.gz; do
    fastqc --memory 4096 "$file"
done
```

#### Assembly Issues
**Problem**: SPAdes assembly fails
```bash
# Check available memory
free -h

# Run with memory limit
spades.py --memory 32 -1 R1.fastq.gz -2 R2.fastq.gz -o output/

# Try different k-mer sizes
spades.py -k 21,33,55 -1 R1.fastq.gz -2 R2.fastq.gz -o output/
```

**Problem**: Poor assembly quality (high fragmentation)
```bash
# Check input data quality first
fastqc input_files.fastq.gz

# Try more aggressive trimming
trimmomatic PE input_R1.fastq.gz input_R2.fastq.gz \
    output_R1.fastq.gz output_R1_unpaired.fastq.gz \
    output_R2.fastq.gz output_R2_unpaired.fastq.gz \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:50

# Use careful mode in SPAdes
spades.py --careful -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz -o careful_assembly/
```

### Tool Installation and Dependencies

#### Conda/Mamba Issues
**Problem**: Environment creation fails
```bash
# Update conda
conda update conda

# Clear package cache
conda clean --all

# Create environment with specific Python version
conda create -n genomics python=3.9

# Use mamba for faster solving
mamba create -n genomics python=3.9
```

**Problem**: Package conflicts
```bash
# Create minimal environment first
conda create -n clean_env python=3.9

# Activate and install packages one by one
conda activate clean_env
conda install -c bioconda fastqc
conda install -c bioconda spades
```

#### Docker/Singularity Issues
**Problem**: Permission denied with Docker
```bash
# Add user to docker group
sudo usermod -aG docker $USER

# Log out and back in, then test
docker run hello-world
```

**Problem**: Singularity image won't run
```bash
# Pull image explicitly
singularity pull docker://biocontainers/fastqc:v0.11.9_cv8

# Run with specific bind paths
singularity exec -B /data:/data image.sif fastqc --version

# Check image integrity
singularity verify image.sif
```

### HPC and Remote Access Issues

#### SSH Connection Problems
**Problem**: Connection timed out
```bash
# Test basic connectivity
ping hostname

# Try different port
ssh -p 2222 username@hostname

# Use verbose mode for debugging
ssh -v username@hostname
```

**Problem**: Key exchange failed
```bash
# Generate compatible key
ssh-keygen -t rsa -b 4096

# Specify key explicitly
ssh -i ~/.ssh/specific_key username@hostname

# Check SSH config
cat ~/.ssh/config
```

#### SLURM Job Issues
**Problem**: Job stuck in queue
```bash
# Check queue status
squeue -u $USER

# Check job details
scontrol show job JOBID

# Check partition availability
sinfo
```

**Problem**: Job fails with memory errors
```bash
# Check job output
cat slurm-JOBID.out

# Increase memory request
#SBATCH --mem=32G

# Use multiple cores if available
#SBATCH --cpus-per-task=8
```

### Data Processing Errors

#### File Format Issues
**Problem**: Unexpected file format
```bash
# Check file type
file filename
head filename

# Convert line endings if needed
dos2unix filename

# Check compression
gunzip -t file.gz
```

**Problem**: Corrupt or truncated files
```bash
# Check file integrity
md5sum file.fastq.gz
# Compare with provided checksum

# Test gzip integrity
gunzip -t file.fastq.gz

# Repair if possible (may lose data)
gzip -d file.fastq.gz
gzip file.fastq
```

#### Large File Handling
**Problem**: Running out of disk space
```bash
# Check disk usage
df -h
du -sh directory/

# Clean up temporary files
rm -rf temp/
rm *.tmp

# Compress large files
gzip *.fastq
tar -czf archive.tar.gz directory/
```

**Problem**: Processing very large files
```bash
# Process in chunks
split -l 4000000 large_file.fastq chunk_
# Process each chunk separately

# Use streaming where possible
zcat file.fastq.gz | head -n 1000000 | tool

# Use efficient tools
seqtk sample file.fastq.gz 10000 > sample.fastq
```

### Analysis and Interpretation Issues

#### Resistance Gene Detection
**Problem**: No resistance genes found (expected some)
```bash
# Check assembly quality
quast.py assembly.fasta

# Try multiple databases
abricate --db resfinder assembly.fasta
abricate --db card assembly.fasta
abricate --db argannot assembly.fasta

# Reduce stringency
abricate --minid 80 --mincov 60 assembly.fasta
```

**Problem**: Too many false positives
```bash
# Increase stringency
abricate --minid 95 --mincov 90 assembly.fasta

# Verify hits manually
blast -query resistance_gene.fasta -subject assembly.fasta

# Check for truncated genes
abricate --mincov 95 assembly.fasta
```

#### Phylogenetic Analysis
**Problem**: Tree looks wrong or unrealistic
```bash
# Check sequence alignment quality
aliview alignment.fasta

# Remove problematic sequences
seqtk subseq sequences.fasta good_ids.txt > clean.fasta

# Try different tree method
FastTree -nt alignment.fasta > tree.newick
iqtree -s alignment.fasta -m TEST
```

**Problem**: Low bootstrap support
```bash
# Increase bootstrap replicates
iqtree -s alignment.fasta -bb 1000

# Check for recombination
gubbins alignment.fasta

# Use only core SNPs
snp-sites -c alignment.fasta > core_snps.fasta
```

### Performance and Resource Issues

#### Memory Management
**Problem**: Out of memory errors
```bash
# Check memory usage
free -h
top

# Limit memory usage
ulimit -v 8000000  # Limit to ~8GB

# Use memory-efficient tools
minimap2 instead of BWA-MEM for large references
```

**Problem**: Process running too slowly
```bash
# Use multiple cores
tool -t 8 input output

# Optimize I/O
# Use local storage instead of network drives
cp data /tmp/
cd /tmp/
# Run analysis
cp results back/to/network/storage
```

#### Storage Management
**Problem**: Quota exceeded
```bash
# Find large files
find . -size +100M -ls

# Clean up intermediate files
rm *.sam  # Keep only BAM files
rm temp_*

# Compress old data
tar -czf old_analysis.tar.gz old_directory/
rm -rf old_directory/
```

## Getting Help

### Before Asking for Help

1. **Check error messages carefully** - Often contain specific solutions
2. **Search documentation** - Tool manuals usually have troubleshooting sections
3. **Try simple test cases** - Use small datasets to isolate problems
4. **Check system resources** - Memory, disk space, permissions

### How to Ask for Help

#### Include Essential Information
- **Exact error message** (copy-paste, don't retype)
- **Command that failed** (exact command with parameters)
- **System information** (OS, tool versions)
- **Input file details** (size, format, sample content)

#### Good Help Request Example
```
Subject: SPAdes assembly fails with error code 1

I'm running SPAdes on paired-end M. tuberculosis data:
Command: spades.py -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o spades_out/

Error message:
"== Error ==  system call for: ['/usr/bin/python3', '/opt/spades/bin/spades_init.py'] finished abnormally, err code: 1"

System: Ubuntu 20.04, SPAdes v3.15.3
Input files: 2x150bp Illumina, ~50x coverage, 2.3GB total
Available memory: 32GB
Disk space: 500GB free

I've tried with --careful flag and different k-mer sizes but get the same error.
```

### Support Resources

#### Course Support
- **Instructors**: Available during course hours
- **Slack Channel**: `#troubleshooting`
- **Office Hours**: Daily 17:00-18:00 (course week)
- **Peer Support**: Encouraged among participants

#### Online Resources
- **Biostars**: General bioinformatics Q&A
- **Stack Overflow**: Programming and command line issues
- **Tool Documentation**: Always check official documentation
- **Galaxy Training**: Alternative tutorials and explanations

#### Emergency Contacts
- **Technical Issues**: tech-support@course.org
- **Data Access Problems**: data-admin@course.org
- **General Questions**: instructors@course.org

### Prevention Tips

#### Best Practices
1. **Test with small datasets first**
2. **Keep detailed logs of commands**
3. **Use version control for scripts**
4. **Regular backups of important results**
5. **Document your workflow steps**

#### Common Pitfalls to Avoid
- Running analysis without checking input quality
- Using inappropriate parameters for your data type
- Ignoring error messages and logs
- Not checking intermediate results
- Working in directories with spaces in names
- Not backing up important data

---

**Remember**: Most bioinformatics problems have been encountered before. Don't hesitate to search online and ask for help - the community is generally very supportive!