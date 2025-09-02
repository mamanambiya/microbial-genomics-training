# Day 6: Nextflow Pipeline Development & Version Control with GitHub

**Date**: September 8, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Pipeline completion and introduction to version control with Git/GitHub

## Overview

Day 6 continues Nextflow pipeline development with a focus on completing the genomic analysis pipeline and introducing basic version control using Git and GitHub. Participants will learn fundamental Git commands to track their pipeline development, create a GitHub repository, and perform basic operations like committing changes and syncing with remote repositories.

## Learning Objectives

By the end of Day 6, you will be able to:

- Complete a functional Nextflow pipeline for genomic analysis
- Initialize a Git repository for your pipeline project
- Make commits to track changes in your pipeline
- Create a GitHub account and repository
- Push your local repository to GitHub
- Pull changes from GitHub to your local machine
- Write a basic README file for your pipeline

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Pipeline development continued and testing* | | Mamana Mbiyavanga |
| **11:30** | **Break** | | |
| **12:00** | *Introduction to Git and GitHub basics* | | Mamana Mbiyavanga |

## Key Topics

### 1. Completing the Nextflow Pipeline
- Finalizing QC, assembly, and annotation processes
- Testing pipeline with sample data
- Debugging common pipeline issues
- Running the complete pipeline end-to-end

### 2. Introduction to Version Control
- Why version control matters for bioinformatics
- Understanding Git concepts (repository, commit, staging)
- Benefits of tracking pipeline development
- Version control best practices

### 3. Basic Git Commands
- `git init` - Initialize a repository
- `git add` - Stage changes
- `git commit` - Save changes with a message
- `git status` - Check repository status
- `git log` - View commit history

### 4. Getting Started with GitHub
- Creating a GitHub account
- Creating your first repository
- Connecting local repository to GitHub
- Understanding remote repositories

### 5. Essential GitHub Operations
- `git push` - Upload changes to GitHub
- `git pull` - Download changes from GitHub
- `git clone` - Copy a repository
- Writing a simple README.md file

## Tools and Software

### Pipeline Development
- **Nextflow** - Workflow orchestration
- **Docker/Singularity** - Container platforms
- **Text editor** - For editing pipeline scripts

### Version Control
- **Git** - Version control system (command line)
- **GitHub** - Online code repository hosting
- **Web browser** - For accessing GitHub website

## Hands-on Exercises

### Exercise 1: Complete the Pipeline (90 minutes)
Finalize and test your Nextflow pipeline.

```groovy
// Complete pipeline structure
workflow {
    // Input channels
    ch_reads = Channel.fromFilePairs(params.reads)
    
    // QC Process
    FASTQC(ch_reads)
    
    // Assembly Process
    SPADES(ch_reads)
    
    // Annotation Process
    PROKKA(SPADES.out.contigs)
    
    // MultiQC Report
    MULTIQC(FASTQC.out.zip.collect())
}
```

Test the pipeline:
```bash
# Run with test data
nextflow run main.nf --reads "data/*_{1,2}.fastq.gz"

# Check outputs
ls -la results/
```

### Exercise 2: Initialize Git Repository (30 minutes)
Set up version control for your pipeline.

```bash
# Navigate to your pipeline directory
cd my-nextflow-pipeline/

# Initialize Git repository
git init

# Check status
git status

# Create .gitignore file
nano .gitignore
# Add these lines:
# work/
# .nextflow/
# .nextflow.log*
# results/
# Save with Ctrl+X, Y, Enter

# Add files to staging
git add main.nf
git add nextflow.config
git add README.md
git add .gitignore

# Make your first commit
git commit -m "Initial commit: Basic Nextflow pipeline"

# View commit history
git log --oneline
```

### Exercise 3: Create GitHub Repository (30 minutes)
Share your pipeline on GitHub.

1. **Create GitHub Account** (if needed):
   - Go to https://github.com
   - Sign up for free account
   - Verify email address

2. **Create New Repository**:
   ```
   - Click "New repository" button
   - Name: my-genomics-pipeline
   - Description: Nextflow pipeline for QC, assembly, and annotation
   - Set as Public or Private
   - DO NOT initialize with README (we already have one)
   - Click "Create repository"
   ```

3. **Connect Local to GitHub**:
   ```bash
   # Add remote repository (replace USERNAME with your GitHub username)
   git remote add origin https://github.com/USERNAME/my-genomics-pipeline.git
   
   # Push to GitHub
   git push -u origin main
   
   # Enter GitHub username and password/token when prompted
   ```

### Exercise 4: Basic Git Workflow (30 minutes)
Practice the edit-add-commit-push cycle.

```bash
# Make changes to your pipeline
nano main.nf
# Add a comment or modify a parameter

# Check what changed
git status
git diff

# Stage the changes
git add main.nf

# Commit with descriptive message
git commit -m "Add parameter for minimum contig length"

# Push to GitHub
git push

# Pull any changes (if working with others)
git pull
```

### Exercise 5: Create README (30 minutes)
Document your pipeline with a README file.

```bash
# Create or edit README
nano README.md
```

Add this content:
```markdown
# My Genomics Pipeline

## Description
A Nextflow pipeline for bacterial genome analysis including:
- Quality control with FastQC
- De novo assembly with SPAdes
- Annotation with Prokka

## Requirements
- Nextflow (version 20.10 or later)
- Docker or Singularity

## Usage
```bash
nextflow run main.nf --reads "data/*_{1,2}.fastq.gz"
```

## Parameters
- `--reads`: Path to paired-end FASTQ files
- `--outdir`: Output directory (default: results)

## Output
- `fastqc/`: Quality control reports
- `assembly/`: Assembled contigs
- `annotation/`: Genome annotations

## Author
Your Name
```

Save and push to GitHub:
```bash
git add README.md
git commit -m "Add README documentation"
git push
```

## Key Concepts

### Version Control Benefits
- **Track Changes**: See what changed, when, and why
- **Backup**: Your code is safe on GitHub
- **Collaboration**: Work with others easily
- **Reproducibility**: Others can use your exact pipeline version

### Git Workflow
```
Working Directory → Staging Area → Local Repository → Remote Repository
     (edit)           (git add)      (git commit)        (git push)
```

### Commit Best Practices
- Write clear, descriptive messages
- Commit logical units of change
- Commit frequently
- Don't commit generated files (use .gitignore)

## Common Issues and Solutions

### Git Configuration
```bash
# Set up your identity (first time only)
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### Authentication Issues
```bash
# If password authentication fails, use personal access token:
# GitHub → Settings → Developer settings → Personal access tokens
# Generate new token with 'repo' permissions
# Use token instead of password when pushing
```

### Common Git Commands
| Command | Purpose | Example |
|---------|---------|---------|
| `git status` | Check repository state | `git status` |
| `git diff` | See changes | `git diff main.nf` |
| `git add .` | Stage all changes | `git add .` |
| `git commit -m` | Save changes | `git commit -m "Fix bug"` |
| `git push` | Upload to GitHub | `git push origin main` |
| `git pull` | Download from GitHub | `git pull origin main` |

## Assessment Activities

### Individual Tasks
- Successfully complete and run the Nextflow pipeline
- Initialize a Git repository for your pipeline
- Create a GitHub account and repository
- Push your pipeline to GitHub
- Create a comprehensive README file

### Group Discussion
- Share GitHub repository URLs with the class
- Discuss pipeline design choices
- Review each other's README files
- Troubleshoot Git/GitHub issues together

## Resources

### Nextflow Resources
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nf-core pipelines](https://nf-co.re/)

### Git/GitHub Resources
- [GitHub Guides](https://guides.github.com/)
- [Git Handbook](https://guides.github.com/introduction/git-handbook/)
- [GitHub Learning Lab](https://lab.github.com/)
- [Pro Git Book (free)](https://git-scm.com/book)

## Looking Ahead

**Day 7 Preview**: Metagenomic Profiling
- Metagenomic sequencing principles
- Quality control for metagenomic data
- Microbiome profiling with R and QIIME2
- Diversity metrics and analysis

---

**Key Learning Outcome**: Completing a functional Nextflow pipeline and establishing version control with Git/GitHub provides the foundation for reproducible, shareable bioinformatics research.