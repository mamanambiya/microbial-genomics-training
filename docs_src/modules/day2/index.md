# Day 2: Introduction to Commandline

**Date**: September 2, 2025  
**Duration**: 09:00-13:00 CAT  
**Focus**: Command line proficiency, M. tuberculosis genomics

## Overview

Day 2 focuses on building strong command line skills essential for bioinformatics work. This day provides the computational foundation needed for all subsequent genomic analyses in the course.

## Learning Objectives

By the end of Day 2, you will be able to:

- Master essential Unix/Linux command line operations for bioinformatics workflows
- Understand M. tuberculosis genomics and co-infection patterns

## Schedule

| Time (CAT) | Topic | Links | Trainer |
|------------|-------|-------|---------|
| **09:00** | *Introduction to command line interface* | [Practical](https://github.com/Arash-Iranzadeh/Microbial-Genomics/blob/main/scripts/practice_unix_commands.sh) | Arash Iranzadeh |
| **11:30** | **Break** | | |
| **12:00** | *Guest talk: MtB and co-infection* | [Speaker Bio](https://docs.google.com/document/d/1ITTQd4jDtgq8gnjYqk-jFHhdCu8JwBjXp_uTDGwBqfY/edit?tab=t.0#heading=h.d9aiuhe8jbzy) | Bethlehem Adnew |

## Key Topics

### 1. Command Line Interface Fundamentals
- Unix/Linux file system navigation
- Essential commands for bioinformatics (grep, awk, sed)
- File manipulation and text processing
- Pipes and command chaining
- Working with compressed files (gzip, tar)
- Shell scripting basics for automation

### 2. M. tuberculosis and Co-infection
- TB genomics and strain typing
- Co-infection patterns and detection
- Clinical implications
- Molecular epidemiology approaches
- Drug resistance mechanisms
- Public health applications

## Tools and Software

### Command Line Tools
- **Bash shell** - Command line interface and scripting
- **GNU coreutils** - Essential Unix utilities (ls, cd, grep, etc.)
- **Text processing** - awk, sed, cut, sort, uniq
- **File compression** - gzip, tar, zip
- **tmux/screen** - Terminal session management
- **rsync/scp** - File transfer and synchronization

## Hands-on Exercises

### Exercise 1: Command Line Fundamentals (90 minutes)
Master essential Unix commands for bioinformatics through practical exercises.

```bash
# Navigate file systems and manipulate files
cd ~/data
ls -la
mkdir analysis_output

# Process text files with Unix tools
grep "^>" sequences.fasta | wc -l  # Count sequences
cat sample.fastq | head -20         # View file contents

# Work with compressed files
gzip large_file.txt
gunzip -c compressed.gz | head

# Use pipes and redirection
cat data.txt | sort | uniq > unique_values.txt
```

## Key Concepts

### Command Line Essentials
- **File system navigation**: Understanding directory structure and paths
- **Text processing**: Using grep, sed, awk for data manipulation
- **Pipes and redirection**: Chaining commands for complex operations
- **Shell scripting**: Automating repetitive tasks
- **Regular expressions**: Pattern matching in bioinformatics data

## Assessment Activities

### Individual Tasks
- Successfully connect to Ilifu HPC system
- Navigate Unix file system and manipulate files
- Complete command line exercises for pathogen genomics data

### Group Discussion
- Share command line tips and tricks
- Discuss HPC resource management strategies
- Troubleshoot connection and job submission issues
- Compare different approaches to batch processing

## Common Challenges

### Command Line Challenges
```bash
# Permission denied errors
chmod +x script.sh    # Make script executable
ls -la                # Check file permissions

# Path issues
echo $PATH            # Check current PATH
export PATH=$PATH:/new/path  # Add to PATH
```

### Command Line Resources
- [Unix for Bioinformatics](https://bioinformatics.ucdavis.edu/research-computing/documentation/unix-basics/)
- [Bash Scripting Guide](https://www.gnu.org/software/bash/manual/)
- [Command Line for Genomics](https://datacarpentry.org/shell-genomics/)

### M. tuberculosis Resources
- [TB-Profiler](https://github.com/jodyphelan/TBProfiler)
- [ReSeqTB](https://platform.reseqtb.org/)
- [TBDB](http://genome.tbdb.org/)

## Guest Lecture: MtB and Co-infection

### Speaker: Bethlehem Adnew

#### Key Topics Covered
- **M. tuberculosis genomics**: Strain diversity and typing methods
- **Co-infection dynamics**: TB-HIV and other respiratory pathogens
- **Diagnostic challenges**: Molecular detection in complex samples
- **Treatment implications**: Drug resistance in co-infected patients
- **Epidemiological insights**: Transmission patterns and control strategies

#### Interactive Discussion Points
- Current challenges in TB diagnosis
- Role of genomics in outbreak investigation
- Future directions in TB research
- Integration of genomic and clinical data

## Looking Ahead

**Day 3 Preview**: 
- Command line proficiency,
- HPC fundamentals
- Quality checking and control with FastQC
- Species identification using Kraken2

---

**Key Learning Outcome**: Mastery of command line operations and HPC infrastructure usage provides the essential computational foundation for all subsequent genomic analyses in the course. 
