# Day 1: Introduction & Setup

**Duration**: Full day (09:00-13:00)  
**Focus**: Course introduction, environment setup, command line basics

## Overview

Welcome to the Microbial Genomics & Metagenomics Training Course! Day 1 introduces you to the course structure, sets up your computing environment, and ensures everyone has the foundational skills needed for the intensive training ahead.

## Learning Objectives

By the end of Day 1, you will be able to:

- Navigate the Unix/Linux command line efficiently
- Use Git for version control in research projects  
- Set up and configure your analysis environment
- Understand the course structure and expectations
- Connect to high-performance computing resources

## Schedule

| Time | Topic | Instructor |
|------|-------|------------|
| 09:00-09:30 | Course Introduction & Logistics | Ephifania Geza |
| 09:30-10:30 | [Command Line Basics](command-line.md) | Arash Iranzadeh |
| 10:30-10:45 | *Coffee Break* | |
| 10:45-11:45 | [Version Control with Git](git-basics.md) | Mamana Mbiyavanga |
| 11:45-12:30 | [Environment Setup](environment.md) | All Instructors |
| 12:30-13:00 | Q&A and Troubleshooting | All Instructors |

## Key Topics

### 1. Course Introduction
- Course objectives and structure
- Introduction to trainers and participants
- Overview of datasets and tools
- Assessment and certification process

### 2. Command Line Fundamentals
- Navigation and file management
- Text processing and pipes
- Process management
- Shell scripting basics

### 3. Version Control with Git
- Git concepts and workflow
- Repository management
- Collaboration and branching
- Integration with research projects

### 4. Environment Configuration
- SSH key setup and HPC access
- Text editor configuration
- Tool installation and testing
- Container setup (Docker/Singularity)

## Prerequisites

Before starting Day 1, ensure you have:

- [ ] Completed the [pre-course setup](../../course/setup.md)
- [ ] Installed Git and configured user information
- [ ] Generated SSH keys
- [ ] Tested basic command line operations

## Hands-on Activities

### Exercise 1: Command Line Navigation
Practice essential Unix commands and explore the course directory structure.

### Exercise 2: Git Repository Setup  
Create your first repository and practice basic Git workflows.

### Exercise 3: HPC Connection
Connect to the course high-performance computing resources.

### Exercise 4: Environment Testing
Verify all tools and connections are working properly.

## Datasets Introduced

- **Sample data** for command line practice
- **Small genomic datasets** for Git exercises
- **Test sequences** for environment verification

## Tools Introduced

- **Bash shell** for command line operations
- **Git** for version control
- **SSH client** for remote connections
- **Text editors** (nano/vim/VS Code)

## Resources

### Essential Reading
- [Unix Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/)
- [Git Tutorial](https://git-scm.com/docs/gittutorial)
- [SSH Essentials](https://www.digitalocean.com/community/tutorials/ssh-essentials-working-with-ssh-servers-clients-and-keys)

### Cheat Sheets
- [Command Line Cheat Sheet](../../appendices/command-reference.md)
- [Git Cheat Sheet](https://education.github.com/git-cheat-sheet-education.pdf)

### Video Tutorials
- Command Line Basics (link TBD)
- Git Fundamentals (link TBD)

## Common Issues

### SSH Connection Problems
```bash
# Check SSH key permissions
chmod 600 ~/.ssh/id_ed25519
chmod 644 ~/.ssh/id_ed25519.pub

# Test SSH agent
ssh-add -l
```

### Git Configuration Issues
```bash
# Verify Git setup
git config --list --global

# Reset if needed
git config --global user.name "Your Name"
git config --global user.email "email@example.com"
```

### File Permission Problems
```bash
# Fix common permission issues
chmod +x script.sh
chmod -R 755 directory/
```

## Assessment

Day 1 assessment includes:

- **Practical exercises** completed during sessions
- **Environment setup checklist** verification
- **Basic command demonstration** at end of day
- **Participation** in group activities

## Looking Ahead

**Day 2 Preview**: We'll dive into data analysis fundamentals, including:
- Quality control of sequencing data
- Read processing and filtering
- Introduction to genome assembly

## Support

### During the Day
- Instructors available for one-on-one help
- Peer support encouraged  
- Shared troubleshooting document

### After Hours
- Slack channel for questions
- Recorded session materials
- Office hours scheduling

---

**Remember**: Don't worry if everything doesn't click immediately. Day 1 is about building foundations – we'll reinforce these concepts throughout the course!