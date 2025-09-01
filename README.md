# Microbial Genomics & Metagenomics Training Course

[![Course Website](https://img.shields.io/badge/Website-Live-brightgreen)](https://cidri-africa.github.io/microbial-genomics-training/)
[![Build Status](https://img.shields.io/badge/Build-Testing-yellow)]
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Course Dates](https://img.shields.io/badge/Course%20Dates-September%201--12,%202025-orange)]()

## 🔬 Course Overview

This comprehensive 10-day intensive training course introduces participants to the core principles and practical applications of microbial genomics and metagenomics in clinical and public health contexts. Using real-world datasets and reproducible workflows, participants will learn to analyze genomic diversity, evolutionary relationships, and antimicrobial resistance (AMR) profiles of clinically significant pathogens.

## 🎯 Target Pathogens

- **_Mycobacterium tuberculosis_** - Drug resistance and transmission dynamics
- **_Vibrio cholerae_** - Outbreak investigation and epidemic analysis

## 📅 Course Information

- **Dates**: September 1-12, 2025
- **Time**: 09:00-13:00 CAT daily
- **Duration**: 10 days intensive training
- **Venue**: MAC room, level 2, Health Science UCT, Barnard Fuller Building
- **Address**: Anzio Rd, Observatory, Cape Town, 7935
- **Format**: Hands-on workshops with lectures

## 🎓 Learning Outcomes

After completing this course, participants will be able to:

- **Understand** the principles of pathogen genomics and its role in infectious disease surveillance and control
- **Describe** the workflow of genomic data generation, from sample collection to sequence analysis
- **Apply** bioinformatics tools to analyze pathogen genomic data
- **Interpret** genomic data for epidemiological insights, including outbreak detection and tracking
- **Understand** the role of genomics in identifying antimicrobial resistance mechanisms
- **Implement** reproducible analysis workflows using modern tools and best practices

## 📚 Course Structure

### Week 1: Foundations (Days 1-5)
- **Day 1**: Introduction and genomic surveillance overview
- **Day 2**: Command line, HPC, and quality control  
- **Day 3**: Genomic characterization and assembly
- **Day 4**: Comparative genomics and phylogenetics
- **Day 5**: Metagenomic profiling

### Week 2: Advanced Applications (Days 6-10)
- **Day 6**: Co-infections and variant calling
- **Day 7-8**: Nextflow workflows and pipeline development
- **Day 9**: Bring your own data analysis
- **Day 10**: Presentations and course wrap-up

## 👨‍🏫 Training Team

| Trainer | Role | Expertise |
|---------|------|-----------|
| **Ephifania Geza** | Lead Instructor | Genomic surveillance, AMR analysis, metagenomics |
| **Arash Iranzadeh** | Technical Instructor | Command line, QC, assembly, phylogenomics |
| **Sindiswa Lukhele** | Technical Instructor | Sequencing technologies, PubMLST |
| **Mamana Mbiyavanga** | HPC/Workflow Specialist | High-performance computing, Nextflow pipelines |

**Guest Speaker**: Bethlehem Adnew (M. tuberculosis co-infection, Day 2)

## 🚀 Getting Started

1. **Visit the Course Website**: [https://cidri-africa.github.io/microbial-genomics-training/](https://cidri-africa.github.io/microbial-genomics-training/)

2. **Complete Pre-Course Setup**: Follow the setup guide to prepare your system

3. **Review Prerequisites**: Ensure you have basic command line knowledge

4. **Join Communication Channels**: Slack invitation will be sent before the course

## 📁 Repository Structure

```
├── data/                     # Sample datasets and links
├── docs/                     # Original course materials  
│   ├── slides/              # Presentation materials
│   └── handouts_reference_material/  # Additional resources
├── docs_src/                 # Website source files
│   ├── course/              # Course information pages
│   ├── modules/             # Daily module content  
│   ├── appendices/          # Reference materials
│   └── stylesheets/         # Custom styling
├── notebooks/               # Analysis notebooks
├── scripts/                 # Utility scripts
├── workflows/               # Nextflow pipelines
├── environment/             # Environment setup
├── mkdocs.yml              # Website configuration
└── requirements.txt        # Python dependencies
```

## 🌐 Website Development

This repository uses MkDocs with Material theme to generate a comprehensive training website similar to [Genomics Aotearoa's metagenomics summer school](https://genomicsaotearoa.github.io/metagenomics_summer_school/).

### Build the Website Locally

```bash
# Install dependencies
pip install -r requirements.txt

# Serve locally
mkdocs serve

# Build for deployment  
mkdocs build
```

### Website Features

- 📱 Responsive design optimized for all devices
- 🔍 Full-text search functionality  
- 💡 Interactive code blocks with copy buttons
- 📊 Progress tracking and navigation
- 🎨 Custom styling matching course branding
- ⚡ Fast loading and optimized performance

## 🖥️ Technical Requirements

### Hardware
- Laptop (Linux/macOS preferred, Git Bash for Windows)
- 8GB+ RAM (16GB recommended)
- 100GB+ free storage
- Stable internet connection

### Computing Resources
- Access to Ilifu HPC facility provided
- Pre-configured analysis environments
- All necessary bioinformatics software installed

## 📊 Datasets

The course utilizes real-world clinical datasets including:

- **M. tuberculosis**: Drug-resistant clinical isolates
- **V. cholerae**: Outbreak investigation data

Participants can also bring their own data for analysis during the final days.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact & Support

- **General Inquiries**: training@cidri-africa.org
- **Technical Issues**: Create an issue in this repository
- **Course Coordinator**: Ephifania Geza
- **Website Maintainer**: CIDRI-Africa Technical Team

## 🏛️ Acknowledgments

This course was developed with input from:
- UCT Computational Biology Division  
- CIDRI-Africa training team
- International collaborators and guest speakers
- Open-source bioinformatics community
- Course participants and alumni feedback

---

## 🎯 Quick Links

- **📋 [Course Schedule](https://cidri-africa.github.io/microbial-genomics-training/course/schedule/)** - Detailed daily program
- **🎯 [Learning Objectives](https://cidri-africa.github.io/microbial-genomics-training/course/objectives/)** - What you'll learn
- **⚙️ [Setup Guide](https://cidri-africa.github.io/microbial-genomics-training/course/setup/)** - Pre-course preparation
- **📊 [Datasets](https://cidri-africa.github.io/microbial-genomics-training/datasets/)** - Training data information

**Ready to dive into microbial genomics?** Visit our [course website](https://cidri-africa.github.io/microbial-genomics-training/) to get started!