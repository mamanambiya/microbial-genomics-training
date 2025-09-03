# Project-Specific Claude Instructions - Microbial Genomics Training

## Project Overview
This is a 10-day intensive training course on microbial genomics and metagenomics for the CIDRI-Africa program at UCT. The training runs from September 1-12, 2025.

## Repository Structure
- `docs_src/` - Main documentation source files
  - `course/` - Course overview, schedule, prerequisites, objectives
  - `modules/` - Individual day modules (day1-day10)
  - `day2/`, `day3/` - Specific day materials and practicals
- `mkdocs.yml` - MkDocs configuration for documentation site

## Important Guidelines

### Module Updates
- **Always keep modules synchronized with the main schedule** (`docs_src/course/schedule.md`)
- Each day module should have matching:
  - Schedule times
  - Topic titles
  - Trainer names
  - Practical/resource links

### Schedule Structure
- Course runs 09:00-13:00 CAT daily
- Coffee break at 11:00 or 11:30
- Day modules link to detailed content in `docs_src/modules/dayX/index.md`

### Git Configuration
- Use `git-work` command to switch to UCT work profile before pushing
- Remote: `git@github-mamanambiya:mamanambiya/microbial-genomics-training.git`
- Main branch: `main`
- Upstream: `CIDRI-Africa/microbial-genomics-training`

### Trainer Materials
External practicals are hosted on GitHub:
- Arash Iranzadeh: `https://github.com/Arash-Iranzadeh/Microbial-Genomics/`
- Local materials in `docs_src/day2/`, `docs_src/day3/`, etc.

### Content Standards
- Use clear headings and consistent formatting
- Include learning objectives for each module
- Provide both theoretical content and practical exercises
- Add links to slides, practicals, and notes where available
- Include assessment activities and common challenges sections

### Day-Specific Focus Areas
- **Day 1**: Introduction and setup
- **Day 2**: Command line and HPC introduction
- **Day 3**: HPC continued, QC, species identification, genome assembly
- **Day 4**: Genome annotation, AMR detection, MLST
- **Day 5**: Nextflow pipeline development
- **Day 6**: Pipeline development continued, Git/GitHub
- **Day 7**: Metagenomic profiling
- **Day 8**: Comparative genomics (pangenomics, phylogenomics)
- **Day 9**: Mobile genetic elements, participant data analysis
- **Day 10**: Wrap-up and presentations

## Documentation Building
- Site uses MkDocs with Material theme
- Navigation defined in `mkdocs.yml`
- Build with: `mkdocs build`
- Serve locally with: `mkdocs serve`

## Quick Commands
```bash
# Switch to work profile
git-work

# Pull latest changes
git pull origin main

# Push changes
git push origin main

# Check module synchronization
grep "09:00" docs_src/course/schedule.md docs_src/modules/*/index.md
```

## Key Contacts
- **Lead Instructor**: Ephifania Geza
- **Technical Instructors**: Arash Iranzadeh, Sindiswa Lukhele
- **HPC/Workflow Specialist**: Mamana Mbiyavanga

## Notes
- Participants bring their own data for analysis
- MAC room, level 2, Health Science UCT, Barnard Fuller Building
- Ilifu HPC cluster is used for computational work