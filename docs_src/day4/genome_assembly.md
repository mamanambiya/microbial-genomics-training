# Objectives

1. To perform *De novo* genome assembly and reference-based genome assembly
2. Assess assembly quality using quality metrics

# Target Organisms
## Mycobacterium tuberculosis

- Genome size: ~4.4 Mb
- Characteristics: High GC content ~65.6%, complex secondary structures
- Genes ~4,000 protein-coding genes
- **Clinical relevance: Major global pathogen, drug resistance concerns**
- **Assembly challenges: Repetitive sequences, PE/PPE gene families**

## Vibrio cholerae

- Characteristics: Dual chromosome structure
- Genome size: ~4.0 Mb (two chromosomes: ~3.0 Mb + ~1.1 Mb)
- GC content: ~47.7%
- Genes: ~3,800 protein-coding genes
- **Clinical relevance: Cholera pandemics, epidemic strain tracking**
- **Assembly challenges: Chromosome separation, mobile genetic elements**

# Genome assembly
Genome assembly involves reconstructing a genome from a set of fragmented DNA sequences (reads) obtained from sequencing technologies. 

It aims to piece together the reads to create a continuous sequence that represents the genome of an organism.

## Key Concepts
 - Reads: Fragmented sequences of DNA obtained from sequencing technologies.
 - Contigs: Continuous sequences formed by overlapping reads.
 - Scaffolds: Ordered and oriented sets of contigs, sometimes with gaps, which represent larger regions of the genome.
 - Coverage: The average number of times each base in the genome is sequenced, which affects the accuracy of the assembly.
 - Assembly can be:
    - *De novo* assembly - the construction of the genome from scratch without a reference.
    - Reference-guided assembly - uses a known reference genome to guide the assembly of the new genome.
    
## Steps in Genome Assembly
1. Preprocessing (quality control, adapter sequences and low-quality bases filtering with FASTQC, TRIMMOMATIC/FASTP/BBMap).
2. Assembly
3. Scaffolding (involves use of long-range information from paired-end or mate-pair reads to order and orient contigs into scaffolds (e.g., SSPACE).
4. Gap Filling (use additional reads or assembly techniques to close gaps within scaffolds (e.g., GapCloser).
5. Polishing (correcting sequencing errors and improving the accuracy of the assembly, e.g., Pilon).
6. Evaluation (Assembly Metrics including N50 (length of the contig such that 50% of the total assembly length is in contigs of this length or longer); genome completeness; and accuracy, e.g., QUAST).

## Genome Assembly Algorithms

### Overlap Layout Consensus (OLC)
- Suitable for long-read sequencing (e.g., PacBio, Oxford Nanopore)
- Determines overlap between reads
- Arranges reads based on ovelaps
- Resolves conflict to build the final sequences
- Handles long reads and complex regions well
- Computationally intensive
        
### De Bruijn Graph (DBG)
- Best for short-read sequencing (e.g., illumina).
- Uses small overlapping sequences (*k-mers*) to build a graph, nodes: *k-mers*, and edges: overlaps.
- Fast and memory efficient for large datasets
- Struggles with repetitive regions

### Hybrid Methods
- Combine DBG and OLC strengths to improve assembly quality, especially for error-prone long reads.​

## Assembler Classification by Read Type

### Short-Read Assemblers​
- Short-read assemblers process high volumes of short sequences using DBG algorithms
- Ideal for technologies like Illumina.​
  - SPAdes (widely used for small genomes)
  - Velvet (older though useful)
  - SOAPdenovo (short-read assembler for larger genomes),
  - ABySS (large genomes).

### Long-Read Assemblers​
- Long-read assemblers handle fewer, longer sequences using OLC algorithms
- Can span entire genes and repetitive regions
- Optimized for platforms like PacBio and ONT.​
- Higher error rates (can be corrected by short reads or polishing)

#### Examples of Long-Read Assemblers for Bacterial Genomes
==============================================

1. FLYE (Recommended for most bacterial genomes)
   - Excellent for PacBio and Oxford Nanopore
   - Good repeat resolution
   - Usage: flye --nano-raw reads.fastq --out-dir output --genome-size 4.5m

2. Canu (High accuracy, slower)
   - Gold standard for accuracy
   - Requires significant computational resources
   - Usage: canu -p prefix -d output genomeSize=4.5m -nanopore reads.fastq

3. Unicycler (Hybrid approach)
   - Combines short and long reads
   - Excellent for complete genomes
   - Usage: unicycler -1 short_R1.fq -2 short_R2.fq -l long_reads.fq -o output

4. Raven (Fast, lightweight)
   - Quick assemblies for preliminary analysis
   - Good for large datasets
   - Usage: raven reads.fastq > assembly.fasta

5. NextDenovo (High accuracy for Nanopore)
   - Specialized for Oxford Nanopore data
   - Good error correction
   - Usage: nextDenovo config.txt

**Polishing Tools:**
- Medaka (Nanopore): medaka_consensus -i reads.fastq -d assembly.fasta -o polished
- Pilon (with short reads): pilon --genome assembly.fasta --frags mapped_reads.bam
- Racon: racon reads.fastq mappings.paf assembly.fasta
  
### Hybrid Assemblers​
- Combine short and long reads, integrating DBG and OLC methods for improved accuracy and contiguity.​
- Balances accuracy (short reads) and completeness (long reads, resolving repeats and structural variants)
- Requires careful data integration
- Unicycler, MaSuRCA

### Reference-guided Assembly
- Works well given a closely related reference genome
- Aligns reads to a reference genome and fills gaps
- Faster and less computationally demanding
- May miss novel sequences or structural variations
        
## NOTE
| **Assembly Strategy** | **Subtypes** | **Read Type Compatibility** | **Examples** | **Notes** |
|------------|-------|-------|---------|---------|
| De novo assembly | DBG, OLC, Hybrid | Short, long, hybrid  | Velvet, SPAdes, Canu, Flye, MaSuRCA, UniCycler| No reference genome used |
| Reference-guided assembly | Mapping-based	|Short, long	| BWA, Bowtie2, Novoalign, Minimap2 |	Aligns reads to a known reference genome |

## Assembler Selection Factors​
Choosing an assembler depends on"
- Sequencing technology,
- Project goals, and
- Available computational resources.​
       
## Best Practices for Genome Assembly
- Start with high-quality, high-coverage sequencing data to improve the accuracy of the assembly.
- Try multiple assemblers and compare results, as different tools may perform better for different data types.
- Use iterative rounds of assembly, scaffolding, and polishing to gradually improve the assembly.
- Validate the final assembly using independent data, such as long-read sequencing or optical mapping.
- Keep detailed records of all parameters and steps used in the assembly process for reproducibility.
