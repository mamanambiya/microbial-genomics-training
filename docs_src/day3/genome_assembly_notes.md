# Objectives
1. To perform *De novo* genome assembly and reference-based genome assembly
2. Assess assembly quality using quality metrics

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
- uses small overlapping sequences (*k-mers*) to build a graph, nodes: *k-mers*, and edges: overlaps.
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
- Canu, Celera, Flye
  
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
| De novo assembly | DBG, OLC, Hybrid | Short, long, hybrid  | Velvet, SPAdes, Canu, Flye, MaSuRCA | No reference genome used |
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
