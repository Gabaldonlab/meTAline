# Description
MeTAline, is a snakemake pipeline for metagenomics analysis. MeTAline, facilitates an efficient workflow to preprocess short reads metagenomics data: from read trimming and filtering, through host read subtraction to taxonomic classification using both k-mer and gene marker-based approaches, and functional profiling of the samples.

---

# Main Features
5 main rules:
- trimming
- host_depletion
- kmer_taxonomy
- RAnalysis
- and BioBakery.

trimming rule:
Essential to ensure high-quality data and is based on trimmomatic, which removes low-quality bases and adapter sequences from the sequencing reads. This rule also includes a quality assessment of the trimmed reads by fastqc and an additional step which concatenates the trimmed reads of the same sample from different lanes into a single file for further processing. This step is necessary for samples sequenced across multiple lanes, a common practice used to avoid lane-specific biases.

host_depletion rule:
Automatically activated when a reference genome is provided as an input, and removes host reads by aligning the trimmed reads to the reference using hisat2. This step outputs a fastq file containing only the unmapped reads, that are subsequently used for the taxonomy assignment. If no reference genome is indicated, all trimmed reads are used for the taxonomy assignment.

kmer_taxonomy rule:
Runs Kraken2 to generate a report with the k-mer based classification using a specific database. Additionally, this rule includes sub-rules that allow the user to visualize taxonomic data in Krona and extract reads of interest such as unclassified reads (the default option).

RAnalysis rule:
Converts Kraken2 reports to Biom format, which is then further transformed into the phyloseq format to ensure compatibility with standard microbiome analysis tools. This step also provides subroutines that run R scripts for calculating alpha diversity metrics and visualizing the relative abundances of the top 25 taxa at the phylum level.

BioBakery rule:
Conducts both taxonomic and functional profiling using gene marker-based tools, specifically MetaPhlAn4 for taxonomic profiling and HUMAnN for functional profiling. MetaPhlan4 rule creates one output for the general microbial profiling and another one containing a viral profiling (“.vsc.txt”). The functional profiling outputs three tabular files: genefamilies.tsv (gene families in the sample), pathabundance.tsv (pathway abundances) and pathcoverage.tsv (presence or absence of the pathways).

---

# New Features
N/A

---

# Modified features
N/A

---

# Bug fixes
N/A

---

# Build from source

## Docker image:
```bash
make docker-image
make docker-container
```

## Singularity image
```bash
make singularity
```

---

# Using pre-built Singularity image
[Singularity usage](https://github.com/Gabaldonlab/meTAline/blob/main/README.md#singularity-image-usage)

---

SHA256 Checksums
```bash

```
