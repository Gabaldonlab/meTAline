[![DOI](https://zenodo.org/badge/431438117.svg)](https://zenodo.org/badge/latestdoi/431438117)

# meTAline: metagenomics Taxonomic Assignment pipeline

meTAline is a taxonomic assignment pipeline implemented in SnakeMake for WGS metagenomics data.

---

## Build and deploy Singularity image

**NOTE: YOU WILL ALWAYS HAVE TO BUILD THE SINGULARITY IMAGES ON YOUR LOCAL MACHINE, BECAUSE IT REQUIRES SUDO AND THAT YOU WON'T HAVE IN YOUR REMOTE HPC ENVIRONMENT!**

To be able to build the Singularity image of MetaLine, you will need to install Singularity first:

1. Download from [HERE](https://github.com/sylabs/singularity/releases/tag/v4.1.5) the corresponding installation package (_.deb or _.rpm) to your operating system. (E.G.: Ubuntu 24.04 needs the "singularity-ce_4.1.5-noble_amd64.deb" file).
2. Install it from the downloaded package. Example command for the above mentioned version:

```bash
sudo apt install ./singularity-ce_4.1.5-noble_amd64.deb
```

3. Test if the installation is OK.

```bash
singularity --version
```

Now do the build:

1. cd into the root directory of MeTAline (where the Dockerfile and the Makefile is)
2. Run the following command (**REQUIRES "sudo" PERMISSIONS! + This might take several minutes!**):
   **NOTE: This will install a SingularityCE provided tool with pip3 (spython) to translate the Dockerfile to Singularity!**

```bash
make singularity
```

3. Copy the output **metaline.sif** file to your cluster with the **rsync** command. Example:

```bash
rsync ./metaline.sif <cluster_user>@<cluster_transfer_address>:</path/to/your/destination>
```

---

## Usage of Singularity image

1. MeTAline uses Snakemake under the hood using .json configuration files, so first you will need to generate the said file with the following command example:

```bash
singularity run --cleanenv metaline.sif metaline-generate-config \
    --configFile my_metaline_config.json \
    --extension fq.gz \
    --basedir </path/to/your/metaline/output/directory> \
    --reads-directory </path/to/your/directory with the raw reads data> \
    --reference-genome </path/to/your/directory with the reference genome data> \
    --krakendb </path/to/your/directory with the kraken database> \
    --sample-barcode my_metaline_job \
    --fastqs <the prefix of your fastq files> \
    --metaphlan_db </path/to/the/metaphlan database> \
    --metaphlan_Index <Index of the metaphlan4 database.> \
    --n_db <Humann database to do the nucleotide search (based on already built annotations.)> \
    --protein_db <Humann database to do the translation search (by deafult this is by-passed).>
```

For further information about the flags, run:
```bash
singularity run --cleanenv metaline.sif metaline-generate-config -h
```

2. Run the pipeline:

```bash
singularity run --cleanenv metaline.sif \
    metaline \
    -r all \
    -j 16 \
    --configfile my_metaline_config.json
```

**OR test it with dry run only!**

```bash
singularity run --cleanenv metaline.sif \
    metaline \
    -r all \
    -j 16 \
    --configfile my_metaline_config.json \
    -n
```

_The -n option is for a dry run, avoding the execution of any job yet. -j is the number of threads provided to the pipeline. In this specific example, we selected the rule all using the parameter -r. You can use different configs by adding the --configfile parameter but it requires that you comment the configfile argument in the main meTAline.smk._

_More information regarding Snakemake and its commands can be found through Snakemake [documentation](https://snakemake.readthedocs.io/en/stable/index.html)._

---

## Target rules

The target rules currently available to use are:

-   rule **all**: From trimming and filtering of the reads, to host removal (if indicated the host genome) to taxonomy assignment based on kraken2 and metaphlan4. Functional prediction based on Humann. It also includes other functuionalities such as representation of the taxonomy assignment by krona, extraction of desired reads, etc.
-   rule **krona_and_reads**: This rule is used if you already have the kraken2 taxonomic assignment and you want to generate the krona images and extract reads (e.g extraction of unclassified reads).
-   rule **taxonomy_assignment:**: This rule is used if you already have the concatenated/unmapped reads and you want to perform the taxonomic assignment without running the initial steps of the pipeline.
-   rule **biom_format**: This rule is used if you already have the Kraken2 assignment and you want to pass from their reports to biom format for further analysis
-   rule **R_Analysis**: This rule is used to make a basic R analysis and plotting alpha diversity.
-   Rule **BioBakery**: This rule is to perform the taxonomy and functional profiling using Biobakery tools based on gene markers: Metaphlan4 and Humann.

---

## Running in HPC environment:

The MeTAline pipeline is intended to be used in HPC environment, because of its high RAM runtime requirement.
You can find an example job definition file for the SLURM schedular (usually this is used in HPCs) in **./example_hpc_sbatch.job** file.

**EXAMPLE'S CONTENT**

```bash
#!/bin/bash

#SBATCH --job-name=metaline_singularity_test # Job’s name (to trace in)
#SBATCH --output=/hpc-cluster-data/projects/my_project/current/metaline_testing/metaline-test-output/err_out/jobnm_%j.out # Output file, where
#SBATCH --error=/hpc-cluster-data/projects/my_project/current/metaline_testing/metaline-test-output/err_out/jobnm_%j.err # File where the error is written

#SBATCH --ntasks=1 # The number of parallel tasks
#SBATCH --cpus-per-task=24 # Number of CPUs per run task
#SBATCH --tasks-per-node=1 # The number of allocated task/node

#SBATCH --qos=your_quality_of_service # The queue for the job
#SBATCH --account=my_account
#SBATCH --time=24:00:00 # The time you request for the job
#SBATCH --constraint=highmem # To run in highmem nodes

module load singularity

singularity run --cleanenv metaline.sif metaline-generate-config \
    --configFile test_kraken2_standard.json \
    --extension fq.gz \
    --basedir /hpc-cluster-data/projects/my_project/metaline_testing/metaline-test-output \
    --reads-directory /hpc-cluster-data/projects/my_project/metaline_testing/raw_reads_data \
    --reference-genome /hpc-cluster-data/projects/my_project/metaline_testing/reference_genome \
    --krakendb /hpc-cluster-data/projects/my_project/metaline_testing/WGS/KRAKEN2_DB/microbiome_db \
    --sample-barcode metaline_test_barcode \
    --fastqs PONSJIG_165.E5-NEB-UDI \
    --metaphlan_db /hpc-cluster-data/projects/my_project/metaline_testing/Metaphlan4/db \
    --metaphlan_Index mpa_testing_index \
    --n_db /hpc-cluster-data/projects/my_project/metaline_testing/WGS/humann_db \
    --protein_db /hpc-cluster-data/projects/my_project/metaline_testing/WGS/uniref

singularity run --cleanenv metaline.sif metaline \
    -r all \
    -j 16 \
    --configfile test_kraken2_standard.json
```

You can run it with:

```bash
sbatch ./example_hpc_sbatch.job
```

---

## Output directory example:

```bash
metaline-test-output/
├── BAM
│   ├── my_metaline_job.hisat2.bam
│   ├── my_metaline_job.hisat2.bam.bai
│   └── my_metaline_job.unmapped.fastq.gz
├── Benchmark
│   ├── 20240815.230537.312432_my_metaline_job.PONSJIG_165.E5-NEB-UDI.trimming.benchmark.txt
│   ├── 20240815.230537.312432_my_metaline_job.alignment.benchmark.txt
│   ├── 20240815.230537.312432_my_metaline_job.biom.benchmark.txt
│   ├── 20240815.230537.312432_my_metaline_job.concat.benchmark.txt
│   ├── 20240815.230537.312432_my_metaline_job.extracted_reads.benchmark.txt
│   ├── 20240815.230537.312432_my_metaline_job.fastqc.benchmark.txt
│   ├── 20240815.230537.312432_my_metaline_job.krona.benchmark.txt
│   ├── 20240815.230537.312432my_metaline_job.R_alpha_div.benchmark.txt
│   ├── 20240815.230537.312432my_metaline_job.R_barplot.benchmark.txt
│   ├── 20240815.230537.312432my_metaline_job.R_conversion.benchmark.txt
│   └── my_metaline_job.kraken2.benchmark.txt
├── EXTRACTED_FASTA
│   ├── my_metaline_job.1.fastq.gz
│   ├── my_metaline_job.2.fastq.gz
│   ├── my_metaline_job1.fastq
│   └── my_metaline_job2.fastq
├── KRAKEN_ASSIGN
│   ├── my_metaline_job.Kraken2.biom
│   ├── my_metaline_job.kraken2.report
│   └── my_metaline_job.kraken2.txt
├── KRONA_HTML
│   ├── my_metaline_job.html
│   └── my_metaline_job.krona
├── METAPHLAN4
│   ├── my_metaline_job.bz2
│   ├── my_metaline_job.vsc.txt
│   ├── my_metaline_job_profiled.txt
│   └── my_metaline_jobsam.bz2
├── R_ANALYSIS
│   ├── my_metaline_job.rds
│   ├── my_metaline_job_Barplot_phyla.jpeg
│   └── my_metaline_job_alpha_div.csv
├── TRIMMOMATIC
│   ├── PONSJIG_165.E5-NEB-UDI_1_paired.fq.gz
│   ├── PONSJIG_165.E5-NEB-UDI_1_unpaired.fq.gz
│   ├── PONSJIG_165.E5-NEB-UDI_2_paired.fq.gz
│   ├── PONSJIG_165.E5-NEB-UDI_2_unpaired.fq.gz
│   ├── my_metaline_job.1.fastq.gz
│   ├── my_metaline_job.2.fastq.gz
│   └── my_metaline_job_qc
│       ├── my_metaline_job.1_fastqc.html
│       ├── my_metaline_job.1_fastqc.zip
│       ├── my_metaline_job.2_fastqc.html
│       └── my_metaline_job.2_fastqc.zip
├── WGS_logs
│   ├── 20240815.230537.312432.PONSJIG_165.E5-NEB-UDI.trimmomatic.log
│   ├── 20240815.230537.312432.R_alpha_div.log
│   ├── 20240815.230537.312432.R_barplot.log
│   ├── 20240815.230537.312432.R_phyloseq.log
│   ├── 20240815.230537.312432.biom.log
│   ├── 20240815.230537.312432.concatenation.log
│   ├── 20240815.230537.312432.extracted_reads.log
│   ├── 20240815.230537.312432.fastqc.log
│   ├── 20240815.230537.312432.illumina_kraken2.log
│   ├── 20240815.230537.312432.krona.log
│   ├── 20240815.230537.312432_my_metaline_job.alignment.log
│   ├── 20240815.230537.312432_my_metaline_job.all.rule.log
│   ├── 20240815.230537.312432_my_metaline_job.bam2fasta.log
│   ├── 20240815.230537.312432_my_metaline_job.index_bam.log
│   └── 20240815.231914.859216.illumina_kraken2.log
└── err_out
    ├── jobnm_4734683.err
    └── jobnm_4734683.out
```

---

## Benchmark

At the end of the procedure you might want to know the resources that were used. A folder called "Benchmark" will be created containing for each of the rules the following parameters:

**s** --> Running time in seconds

**h: m: s** --> Running time in hours, minutes and seconds

**max_rss** --> "Maximum Resident Set Size". Total amount of physical memory.

**max_vms** --> "Maximum Virtual Memory Size". Total amount of virtual memory.

**max_uss** --> "Unique Set Size". Memory that is unique to a process.

**max_pss** --> "Proportional Set Size". Amount of memory shared with other processes.

**io_in** --> MB read.

**io_out** --> MB written.

**mean_load** --> CPU usage over time / Total running time.

---

## Setup Singularity image for debugging

Adding print statements in the source code, then each time recompile the Singularity image can be cumbersome, because an image is **inmutable**.
To solve this issue, we can use a **sandbox** version of the already compiled Singularity image, by running:

```bash
make metaline-debug-sandbox
```

**NOTE: THIS COMMAND REQUIRES sudo!**

This command will unpack the compiled image into a directory ("./metaline-debug"), which then can be modified and run as a normal image, E.G.:

```bash
singularity run --cleanenv ./metalinedebug ...
```

If you want to build again an **inmutable** image from this "./metaline-debug" directory, you can use the following command:

```bash
make metaline-recompile-from-debug-sandbox
```

Unfortunately, IDE debuggers with Singularity images don't work always, but the majority of the pipeline uses Python scripts, so you will still have **pdb (Python DeBugger)** in your toolbox to initialize an interactive debugger session.

To set a breakpoint in one of the scripts, just call the **breakpoint** method:

```bash
breakpoint()
```

---

## 4. Citations / Acknowledgments

### 4.1. Snakemake

[Sustainable data analysis with Snakemake](<(https://doi.org/10.12688/f1000research.29032.1)>)

### 4.2. Kraken2

(Improved metagenomic analysis with Kraken 2)[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0]

### 4.3. KrakenTools

[KrakenTools - Github](https://github.com/jenniferlu717/KrakenTools)

[Lu J, Rincon N, Wood D E, Breitwieser F P, Pockrandt C, Langmead B, Salzberg S L, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)]
(https://www.nature.com/articles/s41596-022-00738-y)

Relevant paper for usage of KrakenTools:

1. [Kraken 1](https://github.com/DerrickWood/kraken)
2. [Kraken 2](https://github.com/DerrickWood/kraken2)
3. [KrakenUniq](https://github.com/fbreitwieser/krakenuniq)
4. [Bracken](https://github.com/jenniferlu717/Bracken)

### 4.4. GNU Parallel

[Tange, O. (2021, August 22). GNU Parallel 20210822 ('Kabul'). Zenodo.](https://doi.org/10.5281/zenodo.5233953)

### 4.5. Trimmomatic

[Trimmomatic: a flexible trimmer for Illumina sequence data](https://doi.org/10.1093/bioinformatics/btu170)

### 4.6 Bracken

[Bracken: estimating species abundance in metagenomics data](https://peerj.com/articles/cs-104/)

### 4.7 Krona

[Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.](http://www.ncbi.nlm.nih.gov/pubmed/21961884)

### 4.8 HTSlib

[HTSlib: C library for reading/writing high-throughput sequencing data; James K Bonfield, John Marshall, Petr Danecek, Heng Li, Valeriu Ohan, Andrew Whitwham, Thomas Keane, Robert M Davies GigaScience, Volume 10, Issue 2, February 2021, giab007](https://doi.org/10.1093/gigascience/giab007)

### 4.9 Samtools

[Twelve years of SAMtools and BCFtools; Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li GigaScience, Volume 10, Issue 2, February 2021, giab008](https://doi.org/10.1093/gigascience/giab008)

### 4.10 Hisat2

[Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://www.nature.com/articles/s41587-019-0201-4)

### 4.11 MetaPhlAn

[Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4.](https://doi.org/10.1038/s41587-023-01688-w)

### 4.12 Humann

[Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3](https://doi.org/10.7554/eLife.65088)

### 4.13 FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
