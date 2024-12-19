[![DOI](https://zenodo.org/badge/431438117.svg)](https://zenodo.org/badge/latestdoi/431438117)
*Paper in preparation.*

# MeTAline: a snakemake pipeline for the study of metagenomes

MeTAline, is a snakemake pipeline for metagenomics analysis. MeTAline, facilitates an efficient workflow to preprocess short reads metagenomics data: from read trimming and filtering, through host read subtraction to taxonomic classification using both k-mer and gene marker-based approaches, and functional profiling of the samples.

<div align="center">
  <img src="https://github.com/user-attachments/assets/1c248a72-625f-480a-97ad-9022714a3dbc" width="60%" height="60%">
</div>
---

## Build and deploy Singularity image

**IMPORTANT: If you are going to run the pipeline in marenostrum5, skip this section!**

*Note: you will always have to build the singularity images on your local machine, because it requires sudo and that you won't have in your remote hpc environment!*

To be able to build the Singularity image of MeTAline, you will need to install Singularity first:

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

#If you run the pipeline in marenostrum5, remember to specify the shared sif image: ~/project/pipelines/meTAline/meTAline-0.8.0-alpha/metaline.sif

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

-   rule **all**: From trimming and filtering of the reads, to host removal (if provided the indexed host genome to make an alignment of the reads) to taxonomy assignment based on both a k-mer and gene marker approach. Functional prediction based on Humann. It also includes other functionalities such as representation of the taxonomy assignment by krona, extraction of desired reads, etc.

- rule **trimming**: Trimming of the reads, quality assessment and concatenation (useful when having samples that are sequenced in different sequencing lanes).
- rule **host_depletion**: This rule is to take trimmed reads and align them to an indexed reference genome, for host substraction.
- rule **kmer_taxonomy**: This rule is used if you already have the concatenated/unmapped reads and you want to perform the taxonomic assignment based on k-mers, without running the initial steps of the pipeline.
- rule **RAnalysis**: This rule is used to make a basic R analysis and plotting alpha diversity.
- rule **BioBakery**: This rule is to perform the taxonomy and functional profiling using Biobakery tools based on gene markers: Metaphlan4 and Humann.

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

Available Kraken2 databases at MN5:

```bash
KRAKEN2_DB_COMPLETE	(v1, 20210308) -->	~/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE
eupathDB_kraken2_prebuilt_db	(v1, 20211004) --> ~/project/pipelines/WGS/KRAKEN2_DB/eupathDB_kraken2_prebuilt_db
Broad_fungal_microbiome_db	(v1, 20230512	--> ~/project/pipelines/WGS/KRAKEN2_DB/Broad_fungal_microbiome_db
Broad_gut_microbiome_db	(v1, 20240111	)	-->	~/project/pipelines/WGS/KRAKEN2_DB/Broad_gut_microbiome_db
HumGut_DB	(v1	, 20210702	)	-->	~/project/pipelines/WGS/KRAKEN2_DB/HumGut_DB
HumGut_DB_plus_human	(v1	, 20210709) -->		~/project/pipelines/WGS/KRAKEN2_DB/HumGut_db_plus_human

```

Other templates to run the pipeline for large-scale datasets (using array of greasies for parallelization) are found in "<https://github.com/Gabaldonlab/meTAline/tree/main/meTAline/Helper_scripts_MN5>"

---

## Test sample and output directory example:

Here, in the test folder (<https://github.com/Gabaldonlab/meTAline/tree/main/test_input>) we provide a test mock sample ( V300091236_L01_100_1.fq.gz , V300091236_L01_100_2.fq.gz) to try the pipeline.

The expected output to obtain running all the rules:

```bash
metaline-test-output/
├── BAM
│   ├── mock.hisat2.bam
│   ├── mock.hisat2.bam.bai
│   └── mock.unmapped.fastq.gz
├── Benchmark
│   ├── 20241108.102621.871345_mock.alignment.benchmark.txt
│   ├── 20241108.102621.871345_mock.biom.benchmark.txt
│   ├── 20241108.102621.871345_mock.concat.benchmark.txt
│   ├── 20241108.102621.871345_mock.extracted_reads.benchmark.txt
│   ├── 20241108.102621.871345_mock.fastqc.benchmark.txt
│   ├── 20241108.102621.871345_mock.krona.benchmark.txt
│   ├── 20241108.102621.871345mock.R_alpha_div.benchmark.txt
│   ├── 20241108.102621.871345mock.R_barplot.benchmark.txt
│   ├── 20241108.102621.871345mock.R_conversion.benchmark.txt
│   ├── 20241108.102621.871345_mock.V300091236_L01_100.trimming.benchmark.txt
│   └── mock.kraken2.benchmark.txt
├── err_out
│   ├── jobnm_11453585.err
│   └── jobnm_11453585.out
├── EXTRACTED_FASTA
│   ├── mock1.fastq
│   ├── mock.1.fastq.gz
│   ├── mock2.fastq
│   └── mock.2.fastq.gz
├── KRAKEN_ASSIGN
│   ├── mock.Kraken2.biom
│   ├── mock.kraken2.report
│   └── mock.kraken2.txt
├── KRONA_HTML
│   ├── mock.html
│   └── mock.krona
├── METAPHLAN4
│   ├── mock
│   │   ├── mock.unmapped_genefamilies.tsv
│   │   ├── mock.unmapped_humann_temp
│   │   │   └── mock.unmapped.log
│   │   ├── mock.unmapped_pathabundance.tsv
│   │   └── mock.unmapped_pathcoverage.tsv
│   ├── mock.bz2
│   ├── mock_profiled.txt
│   ├── mocksam.bz2
│   └── mock.vsc.txt
├── R_ANALYSIS
│   ├── mock_alpha_div.csv
│   ├── mock_Barplot_phyla.jpeg
│   └── mock.rds
├── TRIMMOMATIC
│   ├── mock.1.fastq.gz
│   ├── mock.2.fastq.gz
│   ├── mock_qc
│   │   ├── mock.1_fastqc.html
│   │   ├── mock.1_fastqc.zip
│   │   ├── mock.2_fastqc.html
│   │   └── mock.2_fastqc.zip
│   ├── V300091236_L01_100_1_paired.fq.gz
│   ├── V300091236_L01_100_1_unpaired.fq.gz
│   ├── V300091236_L01_100_2_paired.fq.gz
│   └── V300091236_L01_100_2_unpaired.fq.gz
└── WGS_logs
    ├── 20241108.102621.871345.biom.log
    ├── 20241108.102621.871345.concatenation.log
    ├── 20241108.102621.871345.extracted_reads.log
    ├── 20241108.102621.871345.fastqc.log
    ├── 20241108.102621.871345.illumina_kraken2.log
    ├── 20241108.102621.871345.krona.log
    ├── 20241108.102621.871345_mock.alignment.log
    ├── 20241108.102621.871345_mock.all.rule.log
    ├── 20241108.102621.871345_mock.bam2fasta.log
    ├── 20241108.102621.871345_mock.index_bam.log
    ├── 20241108.102621.871345.R_alpha_div.log
    ├── 20241108.102621.871345.R_barplot.log
    ├── 20241108.102621.871345.R_phyloseq.log
    ├── 20241108.102621.871345.V300091236_L01_100.trimmomatic.log
    └── 20241108.102909.408489.illumina_kraken2.log
```

---

## Benchmark

At the end of the procedure you might want to know the resources that were used. A folder called "Benchmark" will be created containing for each of some of the rules the following parameters:

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
To debug the **meTAline** pipeline efficiently, use the `--bind` flag in Singularity to mount your modified source code into the container without rebuilding the image.

### Accessing the meTAline Source Code in the Image
The **meTAline** source code inside the Singularity image (`metaline.sif`) is located at `/meTAline`.

You can verify this with:
*E.G.:*
```bash
[user@host workspace]$ singularity run --cleanenv <path/to/metaline.sif> ls /meTAline
Helper_scripts_MN5  Illumina_MGI_adapters.fa  README.md  adapter_list_new.txt  conda_envs  debug.greasy.job  lib  meTAline.smk
[user@host workspace]$
```

### Overwriting Source Code with --bind
To test modifications to the pipeline, replace the source code in the image with your local version during execution.
Use the --bind flag to mount your local repository into the container, overwriting the existing directory.

#### Example
For this example you have the current working directory is inside the cloned meTAline repository, which looks like this:

*Note that the meTAline snakemake source code, ./meTAline, is right next to the built ./metaline.sif file!*

```bash
[user@host workspace]$ ls -lisah <path/to/cloned/meTAline/repository>
total 2.7G
75807561554 1.0K drwxrwxr-x  9 user bsc40 4.0K Dec 19 14:27 .
 6980952720 1.0K drwxr-sr-x 25 user bsc40 4.0K Dec 19 10:43 ..
75825726583 1.0K drwxrwxr-x  8 user bsc   4.0K Dec 11 11:22 .git
75825726584 1.0K drwxrwxr-x  3 user bsc   4.0K Dec 11 11:22 .github
75825726590  512 -rw-rw-r--  1 user bsc   3.2K Dec 11 11:22 .gitignore
75825726591  16K -rw-rw-r--  1 user bsc   5.8K Dec 11 11:22 Dockerfile
75825726592  48K -rw-rw-r--  1 user bsc    35K Dec 11 11:22 LICENSE
75825726595  512 -rw-rw-r--  1 user bsc    798 Dec 11 11:22 Makefile
75825726597  32K -rwxrwxr-x  1 user bsc    17K Dec 11 11:22 README.md
75825726599  512 -rw-rw-r--  1 user bsc   1.7K Dec 11 11:22 example_hpc_sbatch.job
75825726600  512 -rw-rw-r--  1 user bsc   1.7K Dec 11 11:22 example_job_def_mn5.job
75825726585 1.0K drwxrwxr-x  3 user bsc   4.0K Dec 11 11:22 external-sources
75825726586 1.0K drwxrwxr-x  5 user bsc   4.0K Dec 19 10:41 meTAline
75825726603  16K -rw-rw-r--  1 user bsc   5.2K Dec 11 11:54 metaline-singularity.def
75810056215 2.7G -rwxr-xr-x  1 user bsc   2.7G Dec 19 11:00 metaline.sif
75825726588 1.0K drwxrwxr-x  2 user bsc   4.0K Dec 11 11:22 scripts
75825726589 1.0K drwxrwxr-x  2 user bsc   4.0K Dec 11 11:22 test_input
[user@host workspace]$
```

Modify the command as follows:

**Original Command:**
```bash
singularity run --cleanenv metaline.sif \
    metaline \
    -r all \
    -j 16 \
    --configfile my_metaline_config.json
```

**Modified command:**
```bash
singularity run --cleanenv --bind ./meTAline:/meTAline metaline.sif \
    metaline \
    -r all \
    -j 16 \
    --configfile my_metaline_config.json
```

This allows you to test modifications to the pipeline without needing to rebuild or recompile the Singularity image, significantly speeding up development and debugging iterations.

#### Debugging tools
While IDE debuggers may not always work with Singularity images, most of the pipeline uses Python scripts, allowing you to utilize pdb (Python Debugger) for interactive debugging sessions.

To set a breakpoint in a Python script, use the **breakpoint()** function:

```bash
x = 10
y = 'Hi'
z = 'Hello'
print(y)

breakpoint()

print(z)
```

For further information about pdb refer to the [official documentation](https://docs.python.org/3/library/pdb.html).

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
