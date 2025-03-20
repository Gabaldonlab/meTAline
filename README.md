[![DOI](https://zenodo.org/badge/431438117.svg)](https://zenodo.org/badge/latestdoi/431438117)

# MeTAline: a snakemake pipeline for the study of metagenomes <a id="metaline" />

MeTAline, is a snakemake pipeline for metagenomics analysis. MeTAline, facilitates an efficient workflow to preprocess short reads metagenomics data: from read trimming and filtering, through host read subtraction to taxonomic classification using both k-mer and gene marker-based approaches, and functional profiling of the samples.

<div align="center">
  <img src="https://github.com/user-attachments/assets/1c248a72-625f-480a-97ad-9022714a3dbc" width="60%" height="60%">
</div>

---

# Table of content

- [Running in HPC environment](#running-in-hpc)
- [Test sample and output directory example](#test-sample-and-output-directory-example)
- [Benchmark](#benchmark)
- [Setup Singularity image for debugging](#singularity-debugging)
    - [Accessing the meTAline Source Code in the Image](#enter-into-image)
    - [Overwriting Source Code with --bind](#overwrite-with-bind-flag)
        - [Example](#example)
    - [Debugging tools](#debugging-tools)
- [4. Citations / Acknowledgments](#citations)
    - [4.1. Snakemake](#snakemake)
    - [4.2. Kraken2](#kraken2)
    - [4.3. KrakenTools](#kraken-tools)
    - [4.4. GNU Parallel](#gnu-parallel)
    - [4.5. Trimmomatic](#trimmomatic)
    - [4.6 Bracken](#bracken)
    - [4.7 Krona](#krona)
    - [4.8 HTSlib](#hts-lib)
    - [4.9 Samtools](#samtools)
    - [4.10 Hisat2](#hisat2)
    - [4.11 MetaPhlAn](#metaphlan)
    - [4.12 Humann](#humann)
    - [4.13 FastQC](#fastqc)

---


---

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

## Test sample and output directory example <a id="test-sample-and-output-directory-example" />

To try the pipeline we provide test mock samples and a selection of external datasets to download:

`test_input` directory:
```sh
tree -L 1 ./test_input/
./test_input/
├── explanation.txt
├── V300091236_L01_100_1.fq.gz
└── V300091236_L01_100_2.fq.gz
```

Download the `test_datasets` directory's datasets by running the following command:
```sh
./scripts/download_test_datasets
```

```sh
tree -L 1 ./test_data
test_data
├── chocophlan.v4_alpha.tar.gz
├── grch38_genome.tar.gz
├── minikraken2_v2_8GB_201904.tgz
├── mpa_vJun23_CHOCOPhlAnSGB_202307_marker_info.txt.bz2
├── mpa_vJun23_CHOCOPhlAnSGB_202307.md5
├── mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt.bz2
├── mpa_vJun23_CHOCOPhlAnSGB_202307.tar
└── uniref90_annotated_v4_alpha_ec_filtered.tar.gz
```

An output sample of the pipeline, after running all the rules:
```sh
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

## Benchmark <a id="benchmark" />

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

## Setup Singularity image for debugging <a id="singularity-debugging" />
To debug the **meTAline** pipeline efficiently, use the `--bind` flag in Singularity to mount your modified source code into the container without rebuilding the image.

### Accessing the meTAline Source Code in the Image <a id="enter-into-image" />
The **meTAline** source code inside the Singularity image (`metaline.sif`) is located at `/meTAline`.

You can verify this with:
```bash
[user@host workspace]$ singularity run --cleanenv <path/to/metaline.sif> ls /meTAline
Helper_scripts_MN5  Illumina_MGI_adapters.fa  README.md  adapter_list_new.txt  conda_envs  debug.greasy.job  lib  meTAline.smk
[user@host workspace]$
```

### Overwriting Source Code with --bind <a id="overwrite-with-bind-flag" />
To test modifications to the pipeline, replace the source code in the image with your local version during execution.
Use the --bind flag to mount your local repository into the container, overwriting the existing directory.

#### Example <a id="example" />
For this example you have the current working directory inside the cloned meTAline repository, which looks like this:

*Note that the meTAline snakemake directory, ./meTAline, is right next to the built ./metaline.sif file!*

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

### Debugging tools <a id="debugging-tools" />
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

For further information about pdb refer to the [official Python 3 documentation](https://docs.python.org/3/library/pdb.html).

---

## 4. Citations / Acknowledgments <a id="citations" />

### 4.1. Snakemake <a id="snakemake" />

[Sustainable data analysis with Snakemake](<(https://doi.org/10.12688/f1000research.29032.1)>)

### 4.2. Kraken2 <a id="kraken2" />

(Improved metagenomic analysis with Kraken 2)[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0]

### 4.3. KrakenTools <a id="kraken-tools" />

[KrakenTools - Github](https://github.com/jenniferlu717/KrakenTools)

[Lu J, Rincon N, Wood D E, Breitwieser F P, Pockrandt C, Langmead B, Salzberg S L, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)]
(https://www.nature.com/articles/s41596-022-00738-y)

Relevant paper for usage of KrakenTools:

1. [Kraken 1](https://github.com/DerrickWood/kraken)
2. [Kraken 2](https://github.com/DerrickWood/kraken2)
3. [KrakenUniq](https://github.com/fbreitwieser/krakenuniq)
4. [Bracken](https://github.com/jenniferlu717/Bracken)

### 4.4. GNU Parallel <a id="gnu-parallel" />

[Tange, O. (2021, August 22). GNU Parallel 20210822 ('Kabul'). Zenodo.](https://doi.org/10.5281/zenodo.5233953)

### 4.5. Trimmomatic <a id="trimmomatic" />

[Trimmomatic: a flexible trimmer for Illumina sequence data](https://doi.org/10.1093/bioinformatics/btu170)

### 4.6 Bracken <a id="bracken" />

[Bracken: estimating species abundance in metagenomics data](https://peerj.com/articles/cs-104/)

### 4.7 Krona <a id="krona" />

[Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.](http://www.ncbi.nlm.nih.gov/pubmed/21961884)

### 4.8 HTSlib <a id="hts-lib" />

[HTSlib: C library for reading/writing high-throughput sequencing data; James K Bonfield, John Marshall, Petr Danecek, Heng Li, Valeriu Ohan, Andrew Whitwham, Thomas Keane, Robert M Davies GigaScience, Volume 10, Issue 2, February 2021, giab007](https://doi.org/10.1093/gigascience/giab007)

### 4.9 Samtools <a id="samtools" />

[Twelve years of SAMtools and BCFtools; Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li GigaScience, Volume 10, Issue 2, February 2021, giab008](https://doi.org/10.1093/gigascience/giab008)

### 4.10 Hisat2 <a id="hisat2" />

[Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://www.nature.com/articles/s41587-019-0201-4)

### 4.11 MetaPhlAn <a id="metaphlan" />

[Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4.](https://doi.org/10.1038/s41587-023-01688-w)

### 4.12 Humann <a id="humann" />

[Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3](https://doi.org/10.7554/eLife.65088)

### 4.13 FastQC <a id="fastqc" />

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
