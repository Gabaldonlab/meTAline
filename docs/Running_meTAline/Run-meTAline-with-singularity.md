## Running meTAline with singularity <a id="running-in-hpc" />

Before starting this step by step process you should have completed the build and setup process. Further information can be found at [Build meTAline singularity image](../Building_images/Build-meTAline-singularity-image.md). Alternatively You can download an already built singularity image [metaline.sif](https://github.com/Gabaldonlab/meTAline/releases/download/v1.0.0/metaline.sif)



> [!IMPORTANT]
> To avoid any unexpected results when finding the sample code prefixes from the fastq.gz files make sure that the `--reads-directory` argument's directory path and its content conforms the following critarias
> - 1. Use ONLY gzipped (with .gz extension) files. Do not `gunzip` them before use!
>
> - 2. FORWARD and REVERSE files of a sample must be in pairs
> ```sh
> # E.G.:
> E100043441_L01_1_1.1.fq.gz
> E100043441_L01_1_1.2.fq.gz
>
> E100043441_L01_1_2.1.fq.gz
> E100043441_L01_1_2.2.fq.gz
>
> # See: The FORWARD is split into 2 files, and the REVERSE also. (N of FORWARD files == N of REVERSE files)
> ```
>
> - 3. Nested directories are allowed, but the FORWARD and REVERSE sequences of a sample must be in the same subdirectory:
> ```sh
> ── TG0002
>    ├── E100043441_L01_1_1.fq.gz
>    └── E100043441_L01_1_2.fq.gz
> ── TG0003
>    ├── E100043441_L01_2_1.fq.gz
>    └── E100043441_L01_2_2.fq.gz
> ── TG0004
>    ├── E100043441_L01_3_1.fq.gz
>    └── E100043441_L01_3_2.fq.gz
> ```
> **OR in case of flat directory**
> ```sh                                                                                                                  │~
> ├── V300091236_L01_100_1.fq.gz
> ├── V300091236_L01_100_2.fq.gz
> ```

**Step by step:**

## 1. Download samples (or use yours):
```sh
cd ./test_input

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/003/ERR1018263/ERR1018263_2.fastq.gz  # ERR1018263 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/003/ERR1018263/ERR1018263_1.fastq.gz  # ERR1018263 read 2

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/004/ERR1018234/ERR1018234_2.fastq.gz  # ERR1018234 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/004/ERR1018234/ERR1018234_1.fastq.gz  # ERR1018234 read 2

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/002/ERR1018192/ERR1018192_1.fastq.gz  # ERR1018192 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/002/ERR1018192/ERR1018192_2.fastq.gz  # ERR1018192 read 2

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1018218/ERR1018218_1.fastq.gz  # ERR1018218 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1018218/ERR1018218_2.fastq.gz  # ERR1018218 read 2
```

## 2. Download test datasets (or use yours):
```sh
make download-test-datasets
```


##  3. Creation of the config file and running of the pipeline
 There are several options at this point:

### If using a hpc cluster:
 Continue using [Run in a HPC environment](Run-meTAline-in-HPC-environment.md) instructions

### If running locally with parallelisation:
 Continue using [Parallelize on local machine](Parallelize-meTAline-on-local-machine.md) instructions

### If running locally without parallelisation:
#### Create the config

```sh
singularity run --cleanenv ./metaline.sif metaline-generate-config `# Path where the downloaded .sif image is.` \
    --config-file ./test_output/test_run.json `# Name you want to give to the output configuration file.` \
    --fastq-extension .fq.gz `# Extension of the raw reads files` \
    --basedir ./test_output `# Base directory in which the pipeline will create output` \
    --reads-directory  ./test_input `# Directory where the fastqs reads are stored` \
    --reference-genome ./test_input/test_datasets/grch38_index/genome `# Reference genome to make host depletion (optional)` \
    --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE `# Kraken2 databse` \
    --sample-barcode test_run `# Identification of the sample in the output files` \
    --fastq-prefix V300091236_L01_100 `# Wildcard parameter corresponding to the basename of the fastq files` \
    --metaphlan-db ./test_input/test_datasets/metaphlan_dbs `# Path where the metaphlan4 database is located` \
    --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 `# Index of the metaphlan4 database` \
    --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz `# Humann database to do the nucleotide search (based on already built annotations)` \
    --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz  `# Humann database to do the translation search (by default this is by-passed)`
```
You can change the required and default parameters your job using the [Parameters documentation](../Available-target-rules.md).

#### Create the required folders and file
```sh
mkdir -p ./Benchmark
mkdir -p ./Benchmark/plots
[ -f ./Benchmark/required.txt ] || touch ./Benchmark/required.txt
```

#### Run the pipeline, indicating the config file with all the details
```sh
singularity run --cleanenv ./metaline.sif metaline \
    kmer_taxonomy \
    -j 16 \
    --configfile ./test_output/test_run.json  # Indicate where the created config file is
```

**OR test it with dry run only!**

```bash
singularity run --cleanenv metaline.sif \
    metaline \
    all \
    -j 16 \
    --configfile my_metaline_config.json \
    -n # Dry run
```

The `-n` option is for a dry run, avoding the execution of any job yet. `-j` is the number of threads provided to the pipeline. In this specific example, we selected the rule all using the parameter `-r`. You can use different configs by adding the `--configfile` parameter but it requires that you comment the configfile argument in the main meTAline.smk.

For further information about the flags, run:

**metaline-generate-config**
```bash
singularity run --cleanenv metaline.sif metaline-generate-config -h
```

**metaline**
```bash
singularity run --cleanenv metaline.sif metaline -h
```
