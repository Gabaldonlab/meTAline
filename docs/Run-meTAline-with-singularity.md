You can download an already built singularity image [metaline.sif](https://github.com/Gabaldonlab/meTAline/releases/download/v1.0.0/metaline.sif)

## Creation of the config file
```sh
singularity run --cleanenv ./metaline.sif metaline-generate-config `# Path where the downloaded .sif image is.` \
    --config-file ./test_output/test_run.json `# Name you want to give to the output configuration file.` \
    --extension fq.gz `# Extension of the raw reads files` \
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

>[!IMPORTANT]
>Please note that while we designed the pipeline to be as generalizable as possible, we identified sample naming as a potential source of issues. To address this, we introduced the `--extension` parameter.We recommend setting the extension to either "fastq.gz" or "fq.gz". By default, meTAline uses "fastq.gz".Additionally, we added the `--prefix-fr` parameter, which allows you to specify the naming convention used to distinguish forward and reverse reads. By default the value will be "_" .
>
>**For example:**
>If _ is used before indicating 1 (for forward) and 2 (for reverse): {file}_1.{extension}
>
>If .R is used as the prefix: {file}.R1.{extension}
>
>**Examples:**
>For sample G12, you have these two files as raw data: **sample_G12_1.fastq.gz** and **sample_G12_2.fastq.gz**
>
>In this case you will indicate `--extension` fastq.gz and `--prefix-fr` "_", but in case the raw files are: "sample_G12.R1.fq.gz" and "sample_G12.R1.fq.gz", then:
>`--extension` fq.gz and `--prefix-fr ".R"`

## Run of the pipeline, indicating the config file with all the details
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
