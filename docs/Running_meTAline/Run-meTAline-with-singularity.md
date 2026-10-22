You can download an already built singularity image [metaline.sif](https://github.com/Gabaldonlab/meTAline/releases/download/v1.0.0/metaline.sif)

## Creation of the config file
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
Required and default parameters for this job can be found at (docs/Available-target-rules.md).

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
