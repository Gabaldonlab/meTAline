## Base Docker command

```sh
docker run --rm \  # Runs the Docker container and removes it after it finishes running
    -v /path/to/input:/input-directory \  # Mounts the host's input directory to the /app/input directory inside the container
    -v /path/to/output:/output-directory \  # Mounts the host's output directory to the /app/output directory inside the container
    cgenomics/metaline:latest \  # The name of your Docker image
    <command_placeholder> /input-directory /output-directory  # The script is executed with /input-directory and /output-directory as arguments
```

Breakdown of the command:
`docker run`: Runs the Docker container.
`--rm`: Removes the container after it finishes running.
`-v /path/to/input-directory:/input-directory`: Mounts the host's input directory to the /input-directory directory inside the container.
`-v /path/to/output-directory:/output-directory`: Mounts the host's output directory to the /output-directory directory inside the container.
`your-docker-image-name`: The name of your Docker image.
`<command_placeholder> /input-directory /output-directory`: The script is executed with /input-directory and /output-directory as arguments. These are the mounted directories inside the container.

Explanation:
`/input-directory`: Inside the container, it refers to the directory mounted from /path/to/input on the host.
`/output-directory`: Inside the container, it refers to the directory mounted from /path/to/output on the host.

Your script will now be able to access the input and output directories using the arguments provided.

## Creation of the config file

```sh
docker run --rm `# Runs the Docker container and removes it after it finishes running` \
    -v ./test_input:/input-directory `# Mounts the host's input directory to the /input-directory directory inside the container` \
    -v ./test_output:/output-directory `# Mounts the host's output directory to the /output-directory directory inside the container` \
    cgenomics/metaline:latest `# The name of your Docker image` \
    metaline-generate-config `# Command to generate the config file` \
        --config-file /output-directory/test_run.json `# Name you want to give to the output configuration file` \
        --extension fq.gz `# Extension of the raw reads files` \
        --basedir /output-directory/ `# Base directory in which the pipeline will create output` \
        --reads-directory  /input-directory/  `# Directory where the fastqs reads are stored` \
        --reference-genome /input-directory/test_datasets/grch38_index/genome  `# Reference genome to make host depletion (optional)` \
        --krakendb /input-directory/test_datasets/minikraken2_v2_8GB_201904_UPDATE  `# Kraken2 database` \
        --sample-barcode test_run `# Identification of the sample in the output files` \
        --fastq-prefix V300091236_L01_100 `# Wildcard parameter corresponding to the basename of the fastq files located in the --reads-directory` \
        --metaphlan-db /input-directory/test_datasets/metaphlan_dbs `# Path where the metaphlan4 database is located` \
        --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 `# Index of the metaphlan4 database` \
        --n-db  /input-directory/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz `# Humann database to do the nucleotide search (based on already built annotations)` \
        --protein-db /input-directory/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz  `# Humann database to do the translation search (by default this is by-passed)`
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
docker run --rm `# Runs the Docker container and removes it after it finishes running` \
    -v ./test_input:/input-directory `# Mounts the host's input directory to the /input-directory directory inside the container` \
    -v ./test_output:/output-directory `# Mounts the host's output directory to the /output-directory directory inside the container` \
    cgenomics/metaline:latest `# The name of your Docker image` \
    metaline `# Command to run the pipeline` \
        kmer_taxonomy `# Rule to execute` \
        -j 16 `# Number of threads provided to the pipeline` \
        --configfile /output-directory/test_run.json  `# Indicate where the created config file (path in --config-file) is`
```

**OR test it with dry run only!**

```bash
docker run --rm \  # Runs the Docker container and removes it after it finishes running
    -v ./test_input:/input-directory \  # Mounts the host's input directory to the /input-directory directory inside the container
    -v ./test_output:/output-directory \  # Mounts the host's output directory to the /output-directory directory inside the container
    cgenomics/metaline:latest \  # The name of your Docker image
    metaline \  # Command to run the pipeline
        kmer_taxonomy \  # Rule to execute
        -j 16 \  # Number of threads provided to the pipeline
        --configfile /output-directory/test_run.json \  # Indicate where the created config file (path in --config-file) is
        -n  # Dry run
```

The `-n` option is for a dry run, avoiding the execution of any job yet. `-j` is the number of threads provided to the pipeline. In this specific example, we selected the rule all using the parameter `-r`. You can use different configs by adding the `--configfile` parameter but it requires that you comment the configfile argument in the main meTAline.smk.

For further information about the flags, run:

**metaline-generate-config**
```bash
docker run --rm cgenomics/metaline:latest metaline-generate-config -h
```

**metaline**
```bash
docker run --rm cgenomics/metaline:latest metaline -h
```
