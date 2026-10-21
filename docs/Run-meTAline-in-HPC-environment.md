## Running in HPC environment <a id="running-in-hpc" />

The MeTAline pipeline is intended to be used in HPC environment, because of its high RAM runtime requirement.
You can find an example job definition file for the SLURM schedular (usually this is used in HPCs) in **./example_hpc_sbatch.job** file.

In the following example SLURM .job file we use some Mare Nostrum 5 (Barcelon Supercomputing Center) specific additional configurations that you can ignore.
Also, the base commands are through the generated Singularity image, which first you need to build on your local machine (because it requires "sudo" permissions, which in most cases not available in HPC machines), then you can copy that built image into the remote HPC machine.

For further information about the build and setup process, visit our Wiki on [Build meTAline singularity image](https://github.com/Gabaldonlab/meTAline/wiki/Build-meTAline-singularity-image)

**EXAMPLE'S CONTENT**


*Later referenced as `example_hpc_sbatch.job`*
```bash
#!/bin/bash

#SBATCH --job-name=metaline_singularity_test # Job’s name (to trace in)
#SBATCH --output=/hpc-cluster-data/projects/my_project/current/metaline_testing/metaline-test-output/err_out/jobnm_%j.out # Output file, where
#SBATCH --error=/hpc-cluster-data/projects/my_project/current/metaline_testing/metaline-test-output/err_out/jobnm_%j.err # File where the error is written

#SBATCH --ntasks=1 # The number of parallel tasks
#SBATCH --cpus-per-task=110 # Number of CPUs per run task
#SBATCH --tasks-per-node=1 # The number of allocated task/node

#SBATCH --qos=your_quality_of_service # The queue for the job
#SBATCH --account=my_account
#SBATCH --time=24:00:00 # The time you request for the job
#SBATCH --constraint=highmem # To run in highmem nodes, in case it's available in your cluster.

module load singularity

singularity run --cleanenv ./metaline.sif metaline-generate-config \
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
    --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz  `# Humann database to do the translation search (by default this is by-passed)` \
    --trimmo-cores 90 \
    --hisat2-cores 90 \
    --kraken2-cores 90

singularity run --cleanenv metaline.sif metaline \
    all \
    -j 16 \
    --configfile test_kraken2_standard.json
```

>[!IMPORTANT]
>Please note that while we designed the pipeline to be as generalizable as possible, we identified sample naming as a potential source of issues. To address this, we introduced the `--extension` parameter.We recommend setting the extension to either "fastq.gz" or "fq.gz". By default, meTAline uses "fastq.gz".Additionally, we added the `--prefix-fr` parameter, which allows you to specify the naming convention used to distinguish forward and reverse reads.
>By default the value will be "_" .
>
>**For example:**
>If "_" is used before indicating 1 (for forward) and 2 (for reverse): {file}_1.{extension}
>
>If .R is used as the prefix: {file}.R1.{extension}
>
>**Examples:**
>For sample G12, you have these two files as raw data: "sample_G12_1.fastq.gz"  and "sample_G12_2.fastq.gz"
>
>In this case you will indicate `--extension` fastq.gz and `--prefix-fr "_"`, but in case the raw files are: **sample_G12.R1.fq.gz** and **sample_G12.R1.fq.gz**, then:
>`--extension` fq.gz and `--prefix-fr ".R"`

You can run it with:

```bash
sbatch ./example_hpc_sbatch.job
```

To see more SLURM job examples (using array of greasies for parallelization), you can visit the [examples directory](https://github.com/Gabaldonlab/meTAline/tree/main/examples) !
