#!/bin/bash

set -e

singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.sif metaline-prepare-greasy-array-job \
    --metaline-rule "all" \
    --basedir ./test_output \
    --metaline_cmd "module load singularity && singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.sif metaline" \
    --greasy_cmd "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy" \
    --reference-genome /gpfs/projects/bsc40/project/pipelines/WGS/reference_genomes/index/T2T/T2T \
    --krakendb /gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE \
    --reads-directory /gpfs/projects/bsc40/current/okhannous/MeTAline_paper/BMC_version/raw_data \
    --fastq_extension fq.gz \
    --metaphlan-db /gpfs/projects/bsc40/current/okhannous/Metaphlan4/db \
    --metaphlan-index  mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/chocophlan \
    --protein-db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/uniref \
    --joblist_size 2 `# Max. number of jobs per each node in the Greasy array.` \
    --cluster-node-cores 110 `# Adjust these cores according to your cluster node's resouces.` \
    --trimmo-cores 90 `# Adjust these cores according to the --cluster-node-cores.` \
    --hisat2-cores 90 `# Adjust these cores according to the --cluster-node-cores.` \
    --kraken2-cores 90 `# Adjust these cores according to the --cluster-node-cores.` \
    --snakemake-cores 16 `# Adjust these cores according to the --cluster-node-cores.`
