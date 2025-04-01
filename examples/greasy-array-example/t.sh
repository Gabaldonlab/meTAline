#!/bin/bash

set -e

singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.fix.sif metaline-prepare-greasy-array-job \
    --basedir /gpfs/projects/bsc40/current/dmajer/metaline-prepare-greasy-array-job/OUT \
    --metaline_cmd "module load singularity && singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.fix.sif metaline -r all -j 16 --configfile" \
    --greasy_cmd "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy" \
    --reference-genome /gpfs/projects/bsc40/project/pipelines/WGS/reference_genomes/index/T2T/T2T \
    --krakendb /gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE \
    --reads-directory /gpfs/projects/bsc40/current/okhannous/MeTAline_paper/BMC_version/raw_data \
    --fastq_extension fq.gz \
    --metaphlan-db /gpfs/projects/bsc40/current/okhannous/Metaphlan4/db \
    --metaphlan-index  mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/chocophlan \
    --protein-db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/uniref \
    --max_workers 4 \
    --joblist_size 2