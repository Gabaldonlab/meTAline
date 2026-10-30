#!/bin/bash

set -e

singularity run --cleanenv /path_to_metaline_directory/metaline.sif metaline-prepare-greasy-array-job \
    --metaline-rule "all" \
    --basedir ./test_output \
    --metaline-cmd "module load singularity && singularity run --cleanenv /path_to_metaline_directory/metaline.sif metaline" \
    --greasy-cmd "path_to_greasy_dir/greasy" \
    --reference-genome path_to_reference_directory \
    --krakendb path_to_kraken_db_directory \
    --reads-directory ./input_data \
    --fastq-extension .fq.gz \
    --metaphlan-db path_to_metaphlan_db_directory \
    --metaphlan-index  mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db path_to_metaphlan_dir/chocophlan \
    --protein-db path_to_metaphlan_dir/metaPhla-db/uniref \
    --sylph_db path_to_sylph_db/sylph_db.syldb \
    --sylph_taxonomy path_to_sylph_taxonomy/taxonomy_sylph.tsv \
    --joblist-size 2 `# Max. number of jobs per each node in the Greasy array.` \
    --cluster-node-cores 24 \
     --trimmo-cores 16 \
     --hisat2-cores 16 \
     --kraken2-cores 16 \
     --snakemake-cores 16
