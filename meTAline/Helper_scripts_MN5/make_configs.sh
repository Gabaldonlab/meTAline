#!/bin/bash

#Raw data directory:
FITreaddir="/gpfs/projects/bsc40/current/okhannous/AECC_MGI_IRB/DATA/MGI_test_run/L02"

#Script to generate the config file, that is find in the shared pipeline:
config="/gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-0.8.0-alpha/meTAline/lib/config/generate_config.py"

#Available human reference genomes:
T2T="/gpfs/projects/bsc40/current/asukhorukova/index/T2T/T2T"
grch38="/gpfs/projects/bsc40/current/asukhorukova/index/grch38/genome"
pan="/gpfs/projects/bsc40/current/asukhorukova/index/pan/pan"

#Kraken2 standard database
kraken="/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE"
kmer="/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE/database150mers.kmer_distrib"

#Kraken2 broad gut microbiome database
bgut="/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/Broad_gut_microbiome_db"

#Base directory
basedir="/gpfs/projects/bsc40/current/okhannous/AECC_MGI_IRB/Preprocessing/OUT/L02"

#Creation of a folder to deposit the config files
mkdir $basedir/configs
#If there was an already created config database: make sure that there are no old config files there: these will be removed
rm $basedir/configs/* 

#For each file that contains the particular regular expression, sort them and take the particular number of interest
for f in $(ls $FITreaddir/*_50ng_1.fq.gz | sort -u )
do
    suffix="_50ng_1.fq.gz"
    f1=${f%"$suffix"}
    fo=${f1#"$FITreaddir"/}
    #The fastqs wildcards
    foo=${fo}_50ng
  
    #Running the config creation python script. Leaving default taxid, output directories, log directory, trimmomatic parameters and using the same number of cores for all tools

    singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-0.8.0-alpha/metaline.sif python3 /meTAline/lib/config/generate_config.py \
    --configFile $basedir/configs/config.$fo.json \
    --extension fq.gz \
    --basedir $basedir/$fo \
    --reads-directory $FITreaddir \
    --reference-genome $T2T \
    --krakendb $bgut \
    --sample-barcode "$fo" \
    --fastq-prefix "$foo" \
    --metaphlan_db /gpfs/projects/bsc40/current/okhannous/Metaphlan4/db \
    --metaphlan_Index mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n_db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/chocophlan \
    --protein_db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/uniref

done

rm $basedir/configs/config.*..json #cleanup
rm $basedir/configs/config.*.*1.f*q.gz.json #cleanup
