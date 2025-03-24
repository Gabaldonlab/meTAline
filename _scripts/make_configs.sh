#!/bin/bash

#Raw data directory:
Reads_dir="/gpfs/projects/bsc40/current/okhannous/MeTAline_paper/BMC_version/raw_data" # DONE

#Script to generate the config file, that is find in the shared pipeline:
config="/gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-0.8.0-alpha/meTAline/lib/config/generate_config.py"  # DONE

#Available human reference genomes:
T2T="/gpfs/projects/bsc40/project/pipelines/WGS/reference_genomes/index/T2T/T2T"  # DONE
grch38="/gpfs/projects/bsc40/current/asukhorukova/index/grch38/genome"  # DONE
pan="/gpfs/projects/bsc40/current/asukhorukova/index/pan/pan"  # DONE

#Kraken2 standard database
kraken="/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE"  # DONE

#Kraken2 broad gut microbiome database
bgut="/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/Broad_gut_microbiome_db"  # DONE

#Base directory
basedir="/gpfs/projects/bsc40/current/okhannous/MeTAline_paper/BMC_version/raw_data/OUT"  # DONE

#Creation of a folder to deposit the config files
mkdir $basedir/configs  # DONE
#If there was an already created config database: make sure that there are no old config files there: these will be removed
rm $basedir/configs/*  # DONE

#For each file that contains the particular regular expression, sort them and take the particular number of interest
for f in $(ls $Reads_dir/*[[:punct:]]*1.f*q.gz | sort -u)
do

    fo=$(basename "$f")

    foo=${fo%[[:punct:]]*1.f*q.gz}

    #Running the config creation python script. Leaving default taxid, output directories, log directory, trimmomatic parameters and using the same number of cores for all tools

    singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-0.8.0-alpha/metaline.sif python3 /meTAline/lib/config/generate_config.py \
    --configFile $basedir/configs/config.$foo.json \
    --extension fq.gz \
    --basedir /$basedir/$foo \
    --reads-directory $Reads_dir \
    --reference-genome $T2T \
    --krakendb $kraken \
    --sample-barcode "$foo" \
    --fastq-prefix "$foo" \
    --metaphlan_db /gpfs/projects/bsc40/current/okhannous/Metaphlan4/db \
    --metaphlan_Index mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n_db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/chocophlan \
    --protein_db /gpfs/projects/bsc40/project/pipelines/WGS/metaPhlan/metaPhla-db/uniref

done

rm $basedir/configs/config.*..json #cleanup
rm $basedir/configs/config.*.*1.f*q.gz.json #cleanup
