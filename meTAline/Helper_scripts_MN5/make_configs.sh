#!/bin/bash

#Raw data directory:
R_dir="/gpfs/projects/bsc40/current/okhannous/Meta-analysis/DATA/PRJEB10878"

#Script to generate the config file, that is find in the shared pipeline:
config="/gpfs/projects/bsc40/project/pipelines/meTAline/meTAline_current/lib/config/generate_config.py"

#Available human reference genomes:
T2T="/gpfs/projects/bsc40/current/asukhorukova/index/T2T/T2T"
grch38="/gpfs/projects/bsc40/current/asukhorukova/index/grch38/genome"
pan="/gpfs/projects/bsc40/current/asukhorukova/index/pan/pan"

#Kraken2 database (this is the path for the standard database)
kraken="/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE"
kmer="/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE/database150mers.kmer_distrib"


#Creation of a folder to deposit the config files
mkdir $R_dir/configs
#If there was an already created config database: make sure that there are no old config files there: these will be removed
rm $R_dir/configs/* 

#For each file that contains the particular regular expression, sort them and take the particular number of interest
for f in $(ls $R_dir/*[[:punct:]]*1.f*q.gz | sort -u)
do

    fo=$(basename "$f")

    foo=${fo%[[:punct:]]*1.f*q.gz}

    #Running the config creation python script. Leaving default taxid, output directories, log directory, trimmomatic parameters and using the same number of cores for all tools
    python $config --configFile $R_dir/configs/config.$foo.json --extension fastq.gz --basedir $R_dir/$foo --reads-directory $R_dir --reference-genome $T2T --krakendb $kraken --kmer_dist $kmer --sample-barcode "$foo" --fastqs "$foo" 
done

rm $R_dir/configs/config.*..json #cleanup
rm $R_dir/configs/config.*.*1.f*q.gz.json #cleanup
