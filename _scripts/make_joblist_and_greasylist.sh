#!/bin/bash

#Base directory
basedir="/gpfs/projects/bsc40/current/okhannous/MeTAline_paper/BMC_version/raw_data/OUT"
#Config files directory
configdir=$basedir/configs
#Path in which we can find the main snakemake file
smk="/gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-0.8.0-alpha/meTAline/meTAline.smk"
#Raw data directory
datadir="/gpfs/projects/bsc40/current/okhannous/MeTAline_paper/BMC_version/raw_data"
cores=24 # make sure that this is the same value as in the make_configs script

rm $basedir/joblist.txt

for file in $(ls $configdir | sort -u)
do
    # extract sample name
    f=$(basename "$file")
    fo=${f#"config."}
    foo=${fo%%.*} 

   # echo $foo

    

    mkdir $foo
    echo "module load singularity && singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-0.8.0-alpha/metaline.sif metaline -r all -j 16 --configfile $configdir/$file" >> $basedir/joblist.txt

done

# make smaller joblists
for i in {1..2}
do
    head -$(($i*2)) $basedir/joblist.txt | tail -2 > $basedir/joblist$i.txt
done

# make greasy list
echo "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy $basedir/joblist1.txt" > $basedir/list_greasy.txt
for i in {2..2}
do
    echo "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy $basedir/joblist$i.txt" >> $basedir/list_greasy.txt
done
