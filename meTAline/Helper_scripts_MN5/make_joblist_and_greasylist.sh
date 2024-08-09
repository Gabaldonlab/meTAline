#!/bin/bash

#Base directory
basedir="/gpfs/projects/bsc40/current/okhannous/Meta-analysis/DATA/PRJEB10878"
#Config files directory
configdir=$basedir/configs
#Path in which we can find the main snakemake file
smk="/gpfs/projects/bsc40/project/pipelines/meTAline/meTAline_current/meTAline.smk"
#Raw data directory
datadir="/gpfs/projects/bsc40/current/okhannous/Meta-analysis/DATA/PRJEB10878"
cores=24 # make sure that this is the same value as in the make_configs script

rm $basedir/joblist.txt

for file in $(ls $configdir | sort -u)
do
    # extract sample name
    f=$(basename "$file")
    fo=${f#"config."}
    foo=${fo%".json"}
   # echo $foo

    

    mkdir $foo
    echo "source /gpfs/projects/bsc40/project/pipelines/anaconda3/etc/profile.d/conda.sh && conda activate meTAline && snakemake -s $smk -r all --cores $cores --configfile $configdir/$file" >> $basedir/joblist.txt
done

# make smaller joblists, change this number according to the number of samples and the array distribution you do
for i in {1..64}
do
    head -$(($i*2)) $basedir/joblist.txt | tail -2 > $basedir/joblist$i.txt
done

# make greasy list
echo "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy $basedir/joblist1.txt" > $basedir/list_greasy.txt
for i in {2..64}
do
    echo "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy $basedir/joblist$i.txt" >> $basedir/list_greasy.txt
done
