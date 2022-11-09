# meTAline: metagenomics Taxonomic Assignation pipeline
meTAline is taxonomic assignation pipeline implemented in SnakeMake for WGS metagenomics data.

Additionally, you can use the bash wrapper located here: https://github.com/Gabaldonlab/meTAline/tree/bash_wrapper
## Getting Started
First it is required to pre-install conda with either anaconda or miniconda3. 

**anaconda** - please follow these [instructions](https://docs.anaconda.com/anaconda/install/)

Alternatively, you can install:

**miniconda3** - please follow these [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Pipeline Installation
After installing conda, download and install the pipeline

```Shell
# First clone the meTAline repository
git clone https://github.com/Dfupa/meTAline.git
# Change to meTAline/conda_envs directory
cd meTAline/conda_envs
# Set the conda environment with all the necessary dependencies
conda env create -f meTAline_env.yml #The main environment 
conda env create -f meTAline_Rrule.yml #You should also create this environment to be activated and used for the R rule
# Activate the conda environment
conda activate meTAline
```

# Input
This pipline uses as an input a config file located in the meTAline/lib/config/ directory.  The user can find in that directory a script called generate_config.py which will build a config file with the parameters provided by the user.

More information is provided by using the help parameter of generate_config.py, provide at least the mandatory arguments (the host reference genome index can be omitted if you are running an environmental sample).

```Shell
# Always assuming the environment is active
python3 lib/config/generate_config.py -h
# An example of the real command
python3 lib/config/generate_config.py --sample-barcode metagenomics_sample --reads-directory /path/to/the/reads/directory/ --reference-genome /path/to/reference_genomes/human_hg38/HUMAN_index --krakendb /path/to/KRAKEN2_DB_COMPLETE --kmer_dist /path/to/KRAKEN2_DB/KRAKEN2_DB_COMPLETE/database150mers.kmer_distrib --basedir /my/desired/directory/ --configFile config.json

```

# Target rules
The target rules currently available to use are:

- rule **all**: Outputs the [Kraken2](https://github.com/DerrickWood/kraken2) taxonomic assignation in [Krona](https://github.com/marbl/Krona) format, as well as extracting the desired reads.
- rule **krona_and_reads**: This rule is used if you already have the kraken2 taxonomic assignation and you want to generate the krona images and extract reads
- rule **taxonomy_assignation:**: This rule is used if you already have the concatenated/unmapped reads and you want to perform the taxonomic assignation¡ without running the initial steps of the pipeline
- rule **biom_format**: This rule is used if you already have the Kraken2 and bracken assignment and you want to pass from their reports to biom format for further analysis
- rule **R_Analysis**: This rule is used to make a basic R analysis and plotting alpha diversity.

# Benchmark

At the end of the procedure you might want to know the resources that were used. A folder called "Benchmark" will be created containing for each of the results the follwing parameters: 

**s**	--> Running time in seconds

**h: m: s**	--> Running time in hours, minutes and seconds

**max_rss**	--> "Maximum Resident Set Size". Total amount of physical memory.

**max_vms** --> "Maximum Virtual Memory Size". Total amount of virtual memory.

**max_uss** --> "Unique Set Size".	Memory that is unique to a process. 

**max_pss** --> "Proportional Set Size". Amount of memory shared with other processes.

**io_in** --> MB read.

**io_out** --> MB written.	

**mean_load** --> CPU usage over time / Total running time.


# How to run snakemake
In order to run the meTAline pipeline just use the following command

```Shell
# Always assuming the environment is active
snakemake -s meTAline.smk -r taxonomy_assignation -j 16 -n

#If you want to run it with thour ouwn configs (after commenting the config parameter in the meTAline.smk file)
snakemake -s meTAline.smk -r taxonomy_assignation -j 16 -n --configfile myconfig.json
```
The -n option is for a dry run, avoding the execution of any job yet. -j is the number of threads provided to the pipeline.  In this specific example, we selected the rule all using the parameter -r. You can use different configs by adding the --configfile parameter but it requires that you comment the configfile argument in the main meTAline.smk.

More information regarding Snakemake and its commands can be found through Snakemake [documentation](https://snakemake.readthedocs.io/en/stable/index.html).


