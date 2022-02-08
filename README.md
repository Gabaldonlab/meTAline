# meTAline: metagenomics Taxonomic Assignation pipeline
meTAline is taxonomic assignation pipeline implemented in SnakeMake for WGS metagenomics data.

## Getting Started
First it is required to pre-install conda with either anaconda or miniconda3. 

**anaconda** - please follow these [instructions](https://docs.anaconda.com/anaconda/install/)

Alternatively, you can install:

**miniconda3** - please follow these [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Pipeline Installation
After installing conda, download and install the pipeline

```Shell
# First clone the meTAgen repository
$ git clone https://github.com/Dfupa/meTAline.git
# Change to meTAline directory
$ cd meTAline/
# Set the conda environment with all the necessary dependencies
$ conda env create -f meTAline_env.yml
# Activate the conda environment
$ conda activate meTAline
```

# Input
This pipline uses as an input a config file located in the meTAline/lib/config/ directory.  The user can find in that directoy a script called get_config.py which will build a config file with the parameters provided by the user.

More information is provided by using the help parameter of generate_config.py

```Shell
# Always assuming the environment is active
$ python3 lib/config/generate_config.py -h
```

# Target rules
The target rules currently available to use are:

- rule **all**: Outputs the [Kraken2](https://github.com/DerrickWood/kraken2) taxonomic assignation in [Krona](https://github.com/marbl/Krona) format, as well as extracting the desired reads.
- rule **krona_and_reads**: This rule is used if you already have the kraken2 taxonomic assignation and you want to generate the krona images and extract reads
- rule **taxonomy_assignation:**: This rule is used if you already have the concatenated/unmapped reads and you want to perform the taxonomic assignation¡ without running the initial steps of the pipeline



# How to run snakemake
In order to run the meTAline pipeline just use the following command

```Shell
# Always assuming the environment is active
$ snakemake -s meTAline.smk -r all -j 16 -n
```
The -n option is for a dry run, avoding the execution of any job yet. -j is the number of threads provided to the pipeline.  In this specific example, we selected the rule all using the parameter -r. 

More information regarding Snakemake and its commands can be found through Snakemake [documentation](https://snakemake.readthedocs.io/en/stable/index.html).
