# meTAgen: metagenomics Taxonomic Assignation pipeline
meTAgen is taxonomic assignation pipeline implemented in SnakeMake for WGS metagenomics data.

## Getting Started
First it is required to pre-install conda with either anaconda or miniconda3. 

**anaconda** - please follow these [instructions](https://docs.anaconda.com/anaconda/install/)

Alternatively, you can install:

**miniconda3** - please follow these [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Pipeline Installation
After installing conda, download and install the pipeline

```Shell
# First clone the meTAgen repository
$ git clone https://github.com/Dfupa/meTAgen.git
# Change to meTAgen directory
$ cd meTAgen/
# Set the conda environment with all the necessary dependencies
$ conda env create -f meTAgen_env.yml
# Activate the conda environment
$ conda activate meTAgen
```

# Input
