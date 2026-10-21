## Overview
This script automates the preparation of a greasy array job for submitting multiple MeTAline jobs in a high-performance computing (HPC) cluster.

## Features
- Generates configuration files for MeTAline.
- Splits jobs into job lists for parallel execution.
- Creates a SLURM job submission script for execution.

## Installation
Clone the repository and navigate to the script directory:
```sh
git clone https://github.com/gabaldonlab/metaline
cp ./meTAline/lib/config/prepare_greasy_array_job.py /usr/local/bin/metaline-prepare-greasy-array-job
chmod +x /usr/local/bin/metaline-prepare-greasy-array-job
```

## Usage
Run the script with the required arguments as follows:

> [!IMPORTANT]
> The arguments for `--metaline-cmd` and `--greasy-cmd` must be quoted (E.G.: `--generate_config_cmd "singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.sif metaline-generate-config"`)

```sh
metaline-prepare-greasy-array-job \
    --basedir ./test_output \
    --metaline-cmd "module load singularity && singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.sif metaline" \
    --greasy-cmd "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy" \
    --reference-genome ./test_input/test_datasets/grch38_index/genome \
    --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE \
    --reads-directory  ./test_input \
    --fastq-extension fq.gz \
    --metaphlan-db ./test_input/test_datasets/metaphlan_dbs \
    --metaphlan-index  mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz \
    --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz \
    --max_workers 4 \
    --joblist-size 2
```

### Required Arguments
- `--basedir` → Base directory for output files.
- `--metaline-cmd` → Command to execute MeTAline.
- `--greasy-cmd` → Command to execute Greasy.
- `--reference-genome` → Path to the indexed reference genome.
- `--krakendb` → Path to the Kraken database.
- `--reads-directory` → Directory containing FASTQ read files.
- `--fastq-extension` → Extension of FASTQ files (default: `fq.gz`).
- `--metaphlan-db` → Path to the Metaphlan database.
- `--metaphlan-index` → Name of the Metaphlan database index.
- `--n-db` → Path to the nucleotide search database.
- `--protein-db` → Path to the protein search database.
- `--max_workers` → Number of parallel workers (default: 4).
- `--joblist-size` → Number of commands per job list file (default: 2).

## Workflow
1. Extracts common prefixes from FASTQ files in the reads directory.
2. Generates MeTAline configuration files for each prefix.
3. Executes the configuration generation commands in parallel.
4. Constructs and writes job list files for greasy execution.
5. Writes a SLURM submission script to manage job execution.

## Output Files
- **Configs Directory**: Stores generated MeTAline config files.
- **Job List Files**: Stores job batch files for `greasy`.
- **`list_greasy.txt`**: List of `greasy` execution commands.
- **`launch_metaline.job`**: SLURM job submission script.

## SLURM Job Submission
> [!IMPORTANT]
> It is important for different users to review and adjust the `launch_metaline.job` file according to their HPC cluster's SLURM setup.

Once the script has generated the necessary files, and the `launch_metaline.job` has been reviewed corresponding to your cluster's Slurm setup, submit the job using:
```sh
sbatch launch_metaline.job
```
