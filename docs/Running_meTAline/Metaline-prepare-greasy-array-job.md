## Overview
This script automates the preparation of a greasy array job for submitting multiple MeTAline jobs in a high-performance computing (HPC) cluster.

## Workflow
1. Extracts common prefixes from FASTQ files in the reads directory.
2. Generates MeTAline configuration files for each prefix.
3. Executes the configuration generation commands in parallel.
4. Constructs and writes job list files for greasy execution.
5. Writes a SLURM submission script to manage job execution.

## Usage
Run the script with the required arguments as follows:

> [!IMPORTANT]
> The arguments for `--metaline-cmd` and `--greasy-cmd` must be quoted (E.G.: `--generate_config_cmd "singularity run --cleanenv /path_to_metaline_directory/metaline.sif metaline-generate-config"`)

> [!NOTE]
> It is advised to keep some of the "--cluster-node-cores" free for greasy and slurm to manage the executions.

```sh
singularity run --cleanenv ./metaline.sif metaline-prepare-greasy-array-job \
    --metaline-rule "all" \
    --basedir ~/metaline-prepare-greasy-array-job/test_output \
    --metaline-cmd "module load singularity && singularity run --cleanenv /path_to_metaline_directory/metaline.sif metaline" \
    --greasy-cmd "/path_to_greasy_directory/greasy" \
    --reference-genome ./test_input/test_datasets/grch38_index/genome \
    --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE \
    --reads-directory  ./test_input \
    --fastq-extension .fq.gz \
    --metaphlan-db ./test_input/test_datasets/metaphlan_dbs \
    --metaphlan-index  mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz \
    --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz \
    --sylph_db ./test_input/test_datasets/sylph_db.syldb \
    --sylph_taxonomy ./test_input/test_datasets/taxonomy_sylph.tsv \
    --joblist-size 2 `# Max. number of jobs per each node in the Greasy array.` \
    --cluster-node-cores 24 `# Adjust these cores according to your cluster node's resouces.` \
    --trimmo-cores 16 `# Adjust these cores according to the --cluster-node-cores.` \
    --hisat2-cores 16 `# Adjust these cores according to the --cluster-node-cores.` \
    --kraken2-cores 16 `# Adjust these cores according to the --cluster-node-cores.` \
    --snakemake-cores 16 `# Adjust these cores according to the --cluster-node-cores.`
```



Required and default arguments can be found in the [Parameters documentation](Available-target-rules.md).

```
## Output Files
- **Configs Directory**: Stores generated MeTAline config files.
- **Job List Files**: Stores job batch files for `greasy`.
- **`list_greasy.txt`**: List of `greasy` execution commands.
- **`launch_metaline.job`**: SLURM job submission script.

```sh
├── configs
│   ├── config.ERR1018192.json
│   ├── config.ERR1018218.json
│   ├── config.ERR1018234.json
│   └── config.ERR1018263.json
├── joblist0.txt
├── joblist1.txt
├── launch_metaline.job
└── list_greasy.txt
```

> [!IMPORTANT]
> It is important for different users to review and adjust the `launch_metaline.job` file according to their HPC cluster's SLURM setup.
