meTAline is written for single input execution, without managing it's own parallelization, which simplifies development, debugging, and maintenance while avoiding parallelization issues.
If you want to parallelize multiple meTAline executions, you can make use of native features of snakemake or external tools like `GREASY arrays with SLURM` provide flexible parallelization, allowing users to control resource allocation dynamically. This approach enhances scalability, enabling execution across multiple computational nodes or even machines without modifying the tool itself. It also improves fault isolation, robustness, and composability within Unix pipelines and automation scripts. By keeping the tool "single-input", you maximize portability and efficiency while leaving parallel execution to specialized tools.

With meTAline we provide the [metaline-prepare-greasy-array-job](https://github.com/Gabaldonlab/meTAline/wiki/Metaline-prepare-greasy-array-job) CLI tool to pre-generate the config files and metaline commands for a greasy array job launch on SLURM based compute nodes.

> [!IMPORTANT]
> To avoid any unexpected results when finding the sample code prefixes from the fastq.gz files make sure that the `--reads-directory` argument's directory path and its content conforms the following critarias
> - 1. Use ONLY gzipped (with .gz extension) files. Do not `gunzip` them before use!
>
> - 2. FORWARD and REVERSE files of a sample must be in pairs
> ```sh
> # E.G.:
> E100043441_L01_1_1.1.fq.gz
> E100043441_L01_1_1.2.fq.gz
>
> E100043441_L01_1_2.1.fq.gz
> E100043441_L01_1_2.2.fq.gz
>
> # See: The FORWARD is split into 2 files, and the REVERSE also. (N of FORWARD files == N of REVERSE files)
> ```
>
> - 3. Nested directories are allowed, but the FORWARD and REVERSE sequences of a sample must be in the same subdirectory:
> ```sh
> ── TG0002
>    ├── E100043441_L01_1_1.fq.gz
>    └── E100043441_L01_1_2.fq.gz
> ── TG0003
>    ├── E100043441_L01_2_1.fq.gz
>    └── E100043441_L01_2_2.fq.gz
> ── TG0004
>    ├── E100043441_L01_3_1.fq.gz
>    └── E100043441_L01_3_2.fq.gz
> ```
> **OR in case of flat directory**
> ```sh                                                                                                                  │~
> ├── V300091236_L01_100_1.fq.gz
> ├── V300091236_L01_100_2.fq.gz
> ```

**EXAMPLE**

1. Download samples (or use yours):
```sh
cd ./test_input

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/003/ERR1018263/ERR1018263_2.fastq.gz  # ERR1018263 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/003/ERR1018263/ERR1018263_1.fastq.gz  # ERR1018263 read 2

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/004/ERR1018234/ERR1018234_2.fastq.gz  # ERR1018234 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/004/ERR1018234/ERR1018234_1.fastq.gz  # ERR1018234 read 2

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/002/ERR1018192/ERR1018192_1.fastq.gz  # ERR1018192 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/002/ERR1018192/ERR1018192_2.fastq.gz  # ERR1018192 read 2

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1018218/ERR1018218_1.fastq.gz  # ERR1018218 read 1
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1018218/ERR1018218_2.fastq.gz  # ERR1018218 read 2
```

2. Download test datasets (or use yours):
```sh
make download-test-datasets
```

3. In the HPC cluster generate the greasy array job setup files:

> [!IMPORTANT]
> The arguments for `--metaline-cmd` and `--greasy-cmd` must be quoted (E.G.: `--generate_config_cmd "singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.sif metaline-generate-config"`)

> [!NOTE]
> It is advised to keep some of the "--cluster-node-cores" free for greasy and slurm to manage the executions.

```sh
singularity run --cleanenv ./metaline.sif metaline-prepare-greasy-array-job \
    --metaline-rule "all" \
    --basedir ~/metaline-prepare-greasy-array-job/test_output \
    --metaline-cmd "module load singularity && singularity run --cleanenv ~/project/pipelines/meTAline/meTAline-1.1.0/metaline.sif metaline" \
    --greasy-cmd "/apps/GPP/GREASY/2.2.4.1/INTEL/IMPI/bin/greasy" \
    --reference-genome ./test_input/test_datasets/grch38_index/genome \
    --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE \
    --reads-directory  ./test_input \
    --fastq-extension .fq.gz \
    --metaphlan-db ./test_input/test_datasets/metaphlan_dbs \
    --metaphlan-index  mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz \
    --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz \
    --joblist-size 2 `# Max. number of jobs per each node in the Greasy array.` \
    --cluster-node-cores 24 `# Adjust these cores according to your cluster node's resouces.` \
    --trimmo-cores 16 `# Adjust these cores according to the --cluster-node-cores.` \
    --hisat2-cores 16 `# Adjust these cores according to the --cluster-node-cores.` \
    --kraken2-cores 16 `# Adjust these cores according to the --cluster-node-cores.` \
    --snakemake-cores 16 `# Adjust these cores according to the --cluster-node-cores.`
```

> [!IMPORTANT]
> Please note that while we designed the pipeline to be as generalizable as possible, we identified sample naming as a potential source of issues. To address this, we introduced the `fastq-extension` parameter.We recommend setting the fastq-extension to either ".fastq.gz" or ".fq.gz". By default, meTAline uses ".fastq.gz".Additionally, we added the `--prefix-fr` parameter, which allows you to specify the naming convention used to distinguish forward and reverse reads.
> By default the value will be "_", but the parameter can take multiple arguments as `--prefix-fr "_" ".R" "."` .
>
> For example:
> If _ is used before indicating 1 (for forward) and 2 (for reverse): {file}_1{fastq-extension}
>
> If .R is used as the prefix: {file}.R1{fastq-extension}
>
> Examples:
> For sample G12, you have these two files as raw data: **sample_G12_1.fastq.gz**  and **sample_G12_2.fastq.gz**
>
> In this case you will indicate `fastq-extension .fastq.gz` and `--prefix-fr "_"`, but in case the raw files are: "sample_G12.R1.fq.gz" and "sample_G12.R1.fq.gz", then:
> `fastq-extension .fq.gz` and `--prefix-fr ".R"`

> [!IMPORTANT]
> It is important for different users to review and adjust the `launch_metaline.job` file according to their HPC cluster's SLURM setup.

4. Once the script has generated the necessary files, and the `launch_metaline.job` has been reviewed corresponding to your cluster's Slurm setup, submit the job using:
```sh
sbatch launch_metaline.job
```

**OUTPUT GREASY JOB SETUP**
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
