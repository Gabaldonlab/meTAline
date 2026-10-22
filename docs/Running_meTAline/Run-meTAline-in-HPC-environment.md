## Running in HPC environment <a id="running-in-hpc" />

The MeTAline pipeline is intended to be used in HPC environment, because of its high RAM runtime requirement.
Before starting this step by step process you should have completed the build and setup process. Further information can be found at [Build meTAline singularity image](docs/Building_images/Build-meTAline-singularity-image.md)


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

**Step by step:**

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

3. There are two options at this point. 

**Parallelisation**
meTAline is written for single input execution, without managing it's own parallelization, which simplifies development, debugging, and maintenance while avoiding parallelization issues.
If you want to parallelize multiple meTAline executions, you can make use of native features of snakemake or external tools like `GREASY arrays with SLURM` provide flexible parallelization, allowing users to control resource allocation dynamically. This approach enhances scalability, enabling execution across multiple computational nodes or even machines without modifying the tool itself. It also improves fault isolation, robustness, and composability within Unix pipelines and automation scripts. By keeping the tool "single-input", you maximize portability and efficiency while leaving parallel execution to specialized tools.

With meTAline we provide the [metaline-prepare-greasy-array-job](docs/Running_meTAline/Metaline-prepare-greasy-array-job.md) CLI tool to pre-generate the config files and metaline commands for a greasy array job launch on SLURM based compute nodes.

**Run jobs individually**
Generate a hpc_sbatch.job file. An example can be found at [example_hpc_sbatch.job](examples/example_hpc_sbatch.job)
In the following example SLURM .job file we use some Mare Nostrum 5 (Barcelona Supercomputing Center) specific additional configurations that you can ignore.
Required and default parameters for this job can be found at (docs/Running_meTAline/Available-target-rules.md).

The base commands of this script are through the generated Singularity image, which first you need to build on your local machine (because it requires "sudo" permissions, which in most cases not available in HPC machines), then you can copy that built image into the remote HPC machine.

If running jobs individually, it is there are mandatory folders and a file that need to be created:
```sh
mkdir -p ./Benchmark
mkdir -p ./Benchmark/plots
[ -f ./Benchmark/required.txt ] || touch ./Benchmark/required.txt
```

4. If you chose to run using parallelisation, once the necessary files have been generated from the prepare_greasy_array_job instructions, and the `launch_metaline.job` has been reviewed corresponding to your cluster's Slurm setup, submit the job using:
```sh
sbatch launch_metaline.job
```

Otherwise if you chose to run the jobs individually, run example_hpc_sbatch.job created in step 3
```bash
sbatch ./example_hpc_sbatch.job
```
To see more SLURM job examples (using array of greasies for parallelization), you can visit the [examples directory](https://github.com/Gabaldonlab/meTAline/tree/main/examples) !

