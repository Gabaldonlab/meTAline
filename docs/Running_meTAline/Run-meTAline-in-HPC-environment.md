## Running in HPC environment <a id="running-in-hpc" />

The MeTAline pipeline is intended to be used in HPC environment, because of its high RAM runtime requirement.


This is a continuation of the step by step instruction in docs/Run-meTAline-with-singularity.md. 


### 3. Creation of the config file and running of the pipeline
There are two options at this point: Parallelisation or running the jobs individually

**Parallelisation**
meTAline is written for single input execution, without managing it's own parallelization, which simplifies development, debugging, and maintenance while avoiding parallelization issues.
If you want to parallelize multiple meTAline executions, you can make use of native features of snakemake or external tools like `GREASY arrays with SLURM` provide flexible parallelization, allowing users to control resource allocation dynamically. This approach enhances scalability, enabling execution across multiple computational nodes or even machines without modifying the tool itself. It also improves fault isolation, robustness, and composability within Unix pipelines and automation scripts. By keeping the tool "single-input", you maximize portability and efficiency while leaving parallel execution to specialized tools.

With meTAline we provide the [metaline-prepare-greasy-array-job](Metaline-prepare-greasy-array-job.md) CLI tool to pre-generate the config files and metaline commands for a greasy array job launch on SLURM based compute nodes.

**Run jobs individually**
Generate a hpc_sbatch.job file. An example can be found at [example_hpc_sbatch.job](../../examples/example_hpc_sbatch.job)
In the following example SLURM .job file we use some Mare Nostrum 5 (Barcelona Supercomputing Center) specific additional configurations that you can ignore.
You can change the required and default parameters your job using the [Parameters documentation](Available-target-rules.md).

The base commands of this script are through the generated Singularity image, which first you need to build on your local machine (because it requires "sudo" permissions, which in most cases not available in HPC machines), then you can copy that built image into the remote HPC machine.

If running jobs individually, it is there are mandatory folders and a file that need to be created:
```sh
mkdir -p ./Benchmark
mkdir -p ./Benchmark/plots
[ -f ./Benchmark/required.txt ] || touch ./Benchmark/required.txt
```
### 4. Launch the job 
If you chose to run using parallelisation, once the necessary files have been generated from the [metaline-prepare-greasy-array-job](Metaline-prepare-greasy-array-job.md) instructions, and the `launch_metaline.job` has been reviewed corresponding to your cluster's Slurm setup, submit the job using:
```sh
sbatch launch_metaline.job
```

Otherwise if you chose to run the jobs individually, run example_hpc_sbatch.job created in step 3:
```bash
sbatch ./example_hpc_sbatch.job
```
To see more SLURM job examples (using array of greasies for parallelization), you can visit the [examples directory]((../../examples/)

