meTAline is written for single input execution, without managing it's own parallelization, which simplifies development, debugging, and maintenance while avoiding parallelization issues. If you want to parallelize multiple meTAline executions, you can use native features of snakemake or use external tools like `GNU Parallel`, `xargs -P` that provide flexible parallelization, allowing users to control resource allocation dynamically. This approach enhances scalability, enabling execution across multiple cores or even machines without modifying the tool itself. It also improves fault isolation, robustness, and composability within Unix pipelines and automation scripts. By keeping the tool "single-input", you maximize portability and efficiency while leaving parallel execution to specialized tools.

## Example with Xargs
1. Create a `commands.txt` file with your commands prepared.

**EXAMPLE**
> [!NOTE]
> To avoid any problems from miswriting newlines, it is suggested to leave the indivial job commands in one line (`metaline-generate-config` and `metaline` commands are joined together with "&&"!)

```sh
singularity run --cleanenv ./metaline.sif metaline-generate-config --config-file test_run.json --extension fq.gz --basedir ./test_output --reads-directory  ./test_input --reference-genome ./test_input/test_datasets/grch38_index/genome --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE --sample-barcode test_run --fastq-prefix V300091236_L01_100 --metaphlan-db ./test_input/test_datasets/metaphlan_dbs --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz && singularity run --cleanenv ./metaline.sif metaline kmer_taxonomy -j 16 --configfile ./test_run.json

singularity run --cleanenv ./metaline.sif metaline-generate-config --config-file test_run.json --extension fq.gz --basedir ./test_output --reads-directory  ./test_input --reference-genome ./test_input/test_datasets/grch38_index/genome --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE --sample-barcode test_run --fastq-prefix V300091236_L01_100 --metaphlan-db ./test_input/test_datasets/metaphlan_dbs --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz && singularity run --cleanenv ./metaline.sif metaline kmer_taxonomy -j 16 --configfile ./test_run.json
```

2. Run the `commands.txt` file with `parallel`

```sh
xargs -I {} -P 2 sh -c "{}" < commands.txt  # Note the -P flag that controlls how many commands should be running at once during execution!
```

## Example with GNU Parallel
1. Install `GNU Parallel`
- `GNU Parallel` is a powerful tool for running commands in parallel. It is not always pre-installed on Unix-based systems, but you can install it easily using package managers.

## **Linux**
`Debian/Ubuntu`
```sh
sudo apt update && sudo apt install parallel
```

`RHEL/CentOS/Fedora`
```sh
sudo dnf install parallel  # For Fedora
```

```sh
sudo yum install parallel  # For CentOS/RHEL
```

`Arch Linux`
```sh
sudo pacman -Sy parallel
```

`macOS`
```sh
brew install parallel
```

`Windows (via Git Bash, Cygwin, or WSL)`
    `Git Bash or Cygwin`: Install parallel from the official GNU sources.
    `WSL (Windows Subsystem for Linux)`: Install it just like on Ubuntu:

    ```sh
    sudo apt install parallel
    ```

> [!NOTE]
> After installing parallel, the first time you run it, you may see this prompt:
```sh
Will cite

GNU Parallel requests citation acknowledgment. To avoid seeing this message every time, you can run:

parallel --citation

and follow the instructions.
Verify Installation
```

2. Create a `commands.txt` file with your commands prepared.

**EXAMPLE**
> [!NOTE]
> To avoid any problems from miswriting newlines, it is suggested to leave the indivial job commands in one line (`metaline-generate-config` and `metaline` commands are joined together with "&&"!)

```sh
singularity run --cleanenv ./metaline.sif metaline-generate-config --config-file test_run.json --extension fq.gz --basedir ./test_output --reads-directory  ./test_input --reference-genome ./test_input/test_datasets/grch38_index/genome --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE --sample-barcode test_run --fastq-prefix V300091236_L01_100 --metaphlan-db ./test_input/test_datasets/metaphlan_dbs --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz && singularity run --cleanenv ./metaline.sif metaline kmer_taxonomy -j 16 --configfile ./test_run.json

singularity run --cleanenv ./metaline.sif metaline-generate-config --config-file test_run.json --extension fq.gz --basedir ./test_output --reads-directory  ./test_input --reference-genome ./test_input/test_datasets/grch38_index/genome --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE --sample-barcode test_run --fastq-prefix V300091236_L01_100 --metaphlan-db ./test_input/test_datasets/metaphlan_dbs --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz && singularity run --cleanenv ./metaline.sif metaline kmer_taxonomy -j 16 --configfile ./test_run.json
```

3. Run the `commands.txt` file with `parallel`

```sh
parallel -j 2 < commands.txt
```
