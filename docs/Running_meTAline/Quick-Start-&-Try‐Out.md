## Quick-Start
Quickly test meTAline with very little setup, and have a look at the structure and formats of the output results.

### Requirements (recommended)
4 CPU
3 Gb RAM
40 Gb free disk space

---

### With Singularity
1. Install Singularity (This example is for Ubuntu 24.04)
- Download from [HERE](https://github.com/sylabs/singularity/releases/tag/v4.1.5) the corresponding installation package (_.deb or _.rpm) to your operating system. (E.G.: Ubuntu 24.04 needs the "singularity-ce_4.1.5-noble_amd64.deb" file).
- Install it from the downloaded package. Example command for the above mentioned version:
```bash
sudo apt install ./singularity-ce_4.1.5-noble_amd64.deb
```
- Test if the installation is OK.
```bash
singularity --version
```

2. Clone the meTAline repository
```sh
git clone https://github.com/gabaldonlab/metaline
```

3. Change your working directory for meTAline
```sh
cd metaline
```

4. Get Singularity image
If you are on Debian based Linux (like Ubuntu) with x86-64 CPU architecture, then you can simply download the provided [Singularity image](https://github.com/Gabaldonlab/meTAline/releases/download/v1.2.0/metaline.sif)

If you are with other kind of CPU architecture or operating system, other Linux flavors or Mac (Windows is not supported by Singularity), then we recommend a fresh build, with the following command:

```sh
make singularity
```
Which will create a `./metaline.sif` file in the root of your cloned repository.

5. Download the provided test datasets
```sh
make download-test-datasets
```
*OR*
```sh
./_scripts/download_test_datasets
```
6. Generate the required files and folders
mkdir -p ./Benchmark
mkdir -p ./Benchmark/plots
[ -f ./Benchmark/required.txt ] || touch ./Benchmark/required.txt

7. Execute meTAline with the following commands:
- Generate the configuration file to run the job:
```sh
singularity run --cleanenv ./metaline.sif metaline-generate-config \
    --config-file test_run.json \
    --fastq-extension .fq.gz \
    --basedir ./test_output \
    --reads-directory  ./test_input \
    --reference-genome ./test_input/test_datasets/grch38_index/genome \
    --krakendb ./test_input/test_datasets/minikraken2_v2_8GB_201904_UPDATE \
    --sample-barcode test_run \
    --fastq-prefix V300091236_L01_100 \
    --metaphlan-db ./test_input/test_datasets/metaphlan_dbs \
    --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db ./test_input/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz \
    --protein-db ./test_input/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz
```

- Run the job with the previously generated configuration file:
```sh
singularity run --cleanenv ./metaline.sif metaline \
    kmer_taxonomy \
    -j 16 \
    --configfile ./test_run.json
```

---

### With Docker
1. Install Docker
Follow the [Official documentation](https://docs.docker.com/engine/install/)

2. Clone the meTAline repository
```sh
git clone https://github.com/gabaldonlab/metaline
```

3. Change your working directory for meTAline
```sh
cd metaline
```

4. Get Docker image
There are two ways to get the meTAline Docker image.

- Pull it from the Dockerhub registry:
```sh
docker pull cgenomics/metaline:latest
```

- Build it from scratch:
```sh
make docker-image  # Which translates to: docker image build -t metaline:latest .
```

5. Download the provided test datasets
```sh
make download-test-datasets
```
*OR*
```sh
./_scripts/download_test_datasets
```

6. Execute meTAline with the following commands:
- Generate the configuration file to run the job:
```sh
docker run --rm \
    -v ./test_input:/input-directory \
    -v ./test_output:/output-directory \
    cgenomics/metaline:latest \
    metaline-generate-config \
    --config-file /output-directory/test_run.json \
    --fastq-extension .fq.gz \
    --basedir /output-directory/ \
    --reads-directory  /input-directory/ \
    --reference-genome /input-directory/test_datasets/grch38_index/genome \
    --krakendb /input-directory/test_datasets/minikraken2_v2_8GB_201904_UPDATE \
    --sample-barcode test_run \
    --fastq-prefix V300091236_L01_100 \
    --metaphlan-db /input-directory/test_datasets/metaphlan_dbs \
    --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202307 \
    --n-db /input-directory/test_datasets/chocophlan_EC_FILTERED.v4_alpha.tar.gz \
    --protein-db /input-directory/test_datasets/uniref90_annotated_v4_alpha_ec_filtered.tar.gz
```

- Run the job with the previously generated configuration file:
```sh
docker run --rm \
    -v ./test_input:/input-directory \
    -v ./test_output:/output-directory \
    cgenomics/metaline:latest \
    metaline \
        kmer_taxonomy \
        -j 16 \
        --configfile /output-directory/test_run.json # Indicate where the created config file is.
```

---

## Related directories
To try the pipeline we provide test mock samples and a selection of external datasets to download:

`test_input` directory:
```sh
tree ./test_input/
./test_input/
├── explanation.txt
├── test_datasets
│   ├── chocophlan_EC_FILTERED.v4_alpha.tar.gz
│   ├── grch38_genome.tar.gz
│   ├── metaphlan_dbs
│   │   ├── mpa_vJun23_CHOCOPhlAnSGB_202307_marker_info.txt.bz2
│   │   ├── mpa_vJun23_CHOCOPhlAnSGB_202307.md5
│   │   ├── mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt.bz2
│   │   └── mpa_vJun23_CHOCOPhlAnSGB_202307.tar
│   ├── minikraken2_v2_8GB_201904.tgz
│   └── uniref90_annotated_v4_alpha_ec_filtered.tar.gz
├── V300091236_L01_100_1.fq.gz
└── V300091236_L01_100_2.fq.gz
```

An output example of the pipeline, after running all the rules:
```sh
metaline-test-output/
├── BAM
│   ├── mock.hisat2.bam
│   ├── mock.hisat2.bam.bai
│   └── mock.unmapped.fastq.gz
├── Benchmark
│   ├── 20241108.102621.871345_mock.alignment.benchmark.txt
│   ├── 20241108.102621.871345_mock.biom.benchmark.txt
│   ├── 20241108.102621.871345_mock.concat.benchmark.txt
│   ├── 20241108.102621.871345_mock.extracted_reads.benchmark.txt
│   ├── 20241108.102621.871345_mock.fastqc.benchmark.txt
│   ├── 20241108.102621.871345_mock.krona.benchmark.txt
│   ├── 20241108.102621.871345mock.R_alpha_div.benchmark.txt
│   ├── 20241108.102621.871345mock.R_barplot.benchmark.txt
│   ├── 20241108.102621.871345mock.R_conversion.benchmark.txt
│   ├── 20241108.102621.871345_mock.V300091236_L01_100.trimming.benchmark.txt
│   └── mock.kraken2.benchmark.txt
├── err_out
│   ├── jobnm_11453585.err
│   └── jobnm_11453585.out
├── EXTRACTED_FASTA
│   ├── mock1.fastq
│   ├── mock.1.fastq.gz
│   ├── mock2.fastq
│   └── mock.2.fastq.gz
├── KRAKEN_ASSIGN
│   ├── mock.Kraken2.biom
│   ├── mock.kraken2.report
│   └── mock.kraken2.txt
├── KRONA_HTML
│   ├── mock.html
│   └── mock.krona
├── METAPHLAN4
│   ├── mock
│   │   ├── mock.unmapped_genefamilies.tsv
│   │   ├── mock.unmapped_humann_temp
│   │   │   └── mock.unmapped.log
│   │   ├── mock.unmapped_pathabundance.tsv
│   │   └── mock.unmapped_pathcoverage.tsv
│   ├── mock.bz2
│   ├── mock_profiled.txt
│   ├── mocksam.bz2
│   └── mock.vsc.txt
├── R_ANALYSIS
│   ├── mock_alpha_div.csv
│   ├── mock_Barplot_phyla.jpeg
│   └── mock.rds
├── TRIMMOMATIC
│   ├── mock.1.fastq.gz
│   ├── mock.2.fastq.gz
│   ├── mock_qc
│   │   ├── mock.1_fastqc.html
│   │   ├── mock.1_fastqc.zip
│   │   ├── mock.2_fastqc.html
│   │   └── mock.2_fastqc.zip
│   ├── V300091236_L01_100_1_paired.fq.gz
│   ├── V300091236_L01_100_1_unpaired.fq.gz
│   ├── V300091236_L01_100_2_paired.fq.gz
│   └── V300091236_L01_100_2_unpaired.fq.gz
```

---
