## Available target rules <a id="target-rules" />

The target rules currently available to use are:

-   rule **all**: From trimming and filtering of the reads, to host removal (if provided the indexed host genome to make an alignment of the reads) to taxonomy assignment based on both a k-mer and gene marker approach. Functional prediction based on Humann. It also includes other functionalities such as representation of the taxonomy assignment by krona, extraction of desired reads, etc.

- rule **trimming**: Trimming of the reads, quality assessment and concatenation (useful when having samples that are sequenced in different sequencing lanes).
- rule **host_depletion**: This rule is to take trimmed reads and align them to an indexed reference genome, for host substraction.
- rule **kmer_taxonomy**: This rule is used if you already have the concatenated/unmapped reads and you want to perform the taxonomic assignment based on k-mers, without running the initial steps of the pipeline.
- rule **RAnalysis**: This rule is used to make a basic R analysis and plotting alpha diversity.
- rule **BioBakery**: This rule is to perform the taxonomy and functional profiling using Biobakery tools based on gene markers: Metaphlan4 and Humann.
- rule **sylph**: This rule is to perform the taxonomy and functional profiling using Biobakery tools based on gene markers: Metaphlan4 and Humann.


> [!IMPORTANT]
> The arguments for `--metaline-cmd` and `--greasy-cmd` must be quoted (E.G.: `--generate_config_cmd "singularity run --cleanenv /gpfs/projects/bsc40/project/pipelines/meTAline/meTAline-1.1.0/metaline.sif metaline-generate-config"`)

### Required Arguments
- `--basedir` → Base directory for output files.
- `--metaline-cmd` → Command to execute MeTAline.
- `--greasy-cmd` → Command to execute Greasy.
- `--krakendb` → Path to the Kraken database.
- `--reads-directory` → Directory containing FASTQ read files.
- `--fastq-extension` → Extension of FASTQ files (default: `.fastq.gz`).
- `--metaphlan-db` → Path to the Metaphlan database.
- `--metaphlan-index` → Name of the Metaphlan database index.
- `--n-db` → Path to the nucleotide search database.
- `--protein-db` → Path to the protein search database.
- `--max_workers` → Number of parallel workers (default: 4).
- `--joblist-size` → Number of commands per job list file (default: 2).
- `--sylph_db` → Sylph database.
- `--sylph_taxonomy` → Sylph taxonomy database.
- `--prefix-fr` → Values between the sample barcode and 1,2 for paired end reads (default: `_`).
- `--fastq_prefix` → This is not required when running as a greasy array

### Optional Arguments
- `--reference-genome` → Path to the indexed reference genome.
- `--reference-phix` → Path to the indexed phix genome.
- `--ani_score_sylph` → ani score for Sylph (default: 70).
- `--ani_score_bacteria_filtering` → filter sylph output by this ani score- bacteria only (default: 95).
- `--ani_score_fungi_filtering` → filter sylph output by this ani score- fungi (default: 80).
- `--cluster-node-cores` → (default: 40).
- `--kraken2-cores` → cores available for kraken2 (default: 10).
- `--snakemake-cores` → cores available for snakemake (default: 5).



>[!IMPORTANT]
>Please note that while we designed the pipeline to be as generalizable as possible, we identified sample naming as a potential source of issues. To address this, we introduced the `--fastq-extension` parameter.We recommend setting the fastq-extension to either "fastq.gz" or "fq.gz". By default, meTAline uses ".fastq.gz".Additionally, we added the `--prefix-fr` parameter, which allows you to specify the naming convention used to distinguish forward and reverse reads.
>By default the value will be "_" .
>
>**For example:**
>If "_" is used before indicating 1 (for forward) and 2 (for reverse): {file}_1{fastq-extension}
>
>If .R is used as the prefix: {file}.R1{fastq-extension}
>
>**Examples:**
>For sample G12, you have these two files as raw data: "sample_G12_1.fastq.gz"  and "sample_G12_2.fastq.gz"
>
>In this case you will indicate `--fastq-extension` .fastq.gz and `--prefix-fr "_"`, but in case the raw files are: **sample_G12.R1.fq.gz** and **sample_G12.R1.fq.gz**, then:
>`--fastq-extension` fq.gz and `--prefix-fr ".R"`


---
