# extract_kraken_output_reads

Extracts reads classified by Kraken as a specified taxonomy ID. Those reads are extracted into a new FASTA file.
This tool is direct refactored derivation from the [KrakenTools repository's extract_kraken_reads scripts](https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py)

## 1. Setups on Unix systems

### 1.1 Install from Github

```bash
pip install https://github.com/Gabaldonlab/extract-kraken-output-reads/archive/main.zip
```

### 1.2 Install from source

```bash
git clone https://github.com/Gabaldonlab/extract-kraken-output-reads
cd extract-kraken-output-reads
make install
```

### 1.3 Install in a fresh virtual environment

```bash
git clone https://github.com/Gabaldonlab/extract-kraken-output-reads
cd extract-kraken-output-reads
source scripts/install-with-fresh-env.sh
```

### 1.4 Uninstall

```bash
make uninstall
```

---

## 2. Usage

## 2.1 Flags

```bash
usage: extract_kraken_output_reads [-h] -k KRAKEN_FILE -s SEQ_FILE1 [-s2 SEQ_FILE2] -t TAXID [TAXID ...] -o OUTPUT_FILE [-o2 OUTPUT_FILE2]
                                   [--append] [--noappend] [--max MAX_READS] [-r REPORT_FILE] [--include-parents] [--include-children]
                                   [--exclude] [--fastq-output]

options:
  -h, --help            show this help message and exit
  -k KRAKEN_FILE        Kraken output file to parse
  -s SEQ_FILE1, -s1 SEQ_FILE1, -1 SEQ_FILE1, -U SEQ_FILE1
                        FASTA/FASTQ File containing the raw sequence letters.
  -s2 SEQ_FILE2, -2 SEQ_FILE2
                        2nd FASTA/FASTQ File containing the raw sequence letters (paired).
  -t TAXID [TAXID ...], --taxid TAXID [TAXID ...]
                        Taxonomy ID[s] of reads to extract (space-delimited)
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output FASTA/Q file containing the reads and sample IDs
  -o2 OUTPUT_FILE2, --output2 OUTPUT_FILE2
                        Output FASTA/Q file containig the second pair of reads [required for paired input]
  --append              Append the sequences to the end of the output FASTA file specified.
  --noappend            Create a new FASTA file containing sample sequences and IDs (rewrite if existing) [default].
  --max MAX_READS       Maximum number of reads to save [default: 100,000,000]
  -r REPORT_FILE, --report REPORT_FILE
                        Kraken report file. [required only if --include-parents/children is specified]
  --include-parents     Include reads classified at parent levels of the specified taxids
  --include-children    Include reads classified more specifically than the specified taxids
  --exclude             Instead of finding reads matching specified taxids, finds all reads NOT matching specified taxids
  --fastq-output        Print output FASTQ reads [requires input FASTQ, default: output is FASTA]
```

## 2.2 Examples

**TODO: Add command, input and output examples.**

---

## 3. Citations / Acknowledgments

### 3.1 KrakenTools

[KrakenTools - Github](https://github.com/jenniferlu717/KrakenTools)

[Lu J, Rincon N, Wood D E, Breitwieser F P, Pockrandt C, Langmead B, Salzberg S L, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)]
(https://www.nature.com/articles/s41596-022-00738-y)

Relevant paper for usage of KrakenTools:

1. [Kraken 1](https://github.com/DerrickWood/kraken)
2. [Kraken 2](https://github.com/DerrickWood/kraken2)
3. [KrakenUniq](https://github.com/fbreitwieser/krakenuniq)
4. [Bracken](https://github.com/jenniferlu717/Bracken)

---
