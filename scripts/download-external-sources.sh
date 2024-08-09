#!/bin/sh

mkdir -p external-sources

# samtools/htslib
wget --directory-prefix=./external-sources "https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2" &

# samtools/samtools
wget --directory-prefix=./external-sources "https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2" &

# usadellab/Trimmomatic
wget --directory-prefix=./external-sources "https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip" &

# fastqc
wget --directory-prefix=./external-sources "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip" &

# DaehwanKimLab/hisat2
wget --directory-prefix=./external-sources "https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz" &

# DerrickWood/kraken2
wget --directory-prefix=./external-sources "https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz" &

# jenniferlu717/KrakenTools
wget --directory-prefix=./external-sources "https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.tar.gz" &

# jenniferlu717/Bracken
wget --directory-prefix=./external-sources "https://github.com/jenniferlu717/Bracken/archive/refs/tags/v2.9.tar.gz" &

# marbl/Krona
wget --directory-prefix=./external-sources "https://github.com/marbl/Krona/archive/refs/tags/v2.8.1.tar.gz" &

# Gabaldonlab/extract-kraken-output-reads
git clone https://github.com/Gabaldonlab/extract-kraken-output-reads ./external-sources/extract_kraken_output_reads && rm -rf ./external-sources/extract_kraken_output_reads/.git

# heloint/MetaPhlAn (Fork of the original: biobakery/MetaPhlAn)
git clone https://github.com/heloint/MetaPhlAn ./external-sources/MetaPhlAn
