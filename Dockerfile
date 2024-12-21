FROM python:3.10.16-slim-bullseye

LABEL base.image="python:3.10.14-bullseye"
LABEL version="0.8.0-alpha"
LABEL software="meTAline"
LABEL software.version="0.8.0-alpha"
LABEL description="Taxonomic assignation pipeline implementation in Snakemake for shotgun genome sequencing."
LABEL website="https://github.com/gabaldonlab/meTAline"
LABEL license="GNU General Public License 3.0"
LABEL maintainer="Daniel Majer (BSC), Diego Fuentes (BSC)"

RUN apt-get update
RUN apt-get install -y \
        wget \
        curl \
        libcurl4-openssl-dev \
        software-properties-common

# Add apt-fast repository to install apt-fast,
# which will make the package installation faster and concurrent.
RUN /bin/bash -c "$(curl -sL https://git.io/vokNn)"

RUN apt-get update
RUN DEBCONF_NOWARNINGS="yes" \
    TZ="Europe/Madrid" \
    DEBIAN_FRONTEND=noninteractive \
    apt-fast install -y \
        pigz \
        perl \
        git \
        g++ \
        make \
        neovim \
        wget \
        parallel \
        automake \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        default-jre \
        r-base=4.0.4-1 \
        r-base-dev=4.0.4-1 \
        r-cran-ggplot2=3.3.3+dfsg-1 \
        r-bioc-phyloseq=1.34.0+dfsg-1 \
        r-bioc-biobase=2.50.0-1

# - htslib
# - samtools
# - Trimmomatic
# - hisat2
# - Kraken2
# - KrakenTools
# - Bracken
# - Krona
# - fastqc
# - extract_kraken_output_reads
COPY ./external-sources /bin/

# Install htslib from source
WORKDIR /bin
RUN tar -vxjf htslib-1.16.tar.bz2
RUN rm -rf htslib-1.16.tar.bz2
RUN cd /bin/htslib-1.16 && make install

# Install Samtools from source
WORKDIR /bin
RUN tar -vxjf samtools-1.13.tar.bz2
RUN cd /bin/samtools-1.13 && make install
RUN rm -rf samtools-1.13.tar.bz2

# Install Hisat2 from source
WORKDIR /bin
RUN tar -zxvf v2.2.1.tar.gz
RUN cd /bin/hisat2-2.2.1/ && make
RUN cp -r /bin/hisat2-2.2.1/hisat2-build* /bin/
RUN cp -r  /bin/hisat2-2.2.1/hisat2-inspect* /bin/
RUN cp -r  /bin/hisat2-2.2.1/hisat2-repeat* /bin/
RUN cp -r  /bin/hisat2-2.2.1/hisat2-align* /bin/
RUN cp -r  /bin/hisat2-2.2.1/hisat2 /bin/
RUN rm -rf v2.2.1.tar.gz
RUN rm -rf /bin/hisat2-2.2.1

# Install Kraken2 from source
WORKDIR /bin
RUN tar -zxvf v2.1.3.tar.gz
RUN rm -rf v2.1.3.tar.gz
RUN cd /bin/kraken2-2.1.3/ && ./install_kraken2.sh .
RUN cp -t /bin /bin/kraken2-2.1.3/kraken2 /bin/kraken2-2.1.3/kraken2-build /bin/kraken2-2.1.3/kraken2-inspect

# Install KrakenTools from source
WORKDIR /bin
RUN tar -zxvf v1.2.tar.gz
RUN rm -rf v1.2.tar.gz
RUN cd /bin && chmod +x KrakenTools-1.2/*
RUN mv -t /bin KrakenTools-1.2/*.py
RUN rm -rf KrakenTools-1.2

# Install Bracken from source
WORKDIR /bin
RUN tar -zxvf v2.9.tar.gz
RUN rm -rf v2.9.tar.gz
RUN cd /bin/Bracken-2.9/ && bash install_bracken.sh
RUN ln -s /bin/Bracken-2.9/src/est_abundance.py /bin/est_abundance.py

# Install Krona from source
WORKDIR /bin
RUN tar -zxvf v2.8.1.tar.gz
RUN rm -rf v2.8.1.tar.gz
RUN cd /bin/Krona-2.8.1/KronaTools
RUN mkdir taxonomy
RUN cd /bin/Krona-2.8.1/KronaTools && ./install.pl --prefix /bin
RUN mv -t /bin/ /bin/bin/*
RUN rm -rf /bin/bin

WORKDIR /bin
RUN unzip Trimmomatic-0.39.zip
RUN rm -rf Trimmomatic-0.39.zip

WORKDIR /bin
RUN unzip fastqc_v0.12.1.zip
RUN rm -rf fastqc_v0.12.1.zip

# # IMPORTANT!
# # To avoid a humongous size of the image, the following line is commented.
# # If you require to use taxonomy for Kronatools, uncomment the line or execute it in your env.
# # RUN cd /bin/Krona-2.8.1/KronaTools && bash updateTaxonomy.sh && bash updateAccessions.sh

# Install and upgrade python dependencies with pip as well as installing snakemake from source.
RUN python3 -m pip install --upgrade pip
RUN pip3 install --upgrade \
        kraken-biom \
        snakemake==7.32.4 \
        pulp==2.6 \
        numpy \
        pysam \
        biopython \
        pandas \
        pandarallel \
        pybedtools \
        MetaPhlAn==4.1.1 \
        humann==3.9

# Install the extract_kraken_output_reads tool derived from the "extract_reads_of_interest.py" script.
RUN cd /bin/extract_kraken_output_reads && make install

WORKDIR /bin
RUN unzip bowtie2-2.3.5.1-linux-x86_64.zip
RUN rm -rf bowtie2-2.3.5.1-linux-x86_64.zip
RUN cp /bin/bowtie2-2.3.5.1-linux-x86_64/bowtie2* /usr/local/bin/
RUN rm -rf /bin/bowtie2-2.3.5.1-linux-x86_64

COPY ./meTAline /meTAline

# Add bash wrappers to make them global executables.
# This is a better way to do, than creating bash aliases.
RUN echo '#!/bin/bash' >> /bin/trimmomatic
RUN echo 'java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' >> /bin/trimmomatic
RUN chmod 777 /bin/trimmomatic

RUN echo '#!/bin/bash' >> /bin/fastqc
RUN echo '/bin/FastQC/fastqc "$@"' >> /bin/fastqc
RUN chmod 777 /bin/fastqc

RUN echo '#!/bin/bash' >> /bin/metaline-generate-config
RUN echo 'python3 /meTAline/lib/config/generate_config.py "$@"' >> /bin/metaline-generate-config
RUN chmod 777 /bin/metaline-generate-config

RUN echo '#!/bin/bash' >> /bin/metaline
RUN echo 'snakemake -s /meTAline/meTAline.smk "$@"' >> /bin/metaline
RUN chmod 777 /bin/metaline

RUN rm -rf /var/lib/apt/lists/* /var/tmp/*
RUN apt-get clean

# Test that every dependency works
RUN R --version
RUN Rscript --version
RUN snakemake --version
RUN hisat2-build -h
RUN samtools --version
RUN kraken2 -h
RUN kraken-biom -h
RUN est_abundance.py -h
RUN kreport2krona.py -h
RUN ktImportText
RUN extract_kraken_output_reads -h
RUN trimmomatic PE -version
RUN fastqc --help
