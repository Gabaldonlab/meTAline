FROM python:3.9.1

LABEL base.image="python:3.9.1"
LABEL version="0.7.0"
LABEL software="meTAline"
LABEL software.version="Latest:0.7.0"
LABEL description="Taxonomic assignation pipeline implementation in Snakemake for shotgun genome sequencing "
LABEL website="https://github.com/Dfupa/meTAline"
LABEL license="GNU General Public License 3.0"
LABEL maintainer="Diego Fuentes (BSC)"

#Set up bash and install basic dependencies
SHELL ["/bin/bash", "-c"]

# Install Java 11.
RUN apt-get update
RUN apt-get install -y openjdk-11-jre-headless
RUN apt-get clean

RUN apt-get update -qq
RUN apt-get install -y \
        perl \
        git \
        g++ \
        make \
        nvim \
        wget \
        parallel \
        automake \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        default-jre \
        python3-pip

WORKDIR /bin

# Download sources for:
# - htslib
# - samtools
# - Trimmomatic
# - hisat2
# - Kraken2
# - KrakenTools
# - Bracken
# - Krona
# - fastqc
RUN parallel wget ::: \
        "https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2" \
        "https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2" \
        "https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip" \
        "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip" \
        "https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz" \
        "https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz" \
        "https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.tar.gz" \
        "https://github.com/jenniferlu717/Bracken/archive/refs/tags/v2.9.tar.gz" \
        "https://github.com/marbl/Krona/archive/refs/tags/v2.8.1.tar.gz"

# UN-TAR .bz2 sources for:
# - htslib
# - samtools
RUN parallel tar -vxjf ::: \
        # htslib
        "htslib-1.16.tar.bz2" \
        # samtools
        "samtools-1.13.tar.bz2"


# UN-TAR .gz sources for:
# - hisat2
# - Kraken2
# - KrakenTools
# - Bracken
# - Krona
RUN parallel tar -zxvf ::: \
        # hisat2
        "v2.2.1.tar.gz" \
        # Kraken2
        "v2.1.3.tar.gz" \
        # KrakenTools
        "v1.2.tar.gz" \
        # Bracken
        "v2.9.tar.gz" \
        # Krona
        "v2.8.1.tar.gz"

# UN-ZIP .zip sources for:
# - Trimmomatic
# - fastqc
RUN parallel unzip ::: \
        # Trimmomatic
        "Trimmomatic-0.39.zip" \
        # fastqc
        "fastqc_v0.12.1.zip"

# Remove all compressed files:
# - htslib
# - samtools
# - Trimmomatic
# - hisat2
# - Kraken2
# - KrakenTools
# - Bracken
# - Krona
# - fastqc
RUN parallel rm -rf ::: \
        # htslib
        "htslib-1.16.tar.bz2" \
        # samtools
        "samtools-1.13.tar.bz2" \
        # hisat2
        "v2.2.1.tar.gz" \
        # Kraken2
        "v2.1.3.tar.gz" \
        # KrakenTools
        "v1.2.tar.gz" \
        # Bracken
        "v2.9.tar.gz" \
        # Krona
        "v2.8.1.tar.gz" \
        # Trimmomatic
        "Trimmomatic-0.39.zip" \
        # fastqc
        "fastqc_v0.12.1.zip"

# Install htslib from source
RUN cd htslib-1.16 && make install

# Install samtools from source
RUN cd samtools-1.13 && make install

# Add "trimmomatic" alias the snakemake rules.
RUN echo "alias trimmomatic='java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar'" >> ~/.bash_aliases
RUN echo "source ~/.bash_aliases" >> ~/.bashrc

# Install hisat2 from source
RUN cd hisat2-2.2.1/ \
        && make

# Install Kraken2 from source
RUN cd kraken2-2.1.3/ \
        && ./install_kraken2.sh . \
        && cp ./kraken2{,-build,-inspect} /bin

# Install KrakenTools from source
RUN chmod +x KrakenTools-1.2/* \
        && mv -t /bin KrakenTools-1.2/*.py

# Install Bracken from source
RUN cd Bracken-2.9/ \
        && bash install_bracken.sh \
        && ln -s /bin/Bracken-2.9/src/est_abundance.py /bin/est_abundance.py

# Install Krona from source
RUN cd Krona-2.8.1/KronaTools \
        && mkdir taxonomy \
        && ./install.pl --prefix /bin \
        && mv -t /bin /bin/bin/* \
        && rmdir /bin/bin

# IMPORTANT!
# To avoid a humongous size of the image, the following line is commented.
# If you require to use taxonomy for Kronatools, uncomment the line or execute it in your env.
# RUN cd /bin/Krona-2.8.1/KronaTools && bash updateTaxonomy.sh && bash updateAccessions.sh

# ================================
# END OF REFACTORED PART! TBC HERE
# ================================

##Install from source htslib, samtools, trimmomatic, hisat2, Kraken2, KrakenTools, Bracken, Krona
# RUN cd /bin && wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && tar -vxjf htslib-1.16.tar.bz2 && rm htslib-1.16.tar.bz2 && cd htslib-1.16 && make install
# RUN cd /bin && wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && tar -vxjf samtools-1.13.tar.bz2  && rm samtools-1.13.tar.bz2 && cd samtools-1.13 && make install
# RUN cd /bin && wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && unzip fastqc_v0.12.1.zip  && rm fastqc_v0.12.1.zip
# RUN cd /bin && wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip && \
#         unzip Trimmomatic-0.39.zip && rm Trimmomatic-0.39.zip && \
#         #Added the following alias for nomenclature sake in the snakemake rules
#         echo "alias trimmomatic='java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar'" >> ~/.bash_aliases && \
#         echo "source ~/.bash_aliases" >> ~/.bashrc
# RUN cd /bin && wget https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz && tar -zxvf v2.2.1.tar.gz && rm v2.2.1.tar.gz && cd hisat2-2.2.1/ && make
# RUN cd /bin && wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz && tar -zxvf v2.1.3.tar.gz && rm v2.1.3.tar.gz && cd kraken2-2.1.3/ && ./install_kraken2.sh . && cp ./kraken2{,-build,-inspect} /bin
# RUN cd /bin && wget https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.tar.gz && tar -zxvf v1.2.tar.gz && rm v1.2.tar.gz && chmod +x KrakenTools-1.2/* && mv -t /bin KrakenTools-1.2/*.py
# RUN cd /bin && wget https://github.com/jenniferlu717/Bracken/archive/refs/tags/v2.9.tar.gz && tar -zxvf v2.9.tar.gz && rm v2.9.tar.gz && cd Bracken-2.9/ && bash install_bracken.sh && ln -s /bin/Bracken-2.9/src/est_abundance.py /bin/est_abundance.py
# RUN cd /bin && wget https://github.com/marbl/Krona/archive/refs/tags/v2.8.1.tar.gz  && tar -zxvf v2.8.1.tar.gz && rm v2.8.1.tar.gz && cd Krona-2.8.1/KronaTools && mkdir taxonomy && ./install.pl --prefix /bin && mv -t /bin /bin/bin/* && rmdir /bin/bin
####IMPORTANT: to avoid a humongous size of the image, the following line is commented. If you require to use taxonomy for Kronatools, uncomment the line or execute it in your env
#RUN cd /bin/Krona-2.8.1/KronaTools && bash updateTaxonomy.sh && bash updateAccessions.sh

##Install and upgrade python dependencies with pip as well as installing snakemake from source
RUN python3 -m pip install --upgrade pip && pip3 install --upgrade kraken-biom snakemake==7.32.4 pulp==2.6 numpy pysam biopython pandas pandarallel pybedtools
##Download R from source, compile and set up dependencies. Tends to be a long process. TODO: add the necessary dependencies for the Phylo_R rules
RUN apt-get upgrade -y gfortran libreadline6-dev libx11-dev libxt-dev libpng-dev libjpeg-dev libcairo2-dev xvfb libzstd-dev libcurl4-openssl-dev texinfo texlive texlive-fonts-extra screen libpcre2-dev
RUN test -d /usr/local/src || mkdir -p /usr/local/src && cd /usr/local/src && wget https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz && tar zxvf R-3.6.1.tar.gz && rm R-3.6.1.tar.gz && cd R-3.6.1/ && ./configure --enable-R-shlib && make && make install 
RUN Rscript -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("Biobase", force=TRUE, ask=FALSE)'
RUN Rscript -e 'BiocManager::install(c("phyloseq", "bioformat"), force=TRUE, ask=FALSE)'

#TODO: add the metaline release instead of the version 
RUN test -d /usr/share || mkdir -p /usr/share && cd /usr/share/ && git clone https://github.com/Gabaldonlab/meTAline.git
#Copy the hisat dependencies to /bin path
RUN cd /bin/hisat2-2.2.1/ && cp ./hisat2-build* ./hisat2-inspect* ./hisat2-repeat* ./hisat2-align* ./hisat2 /bin/
#Add pigz
RUN apt-get install -y pigz
#Purge unnecessary dependencies
RUN apt purge -y git make g++ zlib1g-dev python3-pip automake wget make zlib1g-dev libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev libzstd-dev libreadline6-dev libx11-dev libxt-dev libpng-dev libjpeg-dev libcairo2-dev git ant && rm -rf /var/lib/apt/lists/*
#Export patths and test that every dependency works
RUN export PATH="/bin/"
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
RUN extract_kraken_reads.py -h
#Need to fix the path for trimmomatic adding an entrypoint so the bashrc is used in non-login shells. ie: ENTRYPOINT ["/bin/bash", "-c", "source ~/.bashrc && exec \"$@\"", "-"]
#In the meantime, add full path to the rule using sed.
RUN java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -version
RUN /bin/FastQC/fastqc --help
RUN sed -i 's|trimmomatic|java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar|' /usr/share/meTAline/lib/rules/trimming.smk
RUN sed -i 's|fastqc|/bin/FastQC/fastqc|' /usr/share/meTAline/lib/rules/trimming.smk
