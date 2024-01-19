FROM python:3.9.1
#Image in buster flavor
# metadata
LABEL base.image="python:3.9.1"
LABEL version="0.7.0"
LABEL software="meTAline"
LABEL software.version="Latest:0.7.0"
LABEL description="Taxonomic assignation pipeline implementation in Snakemake for shotgun genome sequencing "
LABEL website="https://github.com/Dfupa/meTAline"
LABEL license="GNU General Public License 3.0"
LABEL maintainer="Diego Fuentes (BSC)"
##Set up bash and install basic dependencies
SHELL ["/bin/bash", "-c"]
#RUN apt update && apt install -y software-properties-common
RUN apt-get update && apt-get install -y openjdk-11-jre-headless && apt-get clean
#add-apt-repository ppa:webupd8team/java
RUN apt-get update -qq && apt-get install -y perl default-jre pigz python3-pip make nano automake wget git g++ zlib1g-dev libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev
##Install from source htslib, samtools, trimmomatic, hisat2, Kraken2, KrakenTools, Bracken, Krona
RUN cd /bin && wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && tar -vxjf htslib-1.16.tar.bz2 && rm htslib-1.16.tar.bz2 && cd htslib-1.16 && make install
RUN cd /bin && wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && tar -vxjf samtools-1.13.tar.bz2  && rm samtools-1.13.tar.bz2 && cd samtools-1.13 && make install
RUN cd /bin && wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && unzip fastqc_v0.12.1.zip  && rm fastqc_v0.12.1.zip  
RUN cd /bin && wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip && \
        unzip Trimmomatic-0.39.zip && rm Trimmomatic-0.39.zip && \
        #Added the following alias for nomenclature sake in the snakemake rules
        echo "alias trimmomatic='java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar'" >> ~/.bash_aliases && \
        echo "source ~/.bash_aliases" >> ~/.bashrc
RUN cd /bin && wget https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz && tar -zxvf v2.2.1.tar.gz && rm v2.2.1.tar.gz && cd hisat2-2.2.1/ && make
#RUN cd /bin/ && ls -lhtp *
RUN cd /bin/hisat2-2.2.1/ && ls -lhtp *  && cp ./hisat2-build* ./hisat2-inspect* ./hisat2-repeat* ./hisat2-align* /bin/
RUN cd /bin && wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz && tar -zxvf v2.1.3.tar.gz && rm v2.1.3.tar.gz && cd kraken2-2.1.3/ && ./install_kraken2.sh . && cp ./kraken2{,-build,-inspect} /bin
RUN cd /bin && wget https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.tar.gz && tar -zxvf v1.2.tar.gz && rm v1.2.tar.gz && chmod +x KrakenTools-1.2/* && mv -t /bin KrakenTools-1.2/*.py
RUN cd /bin && wget https://github.com/jenniferlu717/Bracken/archive/refs/tags/v2.9.tar.gz && tar -zxvf v2.9.tar.gz && rm v2.9.tar.gz && cd Bracken-2.9/ && bash install_bracken.sh && ln -s /bin/Bracken-2.9/src/est_abundance.py /bin/est_abundance.py
RUN cd /bin && wget https://github.com/marbl/Krona/archive/refs/tags/v2.8.1.tar.gz  && tar -zxvf v2.8.1.tar.gz && rm v2.8.1.tar.gz && cd Krona-2.8.1/KronaTools && ./install.pl --prefix /bin
##Install and upgrade python dependencies with pip as well as installing snakemake from source
RUN python3 -m pip install --upgrade pip && pip3 install --upgrade kraken-biom snakemake==7.32.4 pulp==2.6 numpy pysam biopython pandas pandarallel pybedtools
##Download R from source, compile and set up dependencies. Tends to be a long process. TODO: add the necessary dependencies for the Phylo_R rules
RUN apt-get upgrade -y gfortran libreadline6-dev libx11-dev libxt-dev libpng-dev libjpeg-dev libcairo2-dev xvfb libzstd-dev libcurl4-openssl-dev texinfo texlive texlive-fonts-extra screen libpcre2-dev
RUN test -d /usr/local/src || mkdir -p /usr/local/src && cd /usr/local/src && wget https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz && tar zxvf R-3.6.1.tar.gz && rm R-3.6.1.tar.gz && cd R-3.6.1/ && ./configure --enable-R-shlib && make && make install 
RUN Rscript -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("phyloseq", "bioformat"), force=TRUE, ask=FALSE)'

#Clone meTAline github
RUN test -d /usr/share || mkdir -p /usr/share && cd /usr/share/ && git clone https://github.com/Gabaldonlab/meTAline.git

#Purge unnecessary dependencies
RUN apt purge -y git make g++ zlib1g-dev python3-pip automake wget curl make zlib1g-dev libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev libzstd-dev libreadline6-dev libx11-dev libxt-dev libpng-dev libjpeg-dev libcairo2-dev git ant && rm -rf /var/lib/apt/lists/*

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
RUN extract_kraken_reads.py -h

#Need to fix the path for trimmomatic adding an entrypoint so the bashrc is used in non-login shells. ie: ENTRYPOINT ["/bin/bash", "-c", "source ~/.bashrc && exec \"$@\"", "-"]
#In the meantime, add full path to the rule using sed.
RUN java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -version
RUN /bin/FastQC/fastqc --help
RUN sed -i 's|trimmomatic|java -jar /bin/Trimmomatic-0.39/trimmomatic-0.39.jar|' /usr/share/meTAline/lib/rules/trimming.smk
RUN sed -i 's|fastqc|/bin/FastQC/fastqc|' /usr/share/meTAline/lib/rules/trimming.smk
