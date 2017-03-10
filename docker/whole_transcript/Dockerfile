FROM ubuntu:xenial
MAINTAINER ttickle@broadinstitute.org

######################
## Environment
######################

## Constants
### SOFTWARE versions
ENV BOWTIE_VERSION 1.1.2-3
ENV BOWTIE2_VERSION 2.2.6-2
ENV PICARD_VERSION 2.0.1
ENV R_VERSION 3.3.2-1xenial0
ENV RNASEQQC_VERSION v1.1.8
ENV RSEM_VERSION v1.2.31
ENV RESOURCES_VERSION v1
ENV SAMTOOLS_VERSION 1.3.1

### locations
ENV BIN /usr/local/bin
ENV DATA /usr/local/data
ENV SRC /usr/local/src

### URLS
ENV PICARD_URL https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip
ENV RESOURCES_URL https://data.broadinstitute.org/single_cell_portal/rsem/rsem_${RESOURCES_VERSION}.tar.gz

### Update PATHS
ENV PATH=${SRC}/sc_analysis_portal/scripts:${SRC}/sc_analysis_portal/bin:${BIN}/Drop-seq_tools-1.12:${BIN}/RSEM:${BIN}/FastQC:${PATH}
ENV PYTHONPATH=${SRC}/sc_analysis_portal/src/SciEDPiper:${SRC}/sc_analysis_portal/src/SciEDPiper/sciedpiper:${PYTHONPATH}

######################
## Dependencies and Tools
######################
##############
## Helper tools
RUN apt-get update && \
    apt-get install -y apt-utils git unzip wget

##############
## System tools for devtools
## Samtools needs libncurses5-dev libncursesw5-dev
RUN apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev libncurses5-dev libncursesw5-dev

##############
## Install Bowtie
RUN apt-get install -y bowtie=${BOWTIE_VERSION} bowtie2=${BOWTIE2_VERSION}

##############
## C++ for RSEM
RUN apt-get install -y build-essential gcc-multilib

##############
## FastQC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -P ${BIN} && \
    unzip ${BIN}/fastqc_v0.11.5.zip -d ${BIN} && \
    chmod +x ${BIN}/FastQC/fastqc && \
    rm ${BIN}/fastqc_v0.11.5.zip

##############
## JAVA 1.8 for PICARD tools
RUN echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
    echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections && \
    echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee /etc/apt/sources.list.d/webupd8team-java.list && \
    echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee -a /etc/apt/sources.list.d/webupd8team-java.list && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886 && \
    apt-get update && \
    apt-get install -y oracle-java8-installer

##############
## PICARD tools
WORKDIR ${SRC}
RUN wget -P ${SRC} ${PICARD_URL} && \
    unzip ${SRC}/picard-tools-${PICARD_VERSION}.zip -d ${SRC} && \
    mv ${SRC}/picard-tools-${PICARD_VERSION}/* ${BIN} && \
    rm -r ${SRC}/picard-tools-${PICARD_VERSION} && \
    rm ${SRC}/picard-tools-${PICARD_VERSION}.zip

##############
## Samtools
WORKDIR ${BIN}
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 -C ${BIN}
WORKDIR ${BIN}/samtools-${SAMTOOLS_VERSION}
RUN make prefix=/usr/local install

##############
## Perl for RSEM
RUN apt-get install -y perl perl-base

##############
## R for RSEM
## R for scone Need >=3.3
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list && \
    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
    gpg -a --export E084DAB9 | apt-key add - && \
    apt-get update && \
    apt-get install -y r-base=${R_VERSION}

##############
## RSEM
WORKDIR ${BIN}
RUN git clone https://github.com/deweylab/RSEM.git && \
    cd ${BIN}/RSEM && \
    git checkout ${RSEM_VERSION} && \
    make

##############
## RNASEQ QC
WORKDIR ${BIN}
RUN wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_${RNASEQQC_VERSION}.jar && \
    mv ${BIN}/RNA-SeQC_${RNASEQQC_VERSION}.jar ${BIN}/RNA-SeQC.jar

##############
### Pull in the code for the analysis
WORKDIR ${SRC}
RUN git clone --recursive https://github.com/broadinstitute/sc_analysis_portal.git
