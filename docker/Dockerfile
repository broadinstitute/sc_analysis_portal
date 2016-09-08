FROM ubuntu:xenial
MAINTAINER ttickle@broadinstitute.org

######################
## Environment
######################

## Constants
### SOFTWARE versions
ENV PICARD_VERSION 2.0.1
ENV R_VERSION 3.3.1-1xenial0
ENV RNASEQQC_VERSION v1.1.8
ENV RSEM_VERSION v1.2.31
ENV RESOURCES_VERSION v1

### locations
ENV BIN /usr/local/bin
ENV DATA /usr/local/data
ENV SRC /usr/local/src

### URLS
ENV PICARD_URL https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip
ENV RESOURCES_URL https://data.broadinstitute.org/single_cell_portal/rsem/rsem_${RESOURCES_VERSION}.tar.gz

######################
## Dependencies and Tools
######################
##############
## Helper tools
RUN apt-get update && \
    apt-get install -y git && \
    apt-get install -y unzip && \
    apt-get install -y wget

##############
## System tools for devtools
RUN apt-get install -y libssl-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libxml2-dev

##############
## Install Bowtie
RUN apt-get install -y bowtie

##############
## C++ for RSEM
RUN apt-get install -y build-essential gcc-multilib

##############
## FastQC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -P ${BIN} && \
    unzip ${BIN}/fastqc_v0.11.5.zip -d ${BIN} && \
    chmod +x ${BIN}/FastQC/fastqc

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
RUN wget -P ${SRC} ${PICARD_URL} && \
    unzip ${SRC}/picard-tools-${PICARD_VERSION}.zip -d ${SRC} && \
    mv ${SRC}/picard-tools-${PICARD_VERSION}/* ${BIN} && \
    rm -r ${SRC}/picard-tools-${PICARD_VERSION} && \
    rm ${SRC}/picard-tools-${PICARD_VERSION}.zip

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
RUN wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_${RNASEQQC_VERSION}.jar

##############
## Scone
RUN echo "source(\"https://bioconductor.org/biocLite.R\")" > ${SRC}/install_devtools.R && \
    echo "biocLite(\"devtools\")" >> ${SRC}/install_devtools.R && \
    echo "biocLite()" >> ${SRC}/install_devtools.R && \
    echo "biocLite(\"YosefLab/scone\", dependencies=TRUE)" >> ${SRC}/install_devtools.R && \
    Rscript ${SRC}/install_devtools.R
