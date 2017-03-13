# Starting environment
FROM ubuntu:zesty
MAINTAINER ttickle@broadinstitute.org

### Does not come with add-apt_repository
### this installs the command
RUN apt-get update && apt-get install -y software-properties-common python-software-properties debconf-utils

### Other tools generally used
RUN apt-get install -y git apt-transport-https

####### Cromwell Requirements
# Cromwell, requires
## Java 8
### java -version
RUN add-apt-repository -y ppa:webupd8team/java && apt-get update && echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | debconf-set-selections && apt-get install -y oracle-java8-installer

## Scala 2.11.8
### scala -version
RUN apt-get update && apt-get install -y scala

## SBT 0.13.12
### sbt sbtVersion
RUN echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823 && apt-get update && apt-get install -y sbt

# RUN wget http://apt.typesafe.com/repo-deb-build-0002.deb && dpkg -i repo-deb-build-0002.deb && apt-get update

# Need docker running

# Code
WORKDIR /bin
wget https://github.com/broadinstitute/cromwell/releases/download/25/cromwell-25.jar
