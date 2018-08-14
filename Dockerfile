##### 
#
# sppIDer
# - Quinn Langdon's sppIDer bioinformatics pipeline
# - https://hub.docker.com/r/glbrc/sppider/
# - https://github.com/GLBRC/sppIDer
#
#####


### base image
# https://hub.docker.com/_/centos/
FROM centos:7.5.1804


### install prereqs
# epel-release is required first to install the yum repo for R
RUN yum update -y && \
    yum install -y \
    epel-release && \
    yum install -y \
    bzip2 \
    bzip2-devel \
    cairo-devel \
    ed \
    gcc \
    gcc-c++ \
    libcurl \
    libcurl-devel \
    make \
    mariadb-devel \
    ncurses-devel \
    openssl-devel \
    postgresql-devel \
    R \
    wget \
    xz-devel \
    zlib-devel && \
    yum clean all && \
    rm -rf /var/cache/yum


### samtools 1.6
# build and install to /usr/local/bin/
RUN cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 && \
    tar -xjvf samtools-1.6.tar.bz2 && \
    rm -f samtools-1.6.tar.bz2 && \
    cd samtools-1.6/ && \
    make && make install  && \
    rm -rf /tmp/samtools-1.6


### bwa 0.7.12
# build and install to /usr/local/bin/
RUN cd /tmp && \
    wget https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2 && \
    tar -xjvf  bwa-0.7.12.tar.bz2 && \
    rm -f bwa-0.7.12.tar.bz2 && \
    cd bwa-0.7.12 && \
    make && \
    cp bwa /usr/local/bin/ && \
    rm -rf /tmp/bwa-0.7.12


### bedtools 2.26
# build and install to /usr/local/bin/
RUN cd /tmp && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && \
    tar -zxvf bedtools-2.26.0.tar.gz && \
    rm -f bedtools-2.26.0.tar.gz && \
    cd bedtools2 && \
    make && make install && \
    rm -rf /tmp/bedtools2


### R 
# install needed packages
RUN R --vanilla -e 'install.packages("data.table", dependencies=TRUE, repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("doBy", dependencies=TRUE, repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("dplyr", dependencies=TRUE, repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("ggplot2", dependencies=TRUE, repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("modes", dependencies=TRUE, repos="http://cran.us.r-project.org")'


### python 2.7 
# install pip
# install biopython and its prereqs
RUN cd /tmp && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    /usr/bin/python2.7 get-pip.py && \
    pip install biopython && \
    pip install numpy


### package sppIDer script files (*.R, *.py)
RUN mkdir -p /tmp/sppIDer/working/
ADD scripts/*.py scripts/*.R /tmp/sppIDer/


### execute sppIDer by default
WORKDIR /tmp/sppIDer/
ENTRYPOINT ["/usr/bin/python2.7"]
CMD ["docker_helper.py"]


