FROM ubuntu:latest
LABEL authors="akhorgad"

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin
ENV CONDA_DIR=/opt/conda

WORKDIR $SRC

RUN apt-get update --allow-releaseinfo-change
RUN apt-get update -m
RUN apt-get install -y curl git-lfs time datamash wget bzip2 zlib1g-dev libssl-dev lbzip2 make gcc libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev

RUN ln -s /usr/lib/x86_64-linux-gnu/libtiff.so.6 /usr/lib/x86_64-linux-gnu/libtiff.so.5 

## Install Samtools
ENV SAMTOOLS_VERSION 1.13
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && make && make install

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda

# Adding conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# install gsutil
RUN curl https://sdk.cloud.google.com | bash
ENV PATH=/root/google-cloud-sdk/bin:$PATH

#install skera
RUN conda install -y pbskera -c bioconda

# install lima
RUN conda install -y lima -c bioconda

# install ISOSEQ
RUN conda install -y isoseq -c bioconda

# install pysam
RUN pip install pysam

#adding packages to support pb plots
RUN conda install -y matplotlib -c conda-forge
RUN conda install -y tqdm -c conda-forge
RUN conda install -y numpy -c conda-forge

# install other dependencies
RUN conda install -y pandas -c conda-forge
RUN conda install -y seaborn -c conda-forge

RUN git clone -b main https://github.com/MethodsDev/masseq_data_processing.git
