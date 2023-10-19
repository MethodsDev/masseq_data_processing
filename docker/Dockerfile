FROM ubuntu:latest
LABEL authors="akhorgad"

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin
ENV CONDA_DIR=/opt/conda

WORKDIR $SRC


RUN apt-get update --allow-releaseinfo-change
RUN apt-get update -m
RUN apt-get install -y curl git-lfs time datamash wget zlib1g-dev libssl-dev lbzip2 make gcc libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev


# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda

# Adding conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

#adding packages to support pb plots
RUN conda install -c conda-forge -y --quiet matplotlib
RUN conda install -c conda-forge -y --quiet argparse
RUN conda install -c bioconda -y --quiet pysam
RUN conda install -c conda-forge -y --quiet tqdm
RUN conda install -c conda-forge -y --quiet numpy

RUN git clone -b dev_skera-pbplots https://github.com/akankshakhorgade/mdl_data_processing_workflow.git