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

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 \
    && tar xjf samtools-1.14.tar.bz2 \
    && rm samtools-1.14.tar.bz2 \
    && cd samtools-1.14 \
    && ./configure \
    && make \
    && make install

#install isoeseq3
RUN wget -O /usr/local/bin/isoseq3 https://github.com/PacificBiosciences/IsoSeq/releases/download/v3.8.2/isoseq3 \
    && chmod 755 /usr/local/bin/isoseq3

#commented on 15th Aug,2023 to use skera binary for Revio runs
#install skera
#RUN conda install -c bioconda pbskera

#install samtools
RUN conda install -c bioconda samtools

# install gsutil
RUN curl https://sdk.cloud.google.com | bash
ENTRYPOINT ["top", "-b"]

ENV PATH=/root/google-cloud-sdk/bin:$PATH