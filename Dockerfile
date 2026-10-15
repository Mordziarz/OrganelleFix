FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    bash \
    wget \
    bzip2 \
    ca-certificates \
    make \
    gcc \
    g++ \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh

ENV PATH=/opt/conda/bin:$PATH

RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    numpy \
    scipy \
    matplotlib \
    biopython \
    mafft \
    blast \
    iqtree \
    diamond \
    mmseqs2 \
    last \
    r-data.table \
    r-base \
    r-ggplot2 \
    r-dplyr \
    r-seqinr \
    r-rstatix \
    r-ggpubr \
    r-patchwork \
    r-forcats \
    r-stringr \
    r-devtools \
    r-remotes \
    r-circlize \
    r-car \
    bioconductor-complexheatmap \
    bioconductor-ggtree \
    bioconductor-biostrings \
    bioconductor-genomicfeatures \
    bioconductor-genomicranges \
    bioconductor-rsamtools \
    bioconductor-iranges \
    bioconductor-rtracklayer \
    && mamba clean -a -y

RUN R -e "if (!requireNamespace('smplot2', quietly = TRUE)) remotes::install_github('smin95/smplot2', upgrade='never')"
RUN R -e "if (!requireNamespace('phylogram', quietly = TRUE)) install.packages('phylogram', repos='https://cloud.r-project.org')"
RUN R -e "if (!requireNamespace('metablastr', quietly = TRUE)) remotes::install_github('drostlab/metablastr', upgrade='never')"
RUN R -e "if (!requireNamespace('RSCUcaller', quietly = TRUE)) remotes::install_github('Mordziarz/RSCUcaller', upgrade='never')"
RUN R -e "if (!requireNamespace('transgenomes', quietly = TRUE)) remotes::install_github('Mordziarz/transgenomes', upgrade='never')"
