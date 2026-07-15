FROM continuumio/miniconda3

LABEL authors="Joan Sarria, Bruno Contreras Moreira"
LABEL description="Container with align2graph, PHG and selected pangenome graph [Med13, Pan20, Example_Ara]"
LABEL license="GPLv3"
LABEL maintainers="compbio@eead.csic.es"
LABEL org.opencontainers.image.source https://github.com/OWNER/REPO

USER root

# bound in runtime to writable local dir with 'index'
RUN mkdir -m 0777 /gmap_db

RUN mkdir /barleygraph
WORKDIR /barleygraph

## install system deps as root
RUN apt update && apt install --no-install-recommends -y \
    build-essential \
    openjdk-17-jdk-headless \
    nano \
    && apt clean && apt purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# PHG v2.4.72.229 - Download and extract
#RUN wget -qO - https://api.github.com/repos/maize-genetics/phg_v2/releases/tags/2.4.74.229\
#    | awk -F': ' '/browser_download_url/ && /\.tar/ {gsub(/"/, "", $(NF)); system("wget " $(NF))}' \
#    && tar -xvf PHGv2-v2.4.tar

# PHG v2.4 (last 2.4 release) - Download and extract
RUN wget -qO - https://api.github.com/repos/maize-genetics/phg_v2/releases/tags/2.4.75.230 \
    | awk -F': ' '/browser_download_url/ && /\.tar/ {gsub(/"/, "", $(NF)); system("wget " $(NF))}' \
    && tar -xvf *.tar \
    && rm *.tar

ENV PHGPATH /barleygraph/phg/bin

# PHG and python deps
RUN wget -qO environment.yml https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/environment.yml
RUN conda env create -f environment.yml
RUN conda install -n phgv2.4 -c conda-forge -y pandas matplotlib tqdm pyyaml
RUN echo "conda activate phgv2.4" >> ~/.bashrc
ENV PATH /opt/conda/envs/phgv2.4/bin:$PATH
ENV CONDA_DEFAULT_ENV phgv2.4

# GMAP 2013-08-31 makes tiny indexes compared to 2015-12-31 onwards
ENV GVER 2013-08-31
ENV GPATH /barleygraph/gmap-${GVER}/local
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-${GVER}.tar.gz
RUN tar xfz gmap-gsnap-${GVER}.tar.gz && cd gmap-${GVER} && ./configure --prefix=${GPATH} && make && make install && make clean
RUN rm gmap-gsnap-${GVER}.tar.gz

# align2graph 
#RUN wget https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/refs/heads/main/align2graph.py
#Get from local meanwhile, will swap to GitHub once the script is public
COPY docker/scripts/align2graph.py align2graph.py
COPY version.py ./
# Dummy until all updated in github
RUN chmod +x align2graph.py
RUN mv align2graph.py align2graph

ARG graph
LABEL graph=$graph
RUN echo $graph

# graphs, config and bits - download from GitHub repository
#RUN wget -qO graph.yaml https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/graphs/${graph}.yaml
# Getting them locally, will swap to GitHub once the graphs are public
COPY graphs/Pan20_gmap-geno.yaml ./
COPY graphs/Pan20_mmap-pro.yaml ./
ADD ${graph}.tgz .
# Placeholder, until the graph.tgz files are public and available at github
COPY test.fa test.fa

# Sample list files - download from GitHub repository based on graph type
# RUN wget -qO samplelist.tsv https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/graphs/${graph}_samplelist.tsv
# Get it locally untill all public
COPY graphs/${graph}_samplelist.tsv samplelist.tsv

# scripts - download from GitHub repository
RUN mkdir -p python_scripts
RUN wget -qO gmap_build.sh https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/gmap_build.sh
# RUN wget -qO build_imputation_index.sh https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/build_imputation_index.sh
# Local placegolder
COPY docker/scripts/build_imputation_index.sh build_imputation_index.sh
#RUN wget -qO imputation.sh https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/imputation.sh
# Local placeholder
COPY docker/scripts/imputation.sh imputation.sh
RUN wget -qO python_scripts/haplopainting.py https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/haplopainting.py
RUN wget -qO python_scripts/hvcf2bed.py https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/hvcf2bed.py
COPY docker/scripts/index_${graph} index

# Ensure all scripts are executable
RUN chmod +x gmap_build.sh build_imputation_index.sh imputation.sh index python_scripts/haplopainting.py python_scripts/hvcf2bed.py

# Remove .sh extensions from shell scripts
RUN mv build_imputation_index.sh build_imputation_index && \
    mv imputation.sh imputation

# Create wrapper scripts that use conda run
RUN echo '#!/bin/bash' > /barleygraph/haplopainting && \
    echo 'exec conda run -n phgv2.4 python3 /barleygraph/python_scripts/haplopainting.py "$@"' >> /barleygraph/haplopainting && \
    chmod +x /barleygraph/haplopainting

RUN echo '#!/bin/bash' > /barleygraph/hvcf2bed && \
    echo 'exec conda run -n phgv2.4 python3 /barleygraph/python_scripts/hvcf2bed.py "$@"' >> /barleygraph/hvcf2bed && \
    chmod +x /barleygraph/hvcf2bed

ENV PATH="/barleygraph:${GPATH}/bin:${PHGPATH}:${PATH}"
SHELL ["conda", "run", "-n", "phgv2.4", "/bin/bash", "-c"]
