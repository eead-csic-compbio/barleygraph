FROM anaconda/miniconda

LABEL authors="Joan Sarria, Bruno Contreras Moreira"
LABEL description="Container with align2graph and selected tools"
LABEL license="GPLv3"
LABEL maintainers="compbio@eead.csic.es"
LABEL org.opencontainers.image.source="https://github.com/eead-csic-compbio/barleygraph"

USER root

# Bound in runtime to writable local dir. Persistent storage of GMAP db files.
# Required by align2graph only, can be left unset for imputation jobs.
VOLUME /gmap_db

# Bound in runtime to writable local dir. Persistent storage of PHG graph files 
# and possibly results.
VOLUME /graph_db

RUN mkdir /barleygraph
WORKDIR /barleygraph

## install system deps as root
RUN apt update && apt install --no-install-recommends -y \
    build-essential \
    wget \
    nano \
    procps \
    ca-certificates \
    git \
    && apt clean && apt purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Add java via conda
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda install -c conda-forge openjdk=21 -y && \
    conda clean -afy

# PHG v2.4.75.230 - Download and extract (adds 192MB)
RUN URL=$(wget -qO- https://api.github.com/repos/maize-genetics/phg_v2/releases/tags/2.4.75.230 | grep 'browser_download_url' | grep '\.tar' | head -n 1 | cut -d '"' -f 4) && \
    if [ -z "$URL" ]; then echo "ERROR: Could not fetch URL." && exit 1; fi && \
    echo "Downloading from: $URL" && \
    wget -q --show-progress "$URL" && \
    tar -xvf *.tar && \
    rm *.tar

ENV PHGPATH=/barleygraph/phg/bin

# PHG and python deps (adds 2.5GB)
RUN wget -qO environment.yml https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/environment.yml
RUN conda env create -f environment.yml
RUN conda install -n phgv2.4 -c conda-forge -y pandas matplotlib tqdm pyyaml
RUN echo "conda activate phgv2.4" >> ~/.bashrc

# GMAP 2013-08-31 makes tiny indexes compared to 2015-12-31 onwards (adds 121MB)
ENV GVER=2013-08-31
ENV GPATH=/barleygraph/gmap-${GVER}/local
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-${GVER}.tar.gz
RUN tar xfz gmap-gsnap-${GVER}.tar.gz && cd gmap-${GVER} && \
    ./configure --prefix=${GPATH} CFLAGS="-O3 -Dstat64=stat -Wno-implicit-function-declaration" && \
    make && make install && make clean
RUN rm gmap-gsnap-${GVER}.tar.gz

# align2graph 
RUN wget https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/refs/heads/main/align2graph.py
RUN wget https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/refs/heads/main/version.py
RUN chmod +x align2graph.py
RUN mv align2graph.py align2graph

# scripts and config files from GitHub repository
RUN wget -qO gmap_build.sh https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/gmap_build.sh
RUN wget -qO build_imputation_index.sh https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/build_imputation_index.sh
RUN wget -qO imputation.pl https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/imputation.pl
RUN wget -qO setup_graph.pl https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/main/docker/scripts/setup_graph.pl
RUN wget -qO Pan20-mmap-pro.yaml https://github.com/eead-csic-compbio/barleygraph/blob/main/graphs/Pan20-mmap-pro.yaml

# Ensure correct conda paths for anaconda/miniconda image are fully exposed globally
ENV PATH="/opt/miniconda3/bin:/opt/miniconda3/envs/phgv2.4/bin:${PATH}"

# Reset default shell behavior to prevent build fails
SHELL ["/bin/bash", "-c"]

# Install and get phgtools utilities from Github using conda run
RUN conda run -n phgv2.4 pip install git+https://github.com/jsarriaa/PHGv2Tools.git
# Set it in path
RUN echo '#!/bin/bash' > /barleygraph/phgtools && \
    echo 'exec /opt/miniconda3/envs/phgv2.4/bin/phgtools "$@"' >> /barleygraph/phgtools && \
    chmod +x /barleygraph/phgtools
# Ensure writeable cache directory for phgtools
ENV MPLCONFIGDIR=/tmp/matplotlib_cache
ENV XDG_CACHE_HOME=/tmp/xdg_cache

# Ensure scripts are executable
RUN chmod +x gmap_build.sh build_imputation_index.sh imputation.pl setup_graph.pl

# Remove .sh extensions from shell scripts and rename setup_databases.sh to setup_database
RUN mv build_imputation_index.sh build_imputation_index && \
    mv imputation.pl imputation && \
    mv setup_graph.pl setup_graph

# Create wrapper scripts that safely invoke the conda environment via conda run
RUN echo '#!/bin/bash' > /barleygraph/haplopainting && \
    echo 'exec conda run -n phgv2.4 python3 -m phgtools.modules.haplopainting "$@"' >> /barleygraph/haplopainting && \
    chmod +x /barleygraph/haplopainting

RUN echo '#!/bin/bash' > /barleygraph/hvcf2bed && \
    echo 'exec conda run -n phgv2.4 python3 -m phgtools.modules.hvcf2bed "$@"' >> /barleygraph/hvcf2bed && \
    chmod +x /barleygraph/hvcf2bed

# Clean up conda package cache to reduce image size
#RUN rm -rf /opt/miniconda3/pkgs/*

ENV PATH="/barleygraph:${GPATH}/bin:${PHGPATH}:${PATH}"
