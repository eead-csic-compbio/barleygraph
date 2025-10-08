FROM continuumio/miniconda3

LABEL authors="Joan Sarria, Bruno Contreras Moreira"
LABEL description="Container with align2graph, PHG and selected pangenome graph [Med13, Pan20]"
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
	bedtools \
	&& apt clean && apt purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# PHG and python deps; other versions of GMAP and agc are installed later
ENV PHGYAML https://raw.githubusercontent.com/maize-genetics/phg_v2/refs/heads/main/src/main/resources/phg_environment.yml
RUN wget ${PHGYAML}
RUN echo "  - pyyaml" >> phg_environment.yml
RUN conda env create -f phg_environment.yml
RUN echo "conda activate phgv2-conda" >> ~/.bashrc
ENV PATH /opt/conda/envs/phgv2-conda/bin:$PATH
ENV CONDA_DEFAULT_ENV phgv2-conda

# GMAP 2013-08-31 makes indexes smaller than version from PHG
ENV GVER 2013-08-31
ENV GPATH /barleygraph/gmap-${GVER}/local
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-${GVER}.tar.gz
RUN tar xfz gmap-gsnap-${GVER}.tar.gz && cd gmap-${GVER} && ./configure --prefix=${GPATH} && make && make install && make clean
RUN rm gmap-gsnap-${GVER}.tar.gz

# compatible agc 3.0,noavx runs on old linux box, 3.1 does not
ENV AGCVER 3.0
ENV AGCTAR agc-${AGCVER}_x64-linux-noavx.tar.gz
ENV AGCURL https://github.com/refresh-bio/agc/releases/download/v${AGCVER}/${AGCTAR}
ENV AGCPATH /barleygraph/agc-${AGCVER}_x64-linux-noavx/
RUN wget ${AGCURL}
RUN tar xfz ${AGCTAR}
RUN rm ${AGCTAR}

# align2graph 
RUN wget https://raw.githubusercontent.com/eead-csic-compbio/barleygraph/refs/heads/main/align2graph.py
RUN chmod +x align2graph.py
RUN mv align2graph.py align2graph

ARG graph
LABEL graph=$graph
RUN echo $graph

# graphs, config and bits
COPY ${graph}.yaml graph.yaml
ADD ${graph}.tgz $graph
COPY test.fa test.fa

# scripts
COPY gmap_build.sh gmap_build.sh
COPY index_${graph} index

ENV PATH="/barleygraph:${AGCPATH}:${GPATH}/bin:${PATH}"
SHELL ["conda", "run", "--no-capture-output", "-n", "phgv2-conda", "/bin/bash", "-c"]
