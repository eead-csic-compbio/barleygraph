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
	&& apt clean && apt purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# PHG
RUN wget https://github.com/maize-genetics/phg_v2/releases/download/2.4.69.224/PHGv2-v2.4.tar
RUN tar xf PHGv2-v2.4.tar
RUN rm PHGv2-v2.4.tar
ENV PHGPATH /barleygraph/phg/bin

# PHG and python deps
COPY environment.yml environment.yml
RUN conda env create -f environment.yml
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

ENV PATH="/barleygraph:${GPATH}/bin:${PHGPATH}:${PATH}"
SHELL ["conda", "run", "--no-capture-output", "-n", "phgv2.4", "/bin/bash", "-c"]
