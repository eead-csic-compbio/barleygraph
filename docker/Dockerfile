FROM continuumio/miniconda3

LABEL authors="Joan Sarria, Bruno Contreras Moreira"
LABEL description="Container with align2graph, PHG and selected pangenome graph [Med13, Pan20]"
LABEL license="GPLv3"
LABEL maintainers="compbio@eead.csic.es"

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

# PHG and python deps
RUN wget https://raw.githubusercontent.com/maize-genetics/phg_v2/refs/heads/main/src/main/resources/phg_environment.yml
RUN echo "  - pyyaml" >> phg_environment.yml
RUN conda env create -f phg_environment.yml
RUN echo "conda activate phgv2-conda" >> ~/.bashrc
ENV PATH /opt/conda/envs/phgv2-conda/bin:$PATH
ENV CONDA_DEFAULT_ENV phgv2-conda

# GMAP 2013-08-31 makes indexes smaller than version from PHG
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2013-08-31.tar.gz
RUN tar xfz gmap-gsnap-2013-08-31.tar.gz && cd gmap-2013-08-31 && ./configure --prefix=/barleygraph/gmap-2013-08-31/local && make && make install && make clean
RUN rm gmap-gsnap-2013-08-31.tar.gz

# align2graph 
RUN wget https://raw.githubusercontent.com/eead-csic-compbio/eead-csic-compbio.github.io/refs/heads/master/scripts/barleygraph/align2graph.py
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

ENV PATH="/barleygraph:/barleygraph/gmap-2013-08-31/local/bin:${PATH}"
SHELL ["conda", "run", "--no-capture-output", "-n", "phgv2-conda", "/bin/bash", "-c"]

#################################  Useful command lines  ################################# 
#
# 1) Build container
#
# export VERSION=`date '+%Y%m%d'`
# export DREPO=eeadcsiccompbio/barleygraph
# docker build --build-arg graph=Med13 --tag ${DREPO}:Med13-${VERSION} --tag ${DREPO}:Med13-latest --force-rm --compress .
# or
# docker build --build-arg graph=Pan20 --tag ${DREPO}:Pan20-${VERSION} --tag ${DREPO}:Pan20-latest --force-rm --compress .
# 
# 2) Create local folder for GMAP indices, outside the container, as these are bulky
#
# mkdir -m 777 -p ~/local_gmap_db
#  
# 3) Launch container:
# 3.1) make GMAP index if required
# docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest index
#
# 3.2) map FASTA sequences
# docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest align2graph test.fa
