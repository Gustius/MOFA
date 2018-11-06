FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/mofapipeline pipeline"
MAINTAINER Gust Verstraete <gust.verstraete@ugent.be>

COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-mofapipeline-1.0dev/bin:$PATH
