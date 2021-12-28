################# BASE IMAGE #####################
FROM mambaorg/micromamba:0.15.3
##site to test docker configuration files
################## METADATA #######################
LABEL base_image="mambaorg/micromamba"
LABEL version="0.15.3"
LABEL software="k-count-nf"
LABEL software.version="1.1"
LABEL about.summary="Container image containing all requirements for k-count-nf"
LABEL about.home="https://github.com/digenoma-lab/k-count-nf"
LABEL about.documentation="https://github.com/digenoma-lab/k-count-nf/README.md"
LABEL about.license_file="https://github.com/digenoma-lab/k-count-nf/LICENSE.txt"
LABEL about.license="MIT"

################## MAINTAINER ######################
MAINTAINER Alex Di Genova <digenova@gmail.com>
################## INSTALLATION ######################
USER root
#the next command install the ps command needed by nexflow to collect run metrics
RUN apt-get update && apt-get install -y procps
USER micromamba
COPY --chown=micromamba:micromamba environment.yml /tmp/environment.yml
RUN micromamba create -y -n k-count -f /tmp/environment.yml && \
    micromamba clean --all --yes
ENV PATH /opt/conda/envs/k-count/bin:$PATH

