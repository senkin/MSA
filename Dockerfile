################## BASE IMAGE #####################
FROM condaforge/mambaforge:24.9.2-0

################## METADATA #######################

LABEL base_image="condaforge/mambaforge"
LABEL version="24.9.2-0"
LABEL software="MSA"
LABEL software.version="3.0"
LABEL about.summary="Container image containing all requirements for MSA nextflow pipeline"
LABEL about.home="https://gitlab.com/s.senkin/MSA/"
LABEL about.documentation="https://gitlab.com/s.senkin/MSA/README.md"
LABEL about.license_file="https://gitlab.com/s.senkin/MSA/LICENSE.txt"
LABEL about.license="GNU-3.0"
LABEL maintainer="s.senkin@gmail.com"

################## INSTALLATION ######################
COPY environment.yml /tmp/environment.yml
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    mamba env create -n MSA -f /tmp/environment.yml && \
    conda clean -afy && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /root/.cache
RUN echo "conda activate MSA" >> ~/.bashrc
ENV PATH=/opt/conda/envs/MSA/bin:$PATH
WORKDIR /workspace

