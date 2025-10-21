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


############# SYSTEM DEPENDENCIES #################
# Install Java 17 (JRE+JDK) and curl before Nextflow
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update --allow-releaseinfo-change && \
    apt-get -y upgrade && \
    apt-get -y install --no-install-recommends \
        tzdata \
        openjdk-17-jre \
        openjdk-17-jdk \
        curl \
        procps && \
    java -version && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
ENV DEBIAN_FRONTEND=dialog

############ NEXTFLOW INSTALLATION ################
# Install Nextflow in a user-local directory
ENV HOME=/root
RUN mkdir -p $HOME/.local/bin && \
    curl -s https://get.nextflow.io | bash && \
    mv nextflow $HOME/.local/bin/ && \
    chmod +x $HOME/.local/bin/nextflow && \
    $HOME/.local/bin/nextflow -version
ENV PATH="$HOME/.local/bin:${PATH}"

############# MSA ENV INSTALLATION ################
COPY environment.yml /tmp/environment.yml
RUN mamba env create -n MSA -f /tmp/environment.yml && \
    conda clean -afy && \
    rm -rf /root/.cache
# Auto-activate environment when starting container
RUN echo "conda activate MSA" >> ~/.bashrc
ENV PATH=/opt/conda/envs/MSA/bin:$PATH

################## WORKING DIRECTORY ##############
WORKDIR /workspace

################## DEFAULT SHELL ##################
SHELL ["/bin/bash", "--login", "-c"]

