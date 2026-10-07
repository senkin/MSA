################## BASE IMAGE #####################
# pinned by digest (multi-arch index of the 24.9.2-0 tag), so rebuilds start from the same base
FROM condaforge/mambaforge:24.9.2-0@sha256:3b5f55c0b2f447cf94b8f67ca9094f9d2a9d316a4506fd5aa915faaf37ae29d4

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
# Install Java 17 (JRE+JDK; Nextflow 26.04 runs on Java 17 to 26) and curl before Nextflow
ENV DEBIAN_FRONTEND=noninteractive
# (no 'apt-get upgrade': it would make every rebuild pull whatever updates exist that day)
RUN apt-get -y update --allow-releaseinfo-change && \
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
# Install a fixed Nextflow release (self-contained -dist build, nothing downloaded at run
# time) in a user-local directory
ARG NXF_VER=26.04.6
ENV HOME=/root
ENV NXF_VER=${NXF_VER}
RUN mkdir -p $HOME/.local/bin && \
    curl -fsSL -o $HOME/.local/bin/nextflow \
        https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow-${NXF_VER}-dist && \
    chmod +x $HOME/.local/bin/nextflow && \
    $HOME/.local/bin/nextflow -version
ENV PATH="$HOME/.local/bin:${PATH}"

############# MSA ENV INSTALLATION ################
# Exact lock of environment.yml for the image's architecture (see conda-lock/)
ARG TARGETARCH
COPY conda-lock/environment.linux-64.txt conda-lock/environment.linux-aarch64.txt /tmp/conda-lock/
RUN case "${TARGETARCH:-amd64}" in \
        amd64) platform=linux-64 ;; \
        arm64) platform=linux-aarch64 ;; \
        *) echo "unsupported architecture: ${TARGETARCH}" >&2; exit 1 ;; \
    esac && \
    mamba create -y -n MSA --file /tmp/conda-lock/environment.${platform}.txt && \
    conda clean -afy && \
    rm -rf /root/.cache /tmp/conda-lock
# Auto-activate environment when starting container
RUN echo "conda activate MSA" >> ~/.bashrc
ENV PATH=/opt/conda/envs/MSA/bin:$PATH

################## WORKING DIRECTORY ##############
WORKDIR /workspace

################## DEFAULT SHELL ##################
SHELL ["/bin/bash", "--login", "-c"]

