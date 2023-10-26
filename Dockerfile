#To Run: docker run -p 8888:8888 --name testcp1 -it sparced
#To relaunch: docker start -i testcp1

# Ubuntu and conda base images:
FROM ubuntu:22.04
FROM continuumio/miniconda3

# Install system-level dependencies
RUN apt-get update -qq && apt-get install -qq -y curl git python3-dev python3-pip libhdf5-serial-dev libatlas-base-dev libxml2 swig rsync && ln -s /usr/bin/python3 python

# Changing working directory in Docker container
WORKDIR /SPARCED

# Copy data from local into Docker container
ADD . /SPARCED

# Add conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Add .condarc location to path
ENV CONDARC /SPARCED/bin/.condarc

# Create and activate the Conda environment from environment.yml
RUN conda env create -f bin/environment.yml && \
conda init bash

# Set the shell interpreter to use bash instead of /bin/sh
SHELL ["/bin/bash", "-c"]

RUN echo "source activate $(head -1 environment.yml | cut -d' ' -f2)" > ~/.bashrc

ENV PATH /opt/conda/envs/$(head -1 environment.yml | cut -d' ' -f2)/bin:$PATH

# Open directory
WORKDIR ./scripts
RUN source activate sparced && \
python3 createModel_hpc.py
