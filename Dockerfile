#To Build: docker build -t sparced .
#To Run: docker run -p 8888:8888 --name testcp1 -it sparced
#To relaunch: docker start -i testcp1


# Ubuntu and conda base images:
FROM ubuntu:18.04

# Install system-level dependencies
RUN apt-get update -qq && apt-get install -qq -y curl git python3-dev python3-pip libhdf5-serial-dev libatlas-base-dev swig rsync && ln -s /usr/bin/python3 python

# Changing working directory in Docker container
WORKDIR /SPARCED

# Copy data from local into Docker container
ADD . /SPARCED

# Install Miniconda
RUN curl -o miniconda.sh -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    chmod +x miniconda.sh && \
    ./miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Add conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Create and activate the Conda environment from environment.yml
RUN conda env create -f environment.yml
SHELL ["conda", "run", "-n", "sparced_cellPop", "/bin/bash", "-c"]

RUN echo "source activate $(head -1 environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 environment.yml | cut -d' ' -f2)/bin:$PATH

# Now, install Python packages within the Conda environment
RUN pip install --no-cache-dir amici==0.16.1
RUN pip install --no-cache-dir antimony==2.12.0.2
RUN pip install --no-cache-dir numpy==1.23.2

# Open directory
WORKDIR ./scripts
RUN python3 createModel_hpc.py
