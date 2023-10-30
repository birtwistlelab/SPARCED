#To Run: docker run -p 8888:8888 --name testcp1 -it sparced
#To relaunch: docker start -i testcp1

# Ubuntu and conda base images:
FROM ubuntu:22.04

# Install system-level dependencies
RUN apt-get update -qq && apt-get install -qq -y curl git python3-dev \
python3-pip build-essential libhdf5-serial-dev libatlas-base-dev \
libopenblas-dev openmpi-bin openmpi-doc libopenmpi-dev \
libxml2 swig rsync && ln -s /usr/bin/python3 python

# Changing working directory in Docker container
WORKDIR /SPARCED

# Copy data from local into Docker container
ADD . /SPARCED/

#install python dependencies
RUN pip3 install -r requirements.txt

# Open directory
WORKDIR ./scripts
RUN python3 createModel_hpc.py
