# should be built from the SPARCED/docker directory (avoids adding excess files)
# BUILD: docker build -t sparced -f docker/nf.Dockerfile .
# TAG: docker tag sparced ebensma/sparced:latest
# PUSH: docker push ebensma/sparced:latest

FROM ubuntu:18.04

# install system dependencies
RUN apt-get update -qq && apt-get install -qq -y curl git python3-dev python3-pip libhdf5-serial-dev libatlas-base-dev swig rsync && ln -s /usr/bin/python3 python

#changing working directory in Docker container
WORKDIR /app

# copy data from local into Docker container
ADD . /app/

# install python dependencies
RUN pip3 install -r requirements.txt
