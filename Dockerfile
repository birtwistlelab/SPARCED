FROM ubuntu:18.04

# install system dependencies
RUN apt-get update -qq && apt-get install -qq -y curl git python3-dev python3-pip libhdf5-serial-dev libatlas-base-dev swig && ln -s /usr/bin/python3 python

#changing working directory in Docker container
WORKDIR /app

# copy data from local into Docker container
COPY . /app/

# install python dependencies
RUN pip3 install -r requirements.txt

# install system dependencies
RUN apt-get install libatlas-base-dev
RUN apt-get install libhdf5-serial-dev
RUN apt-get install swig

