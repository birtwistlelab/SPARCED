# must be built from root directory using the -f directive
# BUILD: docker build -t sparced-notebook -f docker/jupyter.Dockerfile .
# TEST LOCAL: docker run -p 8888:8888 --name testnb1 -i -t sparced-notebook
# TAG: docker tag sparced-notebook ebensma/sparced-notebook:latest
# PUSH: docker push ebensma/sparced-notebook:latest

FROM ubuntu:18.04

# install system dependencies
RUN apt-get update -qq && apt-get install -qq -y curl git python3-dev python3-pip libhdf5-serial-dev libatlas-base-dev swig rsync && ln -s /usr/bin/python3 python

#changing working directory in Docker container
WORKDIR /app

# copy data from local into Docker container
ADD . /app/

# install python dependencies
RUN pip3 install -r requirements.txt
RUN pip3 install jupyter

# open directory as jupyter notebook
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

# to Run: docker run -p 8888:8888 myaccount/new_project

