FROM debian:stretch-slim
RUN apt-get update && apt-get install -y libidn11 wget libgfortran3 python-requests openmpi-bin python-numpy python3 libxml-simple-perl openssh-client python-fastcluster
COPY . /annotate
ENV PATH="/annotate:${PATH}"
WORKDIR /pwd
