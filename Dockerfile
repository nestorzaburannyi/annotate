FROM debian:stretch-slim
RUN apt-get update && apt-get install -y libidn11 wget
COPY . /annotate
ENV PATH="/annotate:${PATH}"
WORKDIR /cwd
