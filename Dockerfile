FROM debian:stretch-slim
COPY . /annotate
ENV PATH="/annotate:${PATH}"
WORKDIR /cwd
