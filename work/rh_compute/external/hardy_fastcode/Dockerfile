FROM debian:bookworm-slim

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ca-certificates \
        gfortran \
        libopenmpi-dev \
        libpari-dev \
        make \
        openmpi-bin \
        patch \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
