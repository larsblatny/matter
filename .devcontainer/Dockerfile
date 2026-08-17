# Use Ubuntu 24.04 as the base image
FROM ubuntu:24.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install required dependencies
RUN apt-get update && apt-get install -y \
    cmake \
    git \
    build-essential \
    libeigen3-dev \
    libopenvdb-dev \
    libtbb-dev \
    libboost-all-dev \
    libilmbase-dev \
    libopenexr-dev \
    python3-dev \
    python3-pytest \
    python3-numpy

# Set working directory
WORKDIR /workspace
