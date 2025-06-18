# To build and spheral-build-env:
#   sudo env DOCKERBUILDKIT=1 docker build . --target spheral-build-env-local --tag spheral-build-env (--progress=plain)
#   Optional Arguments:
#     --progress=plain        : Prints plain output to terminal instead of windowed version.
#     --build-args SPEC=...   : Specify optional build argument to override. Default = gcc
#                               e.g. --build-args SPEC=clang 

# To build and run a spheral test:
#   sudo env DOCKERBUILDKIT=1 docker build . --target spheral --tag spheral (--progress=plain) (--network none)
#   Optional Arguments:
#     --progress=plain        : Prints plain output to terminal instead of windowed version.
#     --network none          : Simulate a build and run on a system with installed dependencies but no network.
#     --build-args SPEC=...   : Specify the SPEC (must be the same as the spec for spheral-build-env) Default = gcc
#                  JCXX=N     : Specify number of threads to build C++ portion of Spheral build. Default = 8
#                  JPY=N      : Specify number of threads to build python portion of Spheral build. Default = 1

# To use an interactive terminal with built spheral-build-env or spheral img :
#   sudo docker run -it (spheral-build-env/spheral)



# -----------------------------------------------------------------------------
# SPHERAL-BUILD-ENV
# -----------------------------------------------------------------------------
FROM ubuntu:24.04 AS spheral-build-env-local

ARG SPEC=gcc
ARG HOST_CONFIG=docker-$SPEC

# Update Ubuntu and install necessary packages.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y
RUN apt-get upgrade -y

# Required packages
RUN apt-get install -y bzip2 ca-certificates g++ gcc gfortran git gzip
RUN apt-get install -y lsb-release patch python3 tar unzip xz-utils zstd
RUN apt-get install -y libtool curl wget libcurl4-openssl-dev tk-dev autotools-dev
RUN apt-get install -y build-essential python3-dev python3-pip python3-venv

# Recommended packages (MPICH library is broken for 22.04/24.04, use openmpi)
RUN apt-get install -y cmake autoconf automake libopenmpi-dev libreadline-dev
RUN apt-get install -y iputils-ping

# Setup system locale for pip package encoding/decoding 
RUN locale-gen en_US.UTF-8

# Set up TPLs for SPEC
WORKDIR /home/spheral/workspace/
COPY scripts scripts


RUN python3 scripts/devtools/tpl-manager.py --spec spheral%$SPEC --spack-dir /home

COPY . .

# Configure Spheral with SPEC TPLs.
RUN mv *.cmake $HOST_CONFIG.cmake
RUN python3 scripts/devtools/host-config-build.py --host-config $HOST_CONFIG.cmake

# First time install of Spheral pip dependencies
WORKDIR build_$HOST_CONFIG/build
RUN make python_build_env
RUN make python_runtime_env

# Clean workspace once dependencies are installed
WORKDIR /home/spheral/workspace/

RUN rm -rf /home/spheral/workspace/*
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# SPHERAL BUILD & TEST
# -----------------------------------------------------------------------------
FROM spheral-build-env AS spheral

ARG SPEC=gcc
ARG HOST_CONFIG=docker-$SPEC
ARG JCXX=8
ARG JPY=1

WORKDIR /home/spheral/workspace/

# Copy Spheral source and generate host config from tpl-manager (all dependencies should already be installed).
COPY . .
RUN python3 scripts/devtools/tpl-manager.py --spec spheral%$SPEC --spack-dir /home

# Configure Spheral with SPEC TPLs.
RUN mv *.cmake $HOST_CONFIG.cmake
RUN python3 scripts/devtools/host-config-build.py --host-config $HOST_CONFIG.cmake -DSPHERAL_NETWORK_CONNECTED=Off

# Build Spheral
WORKDIR build_$HOST_CONFIG/build
RUN make python_build_env
RUN make python_runtime_env
RUN make -j $JCXX Spheral_CXX
RUN make -j $JPY
RUN make install

# Run ATS testing suite.
WORKDIR ../install

# Instead, we will just run it normally
RUN ./spheral-ats --level 99 tests/integration.ats
# -----------------------------------------------------------------------------
