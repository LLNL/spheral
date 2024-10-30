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
FROM ubuntu:20.04 AS spheral-build-env-local

ARG SPEC=gcc
ARG HOST_CONFIG=docker-$SPEC

# Update Ubuntu and install necessary packages.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y build-essential git gfortran mpich autotools-dev autoconf sqlite pkg-config uuid gettext cmake libncurses-dev libgdbm-dev libffi-dev libssl-dev libexpat-dev libreadline-dev libbz2-dev locales python python3 unzip libtool wget curl tk-dev

# Setup system locale for pip package encoding/decoding 
RUN locale-gen en_US.UTF-8

# Set up TPLs for SPEC
WORKDIR /home/spheral/workspace/
COPY scripts scripts
COPY .uberenv_config.json .
RUN python3 scripts/devtools/tpl-manager.py --spec $SPEC --spheral-spack-dir /home

# Clean workspace once dependencies are installed
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

# Copy Spheral source and generate host config from tpl-manager (all dependencies should already be installed).
COPY . .
RUN python3 scripts/devtools/tpl-manager.py --spec $SPEC --upstream-dir /home/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder_ --spack-url /home/spack

# Configure Spheral with SPEC TPLs.
RUN mv *.cmake $HOST_CONFIG.cmake
RUN python3 scripts/devtools/host-config-build.py --host-config $HOST_CONFIG.cmake

# Build Spheral
WORKDIR build_$HOST_CONFIG/build
RUN make -j $JCXX Spheral_CXX
RUN make -j $JPY
RUN make install

# Run ATS testing suite.
WORKDIR ../install
ENV MPLBACKEND=agg
RUN ./spheral-ats --level 99 --mpiexec /usr/bin/mpirun --npMax $JCXX tests/integration.ats
# -----------------------------------------------------------------------------
