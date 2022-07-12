FROM ubuntu:20.04

# Update Ubuntu and install necessary packages.
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update -y
RUN apt upgrade -y
RUN apt install -y build-essential git gfortran mpich autotools-dev autoconf sqlite pkg-config uuid gettext cmake libncurses-dev libgdbm-dev libffi-dev libssl-dev libexpat-dev libreadline-dev python3 unzip libtool wget

# Initialize script ARGs.
# SPEC can be defined on the command line with --build-args SPEC=<...>
ARG SPEC=gcc
ARG HOST_CONFIG=docker-$SPEC

# Set up TPLs for SPEC
WORKDIR /home/spheral/workspace/
COPY scripts scripts
COPY .uberenv_config.json .
RUN python3 scripts/devtools/tpl-manager.py --spec $SPEC

# Configure Spheral with SPEC TPLs.
RUN mv *.cmake $HOST_CONFIG.cmake
COPY . .
RUN python3 scripts/devtools/host-config-build.py --host-config $HOST_CONFIG.cmake

# Build Spheral
WORKDIR build_$HOST_CONFIG/build
RUN make -j 8 Spheral_CXX
RUN make -j 4 install

# Run ATS testing suite.
WORKDIR ../install
ENV MPLBACKEND=agg
RUN ./spheral-atstest --filter="level<100" tests/integration.ats
