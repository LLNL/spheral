FROM ubuntu:20.04 AS spheral-build-env-local

# Update Ubuntu and install necessary packages.
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y build-essential git gfortran mpich autotools-dev autoconf sqlite pkg-config uuid gettext cmake libncurses-dev libgdbm-dev libffi-dev libssl-dev libexpat-dev libreadline-dev python python3 unzip libtool wget curl

RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python2.7 get-pip.py
RUN rm get-pip.py

# Initialize script ARGs.
# SPEC can be defined on the command line with --build-args SPEC=<...>
ARG SPEC=gcc
ARG HOST_CONFIG=docker-$SPEC

# Set up TPLs for SPEC
WORKDIR /home/spheral/workspace/
COPY scripts scripts
COPY .uberenv_config.json .
COPY cmake/tpl/util/requirements.txt cmake/tpl/util/requirements.txt
RUN python3 scripts/devtools/tpl-manager.py --spec $SPEC --spheral-spack-dir /home
RUN rm -rf /home/spheral/workspace/*


FROM spheral-build-env AS spheral

ARG SPEC=gcc
ARG HOST_CONFIG=docker-$SPEC

COPY . .
RUN python3 scripts/devtools/tpl-manager.py --spec $SPEC --upstream-dir /home/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder_

# Configure Spheral with SPEC TPLs.
RUN mv *.cmake $HOST_CONFIG.cmake
RUN python3 scripts/devtools/host-config-build.py --host-config $HOST_CONFIG.cmake

# Build Spheral
WORKDIR build_$HOST_CONFIG/build
RUN make -j 8 Spheral_CXX
RUN make
#RUN make install
#
## Run ATS testing suite.
#WORKDIR ../install
#ENV MPLBACKEND=agg
#RUN ./spheral-atstest --filter="level<100" tests/integration.ats
