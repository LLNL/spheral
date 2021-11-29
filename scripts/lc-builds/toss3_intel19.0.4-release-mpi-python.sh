#!/usr/bin/env bash

BUILD_SUFFIX=lc_toss3-intel-19.0.4-rel-mpi-py

rm -rf ${BUILD_SUFFIX} 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load intel/19.0.4

cmake \
  ../.. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/intel/intel-19.0.4/bin/icpc \
  -DCMAKE_C_COMPILER=/usr/tce/packages/intel/intel-19.0.4/bin/icc \
  -DENABLE_OPENMP=On \
  -DENABLE_MPI=On \
  -DCMAKE_INSTALL_PREFIX=../install \
  -DSPHERAL_TPL_DIR=/usr/workspace/wsrzd/davis291/SPHERAL/toss_Spheral_intel19.0.4/install \
  -DBUILD_TPL=Off \
  "$@" \
