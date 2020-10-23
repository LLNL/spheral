#!/usr/bin/env bash

BUILD_SUFFIX=lc_toss3-gcc-8.3.1-debug-seq-cxx

rm -rf ${BUILD_SUFFIX} 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load gcc/8.3.1

cmake \
  ../.. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.3.1/bin/g++ \
  -DENABLE_OPENMP=On \
  -DENABLE_MPI=On \
  -DCMAKE_INSTALL_PREFIX=../install \
  -DSPHERAL_TPL_DIR=/usr/workspace/wsrzd/davis291/SPHERAL/toss_Spheral_gcc8/install3/tpl \
  -DBUILD_TPL=Off \
  -DENABLE_STATIC_CXXONLY=On
  "$@" \
