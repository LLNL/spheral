#!/usr/bin/env bash

BUILD_SUFFIX=lc_toss3-clang-9.0.0-release-seq-cxx

rm -rf ${BUILD_SUFFIX} 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load clang/9.0.0

cmake \
  ../.. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/clang/clang-9.0.0/bin/clang++ \
  -DCMAKE_C_COMPILER=/usr/tce/packages/clang/clang-9.0.0/bin/clang \
  -DENABLE_OPENMP=On \
  -DENABLE_MPI=Off \
  -DCMAKE_INSTALL_PREFIX=../install \
  -DSPHERAL_TPL_DIR=/usr/workspace/wsrzd/davis291/SPHERAL/toss_Spheral_clang9/install \
  -DBUILD_TPL=Off \
  -DENABLE_STATIC_CXXONLY=On
  "$@" \
