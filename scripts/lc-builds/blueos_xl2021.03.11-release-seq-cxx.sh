#!/usr/bin/env bash

SCRIPT_PATH=${0%/*}
. "$SCRIPT_PATH/utils/parse-args.sh"

# Inherit build directory name from script name
BUILD_SUFFIX="lc_$(TMP=${BASH_SOURCE##*/}; echo ${TMP%.*})"

rm -rf ${BUILD_SUFFIX} 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load xl/2021.03.11

cmake \
  ${SRC_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/xl/xl-2021.03.11/bin/xlc++ \
  -DCMAKE_Fortran_COMPILER=/usr/tce/packages/xl/xl-2021.03.11/bin/xlf \
  -DENABLE_OPENMP=Off \
  -DENABLE_MPI=Off \
  -DENABLE_ANEOS=Off \
  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DENABLE_STATIC_CXXONLY=On \
  $CMAKE_ARGS \
  #-C ${HOST_CONFIGS_DIR}/lc-builds/toss3/gcc8.3.1_tpl.cmake \
