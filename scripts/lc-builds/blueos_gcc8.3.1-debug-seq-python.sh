#!/usr/bin/env bash

SCRIPT_PATH=${0%/*}
. "$SCRIPT_PATH/utils/parse-args.sh"

# Inherit build directory name from script name
BUILD_SUFFIX="lc_$(TMP=${BASH_SOURCE##*/}; echo ${TMP%.*})"

rm -rf ${BUILD_SUFFIX} 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load gcc/8.3.1

cmake \
  ../.. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.3.1/bin/g++ \
  -C ${HOST_CONFIGS_DIR}/lc-builds/blueos/gcc8.3.1_tpl.cmake \
  -DBLT_CXX_STD=c++14 \
  -DENABLE_OPENMP=On \
  -DENABLE_MPI=Off \
  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  "$@" \
  #-DSPHERAL_TPL_DIR=/usr/workspace/wsrzd/davis291/SPHERAL/blueos_Spheral_gcc8_noMPI/install/tpl \