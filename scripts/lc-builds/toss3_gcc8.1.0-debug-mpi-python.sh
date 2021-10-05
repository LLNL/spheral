#!/usr/bin/env bash

SCRIPT_PATH=${0%/*}
. "$SCRIPT_PATH/utils/parse-args.sh"

# Inherit build directory name from script name
BUILD_SUFFIX="lc_$(TMP=${BASH_SOURCE##*/}; echo ${TMP%.*})"

rm -rf ${BUILD_SUFFIX} 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load gcc/8.1.0

cmake \
  ${SRC_DIR} \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++ \
  -C ${HOST_CONFIGS_DIR}/lc-builds/toss3/gcc8.1.0_tpl.cmake \
  -Dhdf5_BUILD=Off \
  -Dhdf5_DIR="/usr/gapps/Spheral/tpl/$SYS_TYPE/lchdf5" \
  -Dsilo_BUILD=Off \
  -Dsilo_DIR="/usr/gapps/Spheral/tpl/$SYS_TYPE/lcsilo" \
  -DENABLE_OPENMP=On \
  -DENABLE_MPI=On \
  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  $CMAKE_ARGS \
