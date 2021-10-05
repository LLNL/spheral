#!/usr/bin/env bash

SCRIPT_PATH=${0%/*}
. "$SCRIPT_PATH/utils/parse-args.sh"

# Inherit build directory name from script name
BUILD_SUFFIX="lc_$(TMP=${BASH_SOURCE##*/}; echo ${TMP%.*})"

rm -rf ${BUILD_SUFFIX}/build 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load gcc/8.1.0

cmake \
  ${SRC_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++ \
  -DENABLE_OPENMP=On \
  -DENABLE_MPI=On \
  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DBUILD_TPL_ONLY=On \
  -Dhdf5_BUILD=Off \
  -Dhdf5_DIR="/usr/gapps/Spheral/tpl/$SYS_TYPE/lchdf5" \
  -Dsilo_BUILD=Off \
  -Dsilo_DIR="/usr/gapps/Spheral/tpl/$SYS_TYPE/lcsilo" \
  -Wno-dev \
  $CMAKE_ARGS

make -j install

cd -
find ${BUILD_SUFFIX}/ -type d -exec chmod g+rx {} \;
find ${BUILD_SUFFIX}/ -type f -exec chmod g+rx {} \;
find ${BUILD_SUFFIX}/ -name "*egg-info" -exec chgrp wciuser {} \;
