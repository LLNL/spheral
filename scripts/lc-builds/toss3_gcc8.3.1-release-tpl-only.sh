#!/usr/bin/env bash

#-----------------------------------------------------------------------------------
# ***********************
# Handle Named Arguments.
# ***********************

SRC_DIR=
CMAKE_ARGS=
REL_MOD=
while getopts "d:D:C:" name;
do
  case $name in
  d)  SRC_DIR="$OPTARG"
      printf 'Option -d Spheral source  "%s" specified\n' "$SRC_DIR"
      if [ -d $SRC_DIR ]
      then
        printf 'Cannot find "%s": Assuming relative path ...\n' "$SRC_DIR"
        REL_MOD=../../
      fi
      ;;
  D)  CMAKE_ARGS=$CMAKE_ARGS"-D$OPTARG ";;
  C)  CMAKE_ARGS=$CMAKE_ARGS"-C $OPTARG ";;
  ?)  printf "Usage: %s: [-d Spheral source directory] [CMake args ...]\n" $0
        ;;
  esac
done

if [ ! -z "$CMAKE_ARGS" ]; then
  printf "Additional CMake arguments are: %s\n" "$CMAKE_ARGS"
fi
shift $(($OPTIND - 1))
#-----------------------------------------------------------------------------------

BUILD_SUFFIX=lc_toss3-gcc-8.3.1-rel-tpl

rm -rf ${BUILD_SUFFIX} 2>/dev/null
mkdir -p ${BUILD_SUFFIX}/install
mkdir -p ${BUILD_SUFFIX}/build && cd ${BUILD_SUFFIX}/build

module load cmake/3.14.5
module load gcc/8.3.1

cmake \
  ${REL_MOD}${SRC_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.3.1/bin/g++ \
  -DENABLE_OPENMP=On \
  -DENABLE_MPI=On \
  -DCMAKE_INSTALL_PREFIX=${REL_MOD}${SRC_DIR}stall \
  -DBUILD_TPL_ONLY=On \
  $CMAKE_ARGS \
