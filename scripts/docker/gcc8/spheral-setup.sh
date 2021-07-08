#!/bin/sh
set -e
set -x

apt-get -qq update
apt-get -qq install -y g++-8 gfortran-8 zlib1g-dev libssl-dev

apt remove cmake
pip install cmake --upgrade
