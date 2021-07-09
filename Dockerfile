###############################################################################
# Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

FROM davis291/spheral:gcc8-tpl AS gcc8
COPY --chown=axom:axom . /home/axom/workspace
WORKDIR /home/axom/workspace

RUN mkdir -p spheral_release/build spheral_release/install
RUN cd spheral_release/build && cmake ../../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran-8 -DCMAKE_INSTALL_PREFIX="../install" -DENABLE_MPI=Off -DBUILD_TPL=Off -DSPHERAL_TPL_DIR=/home/axom/workspace/spheral_tpl/install/tpl
RUN cd spheral_release/build && make -j4 install
RUN ./spheral_release/install/spheral -c "import Spheral"
