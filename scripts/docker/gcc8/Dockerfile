###############################################################################
# Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

#FROM axom/compilers:ubuntu-18-base AS gcc8
FROM axom/compilers:gcc-8 AS gcc-8
ENV GTEST_COLOR=1

USER root
ADD spheral-setup.sh /root/spheral-setup.sh
RUN /root/spheral-setup.sh

USER axom
