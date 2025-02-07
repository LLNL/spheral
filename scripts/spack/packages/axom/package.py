# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import glob
import os
import shutil
import socket
from os.path import join as pjoin

from spack.package import *
from spack.pkg.builtin.axom import Axom as BuiltinAxom


class Axom(BuiltinAxom):
    """Axom provides a robust, flexible software infrastructure for the development
    of multi-physics applications and computational tools."""

    def initconfig_mpi_entries(self):
        spec = self.spec
        entries = BuiltinAxom.initconfig_mpi_entries(self)

        if "+mpi" in spec:
            # Replace /usr/bin/srun path with srun flux wrapper path on TOSS 4
            # TODO: Remove this logic by adding `using_flux` case in
            #  spack/lib/spack/spack/build_systems/cached_cmake.py:196 and remove hard-coded
            #  path to srun in same file.
            if "toss_4" in self._get_sys_type(spec):
                mpi_exec_index = [
                    index for index, entry in enumerate(entries) if "MPIEXEC_EXECUTABLE" in entry
                ]
                if mpi_exec_index:
                    del entries[mpi_exec_index[0]]

        return entries
