# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import sys

from spack.package import *
from spack.pkg.builtin.hypre import Hypre as BuiltinHypre

class Hypre(BuiltinHypre):
    """Hypre is a library of high performance preconditioners that
    features parallel multigrid methods for both structured and
    unstructured grid problems."""
    variant("pic", default=False, description="Build with -fPIC")

    def setup_build_environment(self, env):
        spec = self.spec
        if "+pic" in spec:
            env.set("CXXFLAGS", "-fPIC")
            env.set("CFLAGS", "-fPIC")
