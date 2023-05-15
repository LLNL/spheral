# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import os


class Polytope(CMakePackage):
    """Polytope is a C++ library for generating polygonal and polyhedral meshes."""

    git = "https://github.com/pbtoast/polytope.git"
    url = "https://github.com/pbtoast/polytope/archive/0.7.3.tar.gz"
    version('0.7.3', tag='0.7.3', submodules=True)

    variant('python', default=True, description='Enable Python Support.')

    extends('python', when='+python')
    depends_on('python@3: +zlib +shared', type=('build', 'run'), when='+python')
    depends_on('py-decorator', type=('build', 'run'), when='+python')
    depends_on('boost', type=('build', 'run'))

    parallel = False      # Should be able to remove this at some point

    def cmake_args(self):
        options = []
        spec = self.spec

        options.append(self.define('USE_MPI', 'Off'))   # Turn back on when polytope fixes parallel generation
        options.append(self.define('BOOST_ROOT', spec['boost'].prefix.include))

        if "+python" in spec:
            options.append(self.define('USE_PYTHON', True))
            options.append(self.define('Python3_EXECUTABLE', os.path.join(self.spec['python'].prefix.bin, 'python3') ) )

        options.append(self.define('TESTING', False))

        return options
