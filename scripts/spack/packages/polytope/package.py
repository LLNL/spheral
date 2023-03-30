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

    # version('0.7.3', sha256='f32817b44d2a3b98407531980b89d0a31b0c14b8b30de37a6a7bc6ec91e48bf1') # missing submodules
    #version('0.7.2', tag='0.7.2', submodules=True)
    #version("2023-03-24", commit="074ec17", submodules=True)

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
        #options.append(self.define('PYB11GENERATOR_ROOT_DIR', 
        #options.append(self.define('PYBIND11_ROOT_DIR', spec['py-pybind11'].prefix.include))

        if "+python" in spec:
            options.append(self.define('USE_PYTHON', True))
            options.append(self.define('Python3_EXECUTABLE', os.path.join(self.spec['python'].prefix.bin, 'python3') ) )

        options.append(self.define('TESTING', False))

        return options
