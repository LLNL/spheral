# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import os


class Polytope(CMakePackage):
    """Polytope is a C++ library for generating polygonal and polyhedral meshes."""

    #url = "https://github.com/pbtoast/polytope/archive/0.7.1.tar.gz"
    #version('0.7.0', sha256='e7be5a5d06a95309b9100f02ab80c2f3401c11ac6304a3204fb1f25052efd77e')

    git = "https://github.com/pbtoast/polytope.git"
    #version('0.7.2', tag='0.7.2', submodules=True)
    version("2023-03-23", commit="77c0f29", submodules=True)

    variant('python', default=True, description='Enable Python Support.')

    extends('python', when='+python')
    depends_on('python@3: +zlib +shared', type=('build', 'run'), when='+python')
    depends_on('py-decorator', type=('build', 'run'), when='+python')
    depends_on('boost', type=('build', 'run'))

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
