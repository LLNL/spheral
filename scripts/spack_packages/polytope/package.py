# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install spheral
#
# You can edit this file again by typing:
#
#     spack edit spheral
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *
import os


class Polytope(CMakePackage):
    """FIXME: Put a proper description of your package here."""

    url = "https://github.com/pbtoast/polytope/archive/0.6.2.tar.gz"
    git = "https://github.com/pbtoast/polytope.git"

    version('0.6.2', sha256='e9ed18c3ebc7b4b231a0563235cc032c26daa7de88839c64141975170e55bcfd')

    patch('polytope-PYB11-CMakeLists.patch')

    variant('python', default=True, description='Enable Python Support.')

    depends_on('python@2.7:2.8 +zlib +shared', type=['build', 'run'], when='+python')
    depends_on('boost')
    depends_on('py-pybind11@2.4.3')
    depends_on('py-pyb11generator')

    def cmake_args(self):
        options = []
        spec = self.spec

        options.append(self.define('PYBIND11_INCLUDE_DIRS', spec['py-pybind11'].prefix.include))
        options.append(self.define('Boost_INCLUDE_DIR', spec['boost'].prefix.include))

        if "+python" in spec:
            options.append(self.define('USE_PYTHON', True))
            options.append(self.define('PYTHON_EXE', os.path.join(self.spec['python'].prefix.bin, 'python') ) )

        options.append(self.define('TESTING', False))

        return options
