# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyPolyclipper(CMakePackage, PythonPackage):
    """Polyclipper"""

    homepage = "https://pypi.org/project/PYB11Generator/"
    url      = "https://github.com/LLNL/PolyClipper/archive/refs/tags/v1.2.3.zip"
    git      = "https://github.com/LLNL/PolyClipper"

    maintainers = ['mdavis36','jmikeowen']

    version('1.2.3', sha256='366e547bc343033c760727b6cdbf34a304c27bc769a208e9bfaeec42c92dba96')

    variant('mpi', default=False, description='Enable MPI Support.')
    variant('openmp', default=True, description='Enable OpenMP Support.')
    variant('docs', default=False, description='Enable building Docs.')

    depends_on('mpi', when='+mpi')
    depends_on('blt')
    depends_on('py-pybind11')
    depends_on('py-pyb11generator')
    depends_on('py-decorator')

    def cmake_args(self):
        spec = self.spec
        args = []

        args.append(self.define('POLYCLIPPER_BLT_DIR', spec['blt'].prefix))
        args.append(self.define('ENABLE_CXXONLY', True))
        args.append(self.define('PYTHON_EXE', spec['python'].prefix+'/bin/python'))
        args.append(self.define('PYBIND11_INCLUDE_PATH', spec['py-pybind11'].prefix+'/include'))
        args.append(self.define('PYB11GEN_PATH', spec['py-pyb11generator'].prefix+'/lib/python2.7/site-packages'))

        args.append(self.define('ENABLE_MPI', '+mpi' in spec))
        if "+mpi" in spec:
            args.append(self.define('MPI_C_COMPILER', spec['mpi'].mpicc) )
            args.append(self.define('MPI_CXX_COMPILER', spec['mpi'].mpicxx) )

        args.append(self.define('ENABLE_OPENMP', '+openmp' in spec))
        args.append(self.define('ENABLE_DOCS', '+docs' in spec))

        return args
