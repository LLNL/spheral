# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import os

class Spheral(CMakePackage, PythonPackage):
    """FIXME: Put a proper description of your package here."""

    homepage = "https://spheral.readthedocs.io/"
    git      = "https://github.com/llnl/spheral.git"
    tags     = ['radiuss', 'simulations', 'hydrodynamics']

    maintainers = ['mdavis36','jmikeowen']

    version('develop', branch='feature/spack', submodules=True)
    version('1.0', tag='FSISPH-v1.0', submodules=True)

    variant('mpi', default=True, description='Enable MPI Support.')
    variant('openmp', default=True, description='Enable OpenMP Support.')
    variant('docs', default=False, description='Enable building Docs.')

    depends_on('mpi', type=['build','run'], when='+mpi')
    depends_on('cmake@3.10.0:', type='build')

    depends_on('zlib@1.2.11 -shared +pic', type='build')

    depends_on('boost@1.74.0 -atomic -container -coroutine -chrono -context -date_time -exception -fiber -graph -iostreams -locale -log -math -mpi -program_options -python -random -regex -serialization -test -thread -timer -wave +pic', type='build')

    # TODO: ANEOS
    # TODO: qhull spack package seems to be broken...
    #depends_on('qhull@2020.2', type='build')
    #TODO: Polyclipper package.

    depends_on('eigen@3.3.7', type='build')
    depends_on('hdf5@1.8.19 ~mpi +hl', type='build')
    depends_on('silo@4.10.2 +hdf5', type='build')

    # Zlib fix has been merged into conduit, using develop until next release.
    #depends_on('conduit@develop +mpi +hdf5 -shared -test', type='build')
    depends_on('conduit@develop +mpi +hdf5 -test', type=['build','run'], when='+mpi')
    depends_on('conduit@develop ~mpi +hdf5 -test', type=['build','run'], when='~mpi')

    depends_on('axom@0.5.0 +mpi +hdf5 -lua -examples -python -fortran -umpire -raja', type=['build','run'], when='+mpi')
    depends_on('axom@0.5.0 ~mpi +hdf5 -lua -examples -python -fortran -umpire -raja', type=['build','run'], when='~mpi')

    depends_on('opensubdiv@3.4.3', type='build')
    depends_on('polytope', type='build')

    extends('python@2.7.16 +zlib +shared', type='build')

    depends_on('py-pip@9.0.1', type='build')
    depends_on('py-pybind11@2.4.3', type='build')
    depends_on('py-pyb11generator@1.0.12', type='build')


    def cmake_args(self):
        options = []
        spec = self.spec

        options.append(self.define('ENABLE_CXXONLY', False))
        options.append(self.define('TPL_VERBOSE', False))
        options.append(self.define('BUILD_TPL', True))

        options.append(self.define('python_BUILD', False))
        options.append(self.define('python_DIR', spec['python'].prefix))

        options.append(self.define('zlib_BUILD', False))
        options.append(self.define('zlib_DIR', spec['zlib'].prefix))

        options.append(self.define('boost_BUILD', False))
        options.append(self.define('boost_DIR', spec['boost'].prefix))

        #options.append(self.define('qhull_BUILD', False))
        #options.append(self.define('qhull_DIR', spec['qhull'].prefix))

        options.append(self.define('hdf5_BUILD', False))
        options.append(self.define('hdf5_DIR', spec['hdf5'].prefix))

        options.append(self.define('conduit_BUILD', False))
        options.append(self.define('conduit_DIR', spec['conduit'].prefix))

        options.append(self.define('axom_BUILD', False))
        options.append(self.define('axom_DIR', spec['axom'].prefix))

        options.append(self.define('silo_BUILD', False))
        options.append(self.define('silo_DIR', spec['silo'].prefix))

        options.append(self.define('eigen_BUILD', False))
        options.append(self.define('eigen_DIR', spec['eigen'].prefix))
        options.append(self.define('eigen_INCLUDES', spec['eigen'].prefix.include.eigen3))

        options.append(self.define('opensubdiv_BUILD', False))
        options.append(self.define('opensubdiv_DIR', spec['opensubdiv'].prefix))

        options.append(self.define('pip_BUILD', False))
        options.append(self.define('pip_DIR', spec['py-pip'].prefix + '/lib/python2.7/site-packages/'))

        options.append(self.define('setuptools_BUILD', False))
        options.append(self.define('setuptools_DIR', spec['py-setuptools'].prefix + '/lib/python2.7/site-packages/'))

        options.append(self.define('pybind11_BUILD', False))
        options.append(self.define('pybind11_DIR', spec['py-pybind11'].prefix))

        options.append(self.define('pyb11generator_BUILD', False))
        options.append(self.define('pyb11generator_DIR', spec['py-pyb11generator'].prefix + '/lib/python2.7/site-packages/'))

        options.append(self.define('polytope_BUILD', False))
        options.append(self.define('polytope_DIR', spec['polytope'].prefix))

        options.append(self.define_from_variant('ENABLE_MPI', 'mpi'))
        if "+mpi" in spec:
            options.extend([
                "-DMPI_C_COMPILER=%s" % spec['mpi'].mpicc,
                "-DMPI_CXX_COMPILER=%s" % spec['mpi'].mpicxx
                ])

        options.append(self.define_from_variant('ENABLE_OPENMP', 'openmp'))
        options.append(self.define_from_variant('ENABLE_DOCS', 'docs'))

        return options
