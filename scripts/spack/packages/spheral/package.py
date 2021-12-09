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

    #depends_on('libevent', type=['build','run'])
    depends_on('mpi', type=['build','run'], when='+mpi')
    depends_on('cmake@3.10.0:', type='build')

    depends_on('zlib@1.2.11 -shared', type='build')

    # TODO make a cleaner way of listing which libs to build here...
    depends_on('boost@1.74.0 -atomic -container -coroutine -chrono -context -date_time -exception -fiber -graph -iostreams -locale -log -math -mpi -program_options -python -random -regex -serialization -test -thread -timer -wave', type='build')

    # TODO: qhull spack package seems to be broken...
    #depends_on('qhull@2020.2', type='build')

    # Eigen will "build" now but oh well
    depends_on('eigen@3.3.7', type='build')

    # hdf5 will need +hl to build the higher_level lib we expect.
    depends_on('hdf5@1.8.19 ~mpi +hl', type='build')
    #depends_on('hdf5@1.10.4 +hl', type='build')

    # Need to depend on hdf5 @44
    depends_on('silo@4.10.2 +hdf5', type='build')

    # Zlib fix has been merged into conduit, using develop until next release.
    #depends_on('conduit@develop +mpi +hdf5 -shared -test', type='build')
    depends_on('conduit@develop +mpi +hdf5 -test', type=['build','run'], when='+mpi')
    depends_on('conduit@develop ~mpi +hdf5 -test', type=['build','run'], when='~mpi')

    # TODO: axom machine broke
    depends_on('axom@0.5.0 +mpi +hdf5 -lua -examples -python -fortran -umpire -raja', type=['build','run'], when='+mpi')
    depends_on('axom@0.5.0 ~mpi +hdf5 -lua -examples -python -fortran -umpire -raja', type=['build','run'], when='~mpi')

    # TODO: ANEOS

    # Upgrading opensubdiv version 3.3.0 is not available via spack
    # TODO: Look into editing or making custom opensubdiv spack package.
    # OpenSubDiv plus all dependencies took over 2 hours to build.
    depends_on('opensubdiv@3.4.3', type='build')

    extends('python@2.7.16 +zlib +shared', type='build')

    depends_on('py-pip@9.0.1', type='build')
    depends_on('py-pybind11@2.4.3', type='build')
    depends_on('py-pyb11generator@1.0.12', type='build')

    #depends_on('py-numpy@1.16.6')
    #depends_on('py-mpi4py')
    #depends_on('py-numpy-stl@2.11.2') # MikeO: What do we use this for? do we need it?
    #depends_on('py-matplotlib@2.2.5') 
    ##depends_on('py-pillow@6.2.0') # Need this to saitsify dep for matplotlib.
    #depends_on('py-gnuplot@1.8')
    #depends_on('py-enum34')
    #depends_on('py-python-utils')

    #depends_on('py-setuptools')
    #depends_on('py-decorator@4.4.2')
    #depends_on('py-h5py@2.10.0')
    #depends_on('py-sphinx@1.8.5')
    #depends_on('py-sphinx-rtd-theme@0.5.0')
    ##depends_on('py-twine@2.0.0', type='build') # TODO: Only v2.0.0 available (Python3), usually request v1.15.0
    #depends_on('py-cython@0.29.21')
    ##depends_on('sobol@0.9', type='build') # MikeO: Do we need this, there is sobol-seq v0.2.0
    #depends_on('py-scipy@1.2.3')
    ###pipreqs==0.4.10 # TODO: DNE
    #depends_on('py-importlib-metadata@2.0.0') # WARNING: v2.1.1 dne. 2.0.0 available
    #depends_on('py-virtualenv') # WARNING: 20.2.2 dne, 16.7.6 is the closest in pkg

    #TODO: ATS


    #TODO: Polytope spack package.
    depends_on('polytope', type='build')
    #TODO: Polyclipper package.

    #def setup_run_environment(self, env):
    #    env('SPACK_PYTHON', os.path.join(self.spec['python'].prefix.bin, 'python') )

    def cmake_args(self):
        options = []
        spec = self.spec

        options.append(self.define('ENABLE_CXXONLY', False))
        options.append(self.define('TPL_VERBOSE', False))
        options.append(self.define('BUILD_TPL', True))
        #options.append(self.define('SPACK_INSTALL', True))

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

        options.append(self.define('pybind11_BUILD', False))
        options.append(self.define('pybind11_DIR', spec['py-pybind11'].prefix))

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
