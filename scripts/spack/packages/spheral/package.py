# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
import spack
import socket
import os

class Spheral(CachedCMakePackage, CudaPackage, ROCmPackage):
    """Spheral++ provides a steerable parallel environment for performing coupled hydrodynamical and gravitational numerical simulations."""

    homepage = "https://spheral.readthedocs.io/"
    git      = "https://github.com/llnl/spheral.git"
    tags     = ['radiuss', 'simulations', 'hydrodynamics']

    maintainers = ['mdavis36','jmikeowen','owen32']

    # -------------------------------------------------------------------------
    # VERSIONS
    # -------------------------------------------------------------------------
    version('develop', branch='develop', submodules=True)
    version('1.0', tag='FSISPH-v1.0', submodules=True)

    # -------------------------------------------------------------------------
    # VARIANTS
    # -------------------------------------------------------------------------
    variant('mpi', default=True, description='Enable MPI Support.')
    variant('openmp', default=True, description='Enable OpenMP Support.')
    variant('docs', default=False, description='Enable building Docs.')
    variant('shared', default=True, description='Build C++ libs as shared.')
    variant('python', default=True, description='Build Python Dependencies.')

    # -------------------------------------------------------------------------
    # DEPENDS
    # -------------------------------------------------------------------------
    depends_on('mpi', when='+mpi')

    depends_on('cmake@3.21.0:', type='build')

    depends_on('boost@1.74.0 +system +filesystem -atomic -container -coroutine -chrono -context -date_time -exception -fiber -graph -iostreams -locale -log -math -mpi -program_options -python -random -regex -test -thread -timer -wave +pic', type='build')

    depends_on('zlib@1.3 +shared +pic', type='build')

    depends_on('qhull@2020.2 +pic', type='build')

    depends_on('m-aneos@1.0')

    depends_on('eigen@3.4.0', type='build')

    depends_on('hdf5 +hl', type='build')

    depends_on('silo +hdf5', type='build')

    depends_on('conduit@0.9.1 +shared +hdf5~hdf5_compat -test ~parmetis', type='build')

    depends_on('axom@0.9.0 +hdf5 -lua -examples -python -fortran', type='build')
    depends_on('axom +shared', when='~cuda', type='build')
    depends_on('axom ~shared', when='+cuda', type='build')

    depends_on('caliper@2.11 ~shared +adiak +gotcha ~libdw ~papi ~libunwind +pic', type='build')

    depends_on("raja@2024.02.0", type="build")

    depends_on('opensubdiv@3.4.3', type='build')

    depends_on('polytope +python', type='build', when='+python')

    # Forward MPI Variants
    mpi_tpl_list = ["hdf5", "conduit", "axom", "caliper", "adiak~shared"]
    for ctpl in mpi_tpl_list:
        for mpiv in ["+mpi", "~mpi"]:
            depends_on(f"{ctpl} {mpiv}", type='build', when=f"{mpiv}")

    # Forward CUDA/ROCM Variants
    gpu_tpl_list = ["raja", "umpire", "axom"]
    for ctpl in gpu_tpl_list:
        for val in CudaPackage.cuda_arch_values:
            depends_on(f"{ctpl} +cuda cuda_arch={val}", type='build', when=f"+cuda cuda_arch={val}")
        for val in ROCmPackage.amdgpu_targets:
            depends_on(f"{ctpl} +rocm amdgpu_target={val}", type='build', when=f"+rocm amdgpu_target={val}")

    # -------------------------------------------------------------------------
    # Conflicts
    # -------------------------------------------------------------------------
    conflicts("+cuda", when="+rocm")

    def _get_sys_type(self, spec):
        sys_type = spec.architecture
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        return sys_type

    # Create a name for the specific configuration being built
    # This name is used to differentiate timings during performance testing
    def _get_config_name(self, spec):
        sys_type = self._get_sys_type(spec)
        config_name = f"{sys_type}_{spec.compiler.name}_{spec.compiler.version}"
        if ("+mpi" in spec):
            config_name += "_" + spec.format("{^mpi.name}_{^mpi.version}")
        if ("+cuda" in spec):
            config_name += "_" + spec.format("{^cuda.name}{^cuda.version}")
        return config_name.replace(" ", "_")

    @property
    def cache_name(self):

        hostname = socket.gethostname()
        if "SYS_TYPE" in env:
            hostname = hostname.rstrip('1234567890')

        envspec = os.environ.get("SPEC")
        spec = self.spec
        if envspec:
          cache_spec = envspec
        else:
          cache_spec = str(spec.compiler.name) + "@" + str(spec.compiler.version)
        return "{0}-{1}.cmake".format(
            str(self._get_sys_type(spec)),
            cache_spec.replace(" ", "_")
        )

    def initconfig_compiler_entries(self):
        spec = self.spec
        entries = super(Spheral, self).initconfig_compiler_entries()
        return entries
    
    def initconfig_mpi_entries(self):
        spec = self.spec
        entries = []
        if "+mpi" in spec:
          entries = super(Spheral, self).initconfig_mpi_entries()
          # When on cray / flux systems we need to tell CMAKE the mpi flag explicitly
          if "cray-mpich" in spec:
            for e in entries:
                if 'MPIEXEC_NUMPROC_FLAG' in e:
                    entries.remove(e)
            entries.append(cmake_cache_string('MPIEXEC_NUMPROC_FLAG', '-n'))
        return entries

    def initconfig_hardware_entries(self):
        spec = self.spec
        entries = super(Spheral, self).initconfig_hardware_entries()

        if '+rocm' in spec:
            entries.append(cmake_cache_option("ENABLE_HIP", True))
            entries.append(cmake_cache_string("ROCM_PATH", spec["hip"].prefix))

        if '+cuda' in spec:
            entries.append(cmake_cache_option("ENABLE_CUDA", True))

            if not spec.satisfies('cuda_arch=none'):
                cuda_arch = spec.variants['cuda_arch'].value
                entries.append(cmake_cache_string(
                    "CUDA_ARCH", 'sm_{0}'.format(cuda_arch[0])))
                entries.append(cmake_cache_string(
                    "CMAKE_CUDA_ARCHITECTURES", '{0}'.format(cuda_arch[0])))
                flag = '-arch sm_{0}'.format(cuda_arch[0])
                entries.append(cmake_cache_string(
                    "CMAKE_CUDA_FLAGS", '{0}'.format(flag)))

            entries.append(cmake_cache_option(
                "ENABLE_DEVICE_CONST", spec.satisfies('+deviceconst')))
        else:
            entries.append(cmake_cache_option("ENABLE_CUDA", False))

        return entries

    def initconfig_package_entries(self):
        spec = self.spec
        entries = []

        entries.append(cmake_cache_option('ENABLE_CXXONLY', False))
        entries.append(cmake_cache_option('TPL_VERBOSE', False))
        entries.append(cmake_cache_option('BUILD_TPL', True))

        entries.append(cmake_cache_string('SPHERAL_SYS_ARCH', self._get_sys_type(spec)))
        entries.append(cmake_cache_string('SPHERAL_CONFIGURATION', self._get_config_name(spec)))

        # TPL locations
        entries.append(cmake_cache_path('caliper_DIR', spec['caliper'].prefix))

        entries.append(cmake_cache_path('adiak_DIR', spec['adiak'].prefix))

        entries.append(cmake_cache_path('boost_DIR', spec['boost'].prefix))

        entries.append(cmake_cache_path('qhull_DIR', spec['qhull'].prefix))

        entries.append(cmake_cache_path('aneos_DIR', spec['m-aneos'].prefix))

        entries.append(cmake_cache_path('hdf5_DIR', spec['hdf5'].prefix))
    
        entries.append(cmake_cache_path('conduit_DIR', spec['conduit'].prefix))

        entries.append(cmake_cache_path('raja_DIR', spec['raja'].prefix))

        entries.append(cmake_cache_path('umpire_DIR', spec['umpire'].prefix))

        entries.append(cmake_cache_path('axom_DIR', spec['axom'].prefix))

        entries.append(cmake_cache_path('silo_DIR', spec['silo'].prefix))

        entries.append(cmake_cache_path('eigen_DIR', spec['eigen'].prefix))
        entries.append(cmake_cache_path('eigen_INCLUDES',spec['eigen'].prefix.include.eigen3))

        entries.append(cmake_cache_path('opensubdiv_DIR', spec['opensubdiv'].prefix))

        entries.append(cmake_cache_path('polytope_DIR', spec['polytope'].prefix))

        entries.append(cmake_cache_option('ENABLE_MPI', '+mpi' in spec))
        if "+mpi" in spec:
            entries.append(cmake_cache_path('-DMPI_C_COMPILER', spec['mpi'].mpicc) )
            entries.append(cmake_cache_path('-DMPI_CXX_COMPILER', spec['mpi'].mpicxx) )

        if "~shared" in spec and "~cuda" in spec:
            entries.append(cmake_cache_option('ENABLE_SHARED', False))

        entries.append(cmake_cache_option('ENABLE_OPENMP', '+openmp' in spec))
        entries.append(cmake_cache_option('ENABLE_DOCS', '+docs' in spec))

        if "+python" in spec:
            entries.append(cmake_cache_path('python_DIR', spec['python'].prefix))

        return entries


    def cmake_args(self):
        options = []
        spec = self.spec

        return options

    @property
    def build_dirname(self):
        """Directory name to use when building the package."""
        return "spack-build-%s" % self.pkg.spec.dag_hash(7)

    @property
    def build_directory(self):
        """Full-path to the directory to use when building the package."""
        spec = self.spec
        if spec.satisfies("@develop"):
            dev_build_dir = "spack-build-" + str(spec.compiler.name) + "-" + str(spec.compiler.version)
            return os.path.join(self.pkg.stage.source_path, build_dirname)
        else:
            return os.path.join(self.pkg.stage.path, self.build_dirname)

