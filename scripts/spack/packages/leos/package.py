from spack.package import *
import os, re
import llnl.util.tty as tty

class Leos(CMakePackage, CudaPackage, ROCmPackage):
    """This is derived from llnl.wci.legacy Leos package, but makes silo optional
       and excludes ascmemory.
    """
    fileLoc = 'file:///usr/gapps/leos/srcs'
    version("8.4.1", sha256="93abeea9e336e3a81cc6cc9de10b2a2fd61eda2a89abece50cac80fef58ec38b",
            url=os.path.join(fileLoc, "leos-8.4.1.tar.gz"))
    version("8.4.2", sha256="08eb87580e30d7a1db72b1e1a457652dda9535df1c0caf7b5badb9cadf39f2a9",
            url=os.path.join(fileLoc, "leos-8.4.2.tar.gz"))

    variant('debug',   default=False, description='Build debug code (-g -O0)')
    variant('filters', default=True , description='Build LEOS filter coding')
    variant('raja',    default=False, description='Build LIP using RAJA')
    variant('cuda',    default=False, description='Build LIP using RAJA + CUDA')
    variant('umpire',  default=False, description='Build LIP using UMPIRE')
    variant('silo',    default=True,  description='Use Silo instead of LEOSPACT')
    variant("cpp14",   default=True, description="Sets the C++ std to C++14")
    variant("yaml",   default=True, description="Enable use of YAML")

    depends_on("mpi")
    depends_on("hdf5+hl")
    depends_on("silo", when="+silo")
    depends_on("boost")
    depends_on("cmake")
    depends_on("zlib")
    depends_on("raja", when="+raja")
    depends_on("raja+cuda", when="+cuda")
    depends_on("cnmem", when="+cuda @:8.2.5")
    depends_on("umpire+cuda", when="+cuda")
    depends_on("umpire", when="+umpire")
    depends_on("camp", when="+umpire")
    depends_on("camp", when="+raja")
    depends_on("blt", type="build")
    depends_on("camp")

    patch('patches/leos-8.2.0-CRAY.patch',  when='@8.0.5:8.2.0 %cce')
    patch('patches/leos-8.4.1-umpire-import.patch',  when='@8.4.1')
    patch('patches/leos-8.4.1-umpire-2024.02-exceptions.patch', when='@8.4.1 ^umpire@2024.02.0:')

    fpic_compilers = ['intel@14.toss',
                      'intel@16.toss',
                      'intel@17.toss',
                      'intel@18.toss',
                      'intel@19.toss',
                      'intel@toss',
                      ]

    c_o2_compilers = ['pgi',
                      ]

    @property
    def libs(self):
        libnames = ['leos', 'leos_F', 'leos_C', 'lip-cpp']
        if (self.spec.satisfies('+yaml')):
            libnames.append('yaml-cpp')
        if (self.spec.satisfies('~silo')):
            libnames.append('leospact')
        libs = ['lib' + x for x in libnames]
        return find_libraries(libs, self.spec.prefix, shared=False, recursive=True)

    def setup_build_environment(self, env):
        spec = self.spec
        if '+cuda' in spec:
            env.set('CUDAHOSTCXX', spack_cxx)

    def cmake_args(self):
        spec = self.spec
        args = []

        if spec.satisfies('+silo'):
            args.extend([
                '-DSILO=ON',
                "-DSILO_PATH=%s" % spec['silo'].prefix,
                "-DUSE_PDB=ON",
            ])
        else:
            args.append('-DSILO=OFF')

        if '+cuda' in spec:
            cuda_flags = None
            if '%xl' in spec:
                cuda_flags = '--expt-extended-lambda -Xptxas -warn-lmem-usage,-warn-spills -DCAMP_USE_PLATFORM_DEFAULT_STREAM'
            elif '%gcc' in spec:
                # Note: %gcc+cuda hasn't been tried yet
                cuda_flags = '--expt-extended-lambda'
            args.extend([
                '-DENABLE_CUDA:BOOL=ON',
                '-DCMAKE_CUDA_COMPILER=%s/nvcc' % spec['cuda'].prefix.bin,
                '-DCMAKE_CUDA_STANDARD=11',
                '-DCMAKE_CUDA_SEPARABLE_COMPILATION=ON',
                # Note this CMake option name for setting CUDA is nonstandard
                '-DCUDA_ARCH=sm_%s' % spec.variants['cuda_arch'].value,
                '-DCMAKE_CUDA_ARCHITECTURES=%s' % spec.variants['cuda_arch'].value,
            ])
            if cuda_flags:
                args.append('-DCMAKE_CUDA_FLAGS={0}'.format(cuda_flags))

        if '+rocm' in spec:
            args.extend([
                "-DENABLE_HIP=ON",
                '-DHIP_CXX_COMPILER=' + str(self.compiler.cxx),
            ])

        cxx_flags = []
        if spec.satisfies('+raja'):
            args.extend([
                "-DENABLE_RAJA=ON",
                "-DRAJA_DIR=%s" % spec['raja'].prefix,
                "-Dcamp_DIR={}/lib/cmake/camp".format(spec["camp"].prefix)
            ])
            if spec.satisfies('+cuda'):
                args.append('-DENABLE_RAJA_CUDA=ON')
                cxx_flags.append(spec['cuda'].libs.ld_flags)

        if '+rocm' in spec:
            # NOTE: This should not be required: FindRaja.cmake and raja-config.cmake
            # should arrange for camp include directory. However, for some reason
            # these files, in conjunction with cmake, are not performing the
            # transitive search
            hipcc_flags = ['-std=c++14']
            hipcc_flags.append('-I{0}'.format(spec['camp'].prefix.include))
            if 'camp+rocm' in spec:
                hipcc_flags.append('-DCAMP_HAVE_HIP')
            args.append('-DHIP_HIPCC_FLAGS={0}'.format(' '.join(hipcc_flags))) 

        if cxx_flags:
            args.append('-DCMAKE_CXX_FLAGS=' + ' '.join(cxx_flags))

        if spec.satisfies('+umpire'):
            umpire_libs = ';'.join(spec['umpire'].libs.libraries)
            args.extend([
                "-DENABLE_UMPIRE=ON",
                "-DENABLE_UMPIRE_IDS=ON",
                "-DUMPIRE_DIR=%s" % spec['umpire'].prefix,
                '-DUMPIRE_LIBRARIES=%s' % umpire_libs,
                "-Dcamp_DIR={}/lib/cmake/camp".format(spec["camp"].prefix)
            ])
        else:
            args.extend([
                '-DENABLE_UMPIRE=OFF',
                '-DENABLE_UMPIRE_IDS=OFF'
            ])

        if spec.satisfies('+filters'):
            args.extend([
                "-DENABLE_FILTERS=ON",
                "-DUSE_LIBXML2=OFF",
                "-DUSE_NOXML=ON",
                "-DFILTER_PLUGIN_DIR={0}".format(spec.prefix.lib)
            ])
        else:
            args.extend([
                "-DENABLE_FILTERS=OFF",
                "-DUSE_LIBXML2=OFF",
                "-DUSE_NOXML=ON"
            ])
        if spec.satisfies('+yaml'):
            args.append("-DENABLE_YAML:BOOL=ON")
        else:
            args.append("-DENABLE_YAML:BOOL=OFF")

        args += [
            '-DBLT_SOURCE_DIR={0}'.format(spec['blt'].prefix),

            '-DLEOS_BUILD_STANDALONE=OFF',

            '-DENABLE_MPI=ON',
            '-DUSE_MPI=ON',
            '-DMPI_C_COMPILER=%s' % self.spec['mpi'].mpicc,
            '-DMPI_CXX_COMPILER=%s' % self.spec['mpi'].mpicxx,
            '-DMPI_Fortran_COMPILER=%s' % self.spec['mpi'].mpifc,

            '-DENABLE_PYTHON=OFF',
            '-DENABLE_ZFP=OFF',
            '-DENABLE_OPENMP=OFF',
            '-DENABLE_CALIPER=OFF',

            '-DENABLE_TESTS=OFF',
            '-DBUILD_TESTS=OFF',
            "-DENABLE_GMOCK=OFF",
            "-DENABLE_GTEST=OFF",
            '-DBUILD_EXAMPLES=OFF',
            '-DBUILD_EXAMPLES_P=OFF',
            '-DENABLE_EXAMPLES=OFF',
            '-DENABLE_PARALLEL_EXAMPLES=OFF',

            "-DBUILD_TOOLS=OFF",
            "-DENABLE_TOOLS=OFF",

            "-DUSE_TCALCFILTER=ON",
            "-DENABLE_V7=ON",
            "-DSTATIC_FILTERS=ON",
            "-DENABLE_COMPRESSION=OFF",
            "-DENABLE_SM_LIP=OFF",
            "-DENABLE_C_LIP=ON",
            "-DENABLE_CPP_LIP=ON",
            "-DLIP_SHARED_MEMORY=OFF",
            "-DENABLE_ASCMEMORY=OFF",
            "-DENABLE_FPE_HANDLING=OFF",
            "-DENABLE_UNCRUSTIFY=OFF",
            "-DENABLE_VALGRIND=OFF",

            "-DUSE_HDF=ON",
            "-DHDF_ROOT=%s" % spec['hdf5'].prefix,
            "-DHDF5_ROOT=%s" % spec['hdf5'].prefix,

            "-DENABLE_FORTRAN=ON",
            "-DENABLE_FORTRAN_INTERFACE=ON",
            "-DENABLE_C=ON",
            "-DENABLE_C_INTERFACE=ON",

            '-DEOS_DATA_ROOT_DIR=/usr/gapps/data/eos',
        ]

        if "+cpp14" in spec:
            args.extend(["-DBLT_CXX_STD=c++14", "-DCMAKE_CXX_STANDARD=14"])

        # If present, pass in language-specific flags to CMake
        cppflags = " ".join(spec.compiler_flags["cppflags"])
        if cppflags:
            # avoid always ending up with ' ' with no flags defined
            cppflags += " "
        cflags = cppflags + " ".join(spec.compiler_flags["cflags"])
        if cflags:
            args.append(f'-DBLT_C_FLAGS="{cflags}"')

        cxxflags = cppflags + " ".join(spec.compiler_flags["cxxflags"])
        if cxxflags:
            args.append(f'-DBLT_CXX_FLAGS="{cxxflags}"')

        fflags = " ".join(spec.compiler_flags["fflags"])
        if spec.satisfies("%cce"):
            fflags += " -ef"
        if spec.satisfies("%rocmcc"):
            fflags += " -O3 -Mstandard -Mfreeform -Mextend -Mbackslash -fno-default-real-8 -cpp"
        if fflags:
            args.append(f'-DBLT_FORTRAN_FLAGS={fflags}')

        return args

    def build(self, spec, prefix):
        fflags_path = join_path(self.build_directory, 'interfaces/fortran/CMakeFiles/leos_F.dir/flags.make')
        tty.msg("Fortran flags path: " + fflags_path)
        filter_file('-DAND', '', fflags_path)
        super(Leos, self).install(spec, prefix)

    def patch(self):
        if self.spec.satisfies('@8.4.1:  ^raja@2024.02.2:'):
            filter_file('loop_', 'seq_', "LIP/cpp_source/raja_api/SimpleRajaDefinitions.hxx", backup=True)

        if self.spec.satisfies('@8.2.0') and self.spec.satisfies('%gcc@10:'):
            filter_file(r'functional>', 'functional>\n#include<string>', 'src/implementation/LEOS_MemoryManager.hxx')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/fortran/LEOS_Function_F.f90')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.h')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.cc')
        elif self.spec.satisfies('@8.1.1p1') and ("%gcc" in self.spec or ("%clang" in self.spec and not self.spec.target.family=='ppc64le')): # gfortran issue
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/fortran/LEOS_Function_F.f90')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.h')
            filter_file(r'LEOS_Function_getMetadata_sS_len', 'LEOS_Function_getMetadata_sZ_len', 'interfaces/C/LEOS_Function_C.cc')
        if self.spec.satisfies('@8.0.5:'):
            filter_file(r'USE_MPI','USE_MPI AND BUILD_TESTS', 'interfaces/C/CMakeLists.txt')
            filter_file(r'USE_MPI','USE_MPI AND BUILD_TESTS', 'interfaces/fortran/CMakeLists.txt')
            filter_file(r'include\(leos-mpi\)','add_definitions(-DLEOS_USE_MPI)', 'CMakeLists.txt')
            if self.spec.satisfies('@8.2:'):
                filter_file(r'NOT \${HDF5_C_HL_LIBRARIES}', 'NOT "${HDF5_C_HL_LIBRARIES}"', 'cmake/Setup3rdParty.cmake')

            if "%pgi" in self.spec:
                filter_file(r' tolower', ' ::tolower', 'src/io/ReaderFactory.cpp')
                filter_file(r' -cpp','', 'interfaces/fortran/CMakeLists.txt')
                if self.spec.satisfies('@8.2:'):
                    filter_file(r'__x86_64__', '__FOO__', 'cereal/include/cereal/external/rapidjson/internal/biginteger.h')
                    filter_file(r'__x86_64__', '__FOO__', 'cereal/include/cereal/external/rapidjson/internal/diyfp.h') # avoid uint128 undefined

            if "%clang@bgq" in self.spec or "%gcc@bgq" in self.spec:
                filter_file(r'-WF,-D', '-D', 'interfaces/fortran/CMakeLists.txt')
        elif self.spec.satisfies('@8.0.0:'):
            if (self.spec.satisfies('platform=linux') and self.spec.target.family=='ppc64le') or ("%arm@sandia" in self.spec and self.spec.satisfies('@8.0.2:')):
                filter_file(r'typedef int LEOS_MPI_Comm', '#include "mpi.h"\n typedef MPI_Comm LEOS_MPI_Comm', 'interfaces/C/LEOS_types_C.h')

            if self.spec.satisfies("@8.0.4:"):
                datafile="src/io/ReaderFactory.cpp"
            else:
                datafile="src/io/databaseNames.h"
            filter_file(r'/usr/gapps/data/eos', self.datadir, datafile)
        if self.spec.satisfies('@8.0.5:8.1.0') and "%cce" in self.spec:
            filter_file(r'\${TMP_VERS}','"#define BOOST_LIB_VERSION \"1_54\""', 'CMAKE/leos-boost.cmake')
        if self.spec.satisfies('@8.0:8.2') and "%pgi" in self.spec:
            filter_file(r'LIP_setup.h \*/', 'LIP_setup.h */\n#ifndef LIPSETUPH\n#define LIPSETUPH', 'LIP/source/LIP_setup.h')
            with open("LIP/source/LIP_setup.h", 'a') as f:
                f.write("#endif")
            filter_file(r'LIP_proto.h \*/', 'LIP_proto.h */\n#ifndef LIPPROTOH\n#define LIPPROTOH', 'LIP/source/LIP_proto.h')
            filter_file(r'\}\n', '}\n#endif\n', 'LIP/source/LIP_proto.h')
            filter_file(r'LIP_setup.h', 'LIP_setup.h"\n#include "LIP_PBH.h', 'LIP/source/lip_setup_interp_shmem.cc')

    def old_flag_handler(self, name, flags):
        my_flags = flags
        if name in ('cflags', 'cxxflags'):
            if any(self.spec.compiler.satisfies(c) for c in self.fpic_compilers):
                my_flags.append(self.compiler.cc_pic_flag)
                my_flags.append('-g -fp-model precise -fp-model source -unroll-aggressive -finline-functions -nolib-inline -axCORE-AVX2 -xAVX -no-fma -ip -no-ftz -prec-div -prec-sqrt -diag-disable cpu-dispatch')
            else:
                my_flags.append('-g')
        if name == 'cxxflags' and any(self.spec.compiler.satisfies(c) for c in self.c_o2_compilers):
            my_flags = list(map(lambda x: "" if x == "-O2" else x, my_flags))
        return (None, my_flags, None)

    lanl_compilers = ['intel@lanl',
                      'intel@trinity',
                      'cce@trinity']

    sandia_compilers = ['intel@sandia',
                        'arm@sandia']

    @property
    def datadir(self):
        datadir = '/usr/gapps/data/eos'
        if any(self.spec.compiler.satisfies(c) for c in self.lanl_compilers):
            datadir = '/usr/projects/llnl_b/data/eos'
        if any(self.spec.compiler.satisfies(c) for c in self.sandia_compilers):
            datadir = '/projects/llnl_b/data/eos'
        return datadir

    def old_install(self, spec, prefix):
        compiler_install_path = "unset";

        std_cmake_args.append("-DCMAKE_C_FLAGS=%s" % env['CFLAGS'])
        std_cmake_args.append("-DCMAKE_CXX_FLAGS=%s" % env['CXXFLAGS'])

        if "%pgi@19.toss" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-19-mvapich2-2.3"

        if "%pgi@20.toss" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-20-mvapich2-2.3_NEW"

        if "%pgi@21.toss" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-21-mvapich2-2.3"

        if "%pgi@toss" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-21-mvapich2-2.3"

        if "%intel@14.toss" in spec:
            compiler_install_path = "intel-14-mvapich2-2.2"

        if "%intel@16.toss" in spec:
            compiler_install_path = "intel-16-mvapich2-2.2"

        if "%intel@17.toss" in spec:
            compiler_install_path = "intel-17-mvapich2-2.2"

        if "%intel@18.toss" in spec:
            compiler_install_path = "intel-18-mvapich2-2.3"

        if "%intel@19.toss" in spec:
            compiler_install_path = "intel-19-mvapich2-2.3"

        if "%intel@toss" in spec:
            compiler_install_path = "intel-19-mvapich2-2.3"

        if "%gcc@4.9.3.toss" in spec:
            compiler_install_path = "gnu-4.9-mvapich2-2.3"

        if "%gcc@7.toss" in spec:
            compiler_install_path = "gnu-7-mvapich2-2.3"

        if "%gcc@8.1.toss" in spec:
            compiler_install_path = "gnu-8-mvapich2-2.3"

        if "%gcc@8.3.toss" in spec:
            compiler_install_path = "gnu-8.3-mvapich2-2.3"

        if "%gcc@8.toss" in spec:
            compiler_install_path = "gnu-8.3-mvapich2-2.3"

        if "%gcc@toss" in spec:
            compiler_install_path = "gnu-8.3-mvapich2-2.3"

        if "%gcc@10.toss" in spec:
            compiler_install_path = "gnu-10-mvapich2-2.3"

        if "%clang@9.toss" in spec:
            compiler_install_path = "clang-9-mvapich2-2.3"

        if "%clang@10.toss" in spec:
            compiler_install_path = "clang-10-mvapich2-2.3"

        if "%clang@11.toss" in spec:
            compiler_install_path = "clang-11-mvapich2-2.3"

        if "%clang@toss" in spec:
            compiler_install_path = "clang-11-mvapich2-2.3"

        if "%xl_r@bgq" in spec:
            env['LLNL_CHECK_COMPILE_LINE'] = "None"
            env['CC']="/usr/local/bin/mpixlc_r-fastmpi" # set for detect to run on front end
            env['CXX']="/usr/local/bin/mpixlcxx_r-fastmpi"
            compiler_install_path = "xlc-12.1-opt"

        if "%clang@bgq" in spec:
            env['LLNL_CHECK_COMPILE_LINE'] = "None"
            env['CC']="/usr/apps/gnu/clang/2017.06.06/bin/mpiclang-fastmpi" # set for detect to run on front end
            env['CXX']="/usr/apps/gnu/clang/2017.06.06/bin/mpiclang++11-fastmpi"
            compiler_install_path = "bgqclang++11-opt"

        if "%gcc@bgq" in spec:
            env['LLNL_CHECK_COMPILE_LINE'] = "None"
            env['CC']="/usr/local/tools/compilers/ibm/mpicc-4.8.4-fastmpi" # set for detect to run on front end
            env['CXX']="/usr/local/tools/compilers/ibm/mpicxx-4.8.4-fastmpi"
            compiler_install_path = "g++11-472-opt"

        if "%intel@sandia" in spec or "%arm@sandia" in spec:
            #env['CC'] = "%s" % env['SPACK_CC']
            #env['CXX'] = "%s" % env['SPACK_CXX']
            #std_cmake_args.append("-DCMAKE_C_FLAGS=%s -g" % env['CFLAGS'])
            #std_cmake_args.append("-DCMAKE_CXX_FLAGS=%s -g" % env['CXXFLAGS'])
            filter_file(r'typedef int LEOS_MPI_Comm;', 'struct ompi_communicator_t;\ntypedef ompi_communicator_t * LEOS_MPI_Comm;', 'interfaces/C/LEOS_types_C.h') # openmpi
            std_cmake_args.append("-DENABLE_FIND_MPI=OFF")
            std_cmake_args.append("-DMPI_FOUND=TRUE")


        if "%cce" in spec:
            #std_cmake_args.append("-DCMAKE_C_FLAGS=%s -g" % env['CFLAGS'])
            #std_cmake_args.append("-DCMAKE_CXX_FLAGS=%s -g" % env['CXXFLAGS'])
            #std_cmake_args.append("-DUSE_BOOST=OFF") # Currently doesnt like boost
            #std_cmake_args.append("-DUSE_MPI=OFF") # MPI version needs boost
            env['FORMAT_TYPE_CHECKING'] = "RELAXED"
            env['FFLAGS'] += " -ef -e Z -f free"
            compiler_install_path = "cray-cce-9-llvm"

        if "%gcc" in spec and spec.target.family=='ppc64le':
            compiler_install_path = "gcc-8-cuda-10"

        if "%clang" in spec and spec.target.family=='ppc64le':
            env['FC'] = "xlf2003" # Need 2003
            compiler_install_path = "clang-11-cuda-11.2-gcc831"

        if "%clang@old.blueos":
            compiler_install_path = "clang-9-cuda-10"

        if "%clang@new.blueos":
            compiler_install_path = "clang-9-cuda-10_NEW"

        if "%xl" in spec and spec.target.family=='ppc64le':
            #if spec.satisfies('+cuda'):
            #    env['CXXFLAGS'] = "-std=c++11"
            compiler_install_path = "xlc-16-cuda-11.2-gcc831"

        if "%xl@old.blueos" in spec:
            compiler_install_path = "xlc-16-cuda-10"

        if "%xl@new.blueos" in spec:
            compiler_install_path = "xlc-16-cuda-10_NEW"

        if "%gcc@4.9.3.blueos" in spec:
            compiler_install_path = "gcc-4.9-cuda-10"

        if "%gcc@7.blueos" in spec:
            compiler_install_path = "gcc-7-cuda-10"

        if "%gcc@8.blueos" in spec:
            compiler_install_path = "gcc-8-cuda-10"

        if "%gcc@blueos" in spec:
            compiler_install_path = "gcc-8-cuda-10"

        if "%pgi@19.blueos" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-19-cuda-10"

        if "%pgi@20.blueos" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-20-cuda-10_NEW"

        if "%pgi@21.blueos" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-21-cuda-10"

        if "%pgi@blueos" in spec:
            std_cmake_args.append("-DCMAKE_BUILD_TYPE:STRING=Debug") # opt fails
            compiler_install_path = "pgi-21-cuda-10"

        os.system("mkdir -p leos-build")
        with working_dir('leos-build'):
            # Useful Cmake arguments: -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=c++
            # Note: -DCMAKE_BUILD_TYPE=Release (sets -O3), RelWithDebInfo (sets -g -O2), Debug (sets -g -O0)
            # ***** -DCMAKE_BUILD_TYPE OVERRIDES -DCMAKE_C_FLAGS= & -DCMAKE_CXX_FLAGS= SETTINGS!!!!
            # Append extra or overriding args at end by:  std_cmake_args.append("-DLEOS_USE_MPI=OFF")

            if spec.satisfies('+debug'):
                for narg in range(len(std_cmake_args)):
                    std_cmake_args[narg] = re.sub("^-DCMAKE_BUILD_TYPE:STRING=.*", "-DCMAKE_BUILD_TYPE:STRING=Debug", std_cmake_args[narg])
                    if re.match("^-DCMAKE_.*_FLAGS=.*", std_cmake_args[narg]):
                        std_cmake_args[narg] = re.sub(r"-O\d?\b", "-O0", std_cmake_args[narg])
                        std_cmake_args[narg] += " -g -O0"

            if spec.satisfies('+silo'):
                if spec.satisfies('@8.0.5p2:') or spec.satisfies('@8.0.5p1-cpp-lip'):
                    std_cmake_args.append("-DSILO=ON")
                    std_cmake_args.append("-DSILO_PATH=%s" % spec['silo'].prefix)
                    std_cmake_args.append("-DUSE_PDB=ON")
                else:
                    std_cmake_args.append("-DSILO=OFF")

            if spec.satisfies('+raja') or spec.satisfies('+cuda'):
                std_cmake_args.append("-DENABLE_RAJA=ON")
                std_cmake_args.append("-DRAJA_DIR=%s" % spec['raja'].prefix)
                if spec.satisfies('@:8.2.1'):
                    filter_file(r'0-9\]','0-9]+','../cmake/FindRAJA.cmake')

            if spec.satisfies('+cuda'):
                std_cmake_args.append("-DENABLE_RAJA_CUDA=ON")
                std_cmake_args.append("-DENABLE_CUDA=ON")
                if spec.satisfies('platform=linux') and spec.target.family=='ppc64le':
                    std_cmake_args.append("-DCUDA_TOOLKIT_ROOT_DIR=%s" % env['CUDA_HOME'])
                else:
                    std_cmake_args.append("-DCUDA_TOOLKIT_ROOT_DIR=/opt/cudatoolkit-8.0")
                if spec.satisfies('@:8.2.5'):
                    std_cmake_args.append("-DENABLE_CNMEM=ON")
                    std_cmake_args.append("-DCNMEM_DIR=%s" % spec['cnmem'].prefix)
                else:
                    std_cmake_args.append("-DENABLE_CNMEM=OFF")
                    filter_file(r' NOT ENABLE_CUDA', ' NOT ENABLE_FOO', '../interfaces/CMakeLists.txt')

            if spec.satisfies('+umpire') or spec.satisfies('+cuda'):
                std_cmake_args.append("-DENABLE_UMPIRE_IDS=ON")
                std_cmake_args.append("-DENABLE_UMPIRE=ON")
                std_cmake_args.append("-DUMPIRE_DIR=%s" % spec['umpire'].prefix)
                std_cmake_args.append('-DUMPIRE_LIBRARIES=%s' % spec['umpire'].prefix.lib)
                if spec.satisfies('@:8.1.999') and spec.satisfies('^umpire@0.3.4:'):
                    filter_file(r'umpire_op umpire_resource umpire_strategy umpire_util umpire_tpl_judy', '', '../cmake/FindUMPIRE.cmake')
                if spec.satisfies('^umpire@3.0.0:'):
                    filter_file(r'0.0', '3.0', '../cmake/FindUMPIRE.cmake')
            else:
                std_cmake_args.append("-DENABLE_UMPIRE_IDS=OFF")
                std_cmake_args.append("-DENABLE_UMPIRE=OFF")

            std_cmake_args.append("-DENABLE_FILTERS=OFF")
            std_cmake_args.append("-DUSE_LIBXML2=OFF")
            std_cmake_args.append("-DUSE_NOXML=ON")

            # Common cmake configuration args for this package
            std_cmake_args.append("-DCMAKE_C_COMPILER=%s" % env['CC'])
            std_cmake_args.append("-DCMAKE_CXX_COMPILER=%s" % env['CXX'])
            std_cmake_args.append("-DCMAKE_Fortran_COMPILER=%s" % env['FC'])
            std_cmake_args.append("-DCMAKE_INSTALL_PREFIX=%s" % prefix)
            if spec.satisfies('@8.2.0:') and not (spec.satisfies('%pgi') and 'toss' in str(spec.compiler)):
                if spec.satisfies("@8.3.6:"):
                    std_cmake_args.append("-DENABLE_CEREAL=ON")
                std_cmake_args.append("-DUSE_CEREAL=ON")
                std_cmake_args.append("-DUSE_BOOST=OFF")
            else:
                if spec.satisfies("@8.3.6:"):
                    std_cmake_args.append("-DENABLE_CEREAL=OFF")
                std_cmake_args.append("-DUSE_CEREAL=OFF")
                std_cmake_args.append("-DBOOST_ROOT=%s" % spec['boost'].prefix)
                std_cmake_args.append("-DUSE_BOOST=ON")
            std_cmake_args.append("-DENABLE_MPI=ON")
            std_cmake_args.append("-DUSE_MPI=ON")
            std_cmake_args.append("-DENABLE_PYTHON=OFF")
            std_cmake_args.append("-DENABLE_ZFP=OFF")
            std_cmake_args.append("-DBUILD_TESTS=OFF")
            std_cmake_args.append("-DBUILD_EXAMPLES=OFF")
            std_cmake_args.append("-DBUILD_EXAMPLES_P=OFF")
            if spec.satisfies('@8.1.0:'):
                std_cmake_args.append("-DBUILD_TOOLS=OFF")
                std_cmake_args.append("-DENABLE_V7=ON")
                std_cmake_args.append("-DSTATIC_FILTERS=ON")
                std_cmake_args.append("-DENABLE_COMPRESSION=OFF")
                std_cmake_args.append("-DENABLE_SM_LIP=OFF")
            else:
                std_cmake_args.append("-DBUILD_TOOLS=ON")
            std_cmake_args.append("-DENABLE_C_LIP=ON")
            std_cmake_args.append("-DENABLE_CPP_LIP=ON")
            std_cmake_args.append("-DLIP_SHARED_MEMORY=OFF")
            std_cmake_args.append("-DENABLE_ASCMEMORY=OFF")

            std_cmake_args.append("-DUSE_HDF=ON")
            std_cmake_args.append("-DHDF_ROOT=%s" % spec['hdf5'].prefix)
            std_cmake_args.append("-DHDF5_ROOT=%s" % spec['hdf5'].prefix)
            std_cmake_args.append("-DCONFIG=spack")

            if spec.satisfies('@8.0.5:'):
                std_cmake_args.append("-DENABLE_FORTRAN=ON")
                std_cmake_args.append("-DENABLE_C=ON")

            # FIXME: Modify the configure line to suit your build system here.
            print("NOTE: cc=%s  c++=%s  fc=%s" % (env['SPACK_CC'],env['SPACK_CXX'],env['SPACK_FC']))
            # Traditionally built with -DCMAKE_BUILD_TYPE=Release (-O3) not RelWithDebInfo (-O2)
            # Optimization level discrepancies show up as test curve differences

            if spec.satisfies("%intel@sandia") or spec.satisfies("%arm@sandia"):
                env['SYS_TYPE'] = 'sandia'
                if spec.satisfies("%arm@sandia"):
                    env['FFLAGS'] = "%s -std=f2003" % env['FFLAGS']
            elif spec.satisfies('platform=cray') or spec.satisfies("%intel@trinity") or spec.satisfies("%cce@trinity"):
                env['SYS_TYPE'] = 'trinity'
            cmake("..",
                  "-DEOS_DATA_ROOT_DIR=%s" % self.datadir,
                  *std_cmake_args)

            if spec.satisfies('~silo') or spec.satisfies('@:8.0.5p2'):
                if any(self.spec.compiler.satisfies(c) for c in self.lanl_compilers + self.sandia_compilers):
                    filter_file(r'Bdynamic', 'Bstatic', 'LEOSPACT/CMakeFiles/detect.dir/link.txt')

                #if spec.satisfies('%intel platform=cray'):
                #    filter_file(r'gnu11', 'gnu11 -mavx -axCORE-AVX2,MIC-AVX512 -mtune=core-avx2', 'LEOSPACT/CMakeFiles/detect.dir/build.make')

            # FIXME: Add logic to build and install here
            make("VERBOSE=1")
            make("install")
