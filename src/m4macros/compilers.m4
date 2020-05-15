# =======================================================================
# Macros dealing with c/c++ compilers and options
# =======================================================================
# Note that the options for MPICCCOMPILERS, MPICXXCOMPILERS, CCCOMPILERS,
# and CXXCOMPILERS are all set in mpi.m4.

# =======================================================================
# SETUP_COMPILERS_THE_WAY_I_WANT sets the CC and CXX variables.
# from a command line option, if evironment not set.
# =======================================================================
AC_DEFUN([SETUP_COMPILERS_THE_WAY_I_WANT],[
# Begin macro SETUP_COMPILERS_THE_WAY_I_WANT.

AC_SUBST(OSNAME)

AC_SUBST(CC)
AC_SUBST(CXX)
AC_SUBST(FORT)
AC_SUBST(FORTLINK)

AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(LD)
AC_SUBST(DYLIBEXT)
AC_SUBST(SHLIBEXT)

AC_SUBST(PYTHONCC)
AC_SUBST(PYTHONCXX)

AC_SUBST(GCCXMLCC)
AC_SUBST(GCCXMLCXX)

AC_SUBST(CMAKECC)
AC_SUBST(CMAKECXX)

AC_SUBST(PARMETISCC)
AC_SUBST(MPI4PYCC)

AC_SUBST(CXXCOMPILERTYPE)
AC_SUBST(JAMOPTS)
AC_SUBST(JAMTOOLSET)
AC_SUBST(JAMTOOLSETOPTS)
AC_SUBST(BOOSTEXT)

AC_SUBST(LD)
AC_SUBST(LDFLAGS)
AC_SUBST(LDRPATH)
AC_SUBST(LDINSTALLNAME)
AC_SUBST(LIBTARGETFLAGS)
AC_SUBST(LIBS)
AC_SUBST(DEPENDRULES)
AC_SUBST(DYNLIBFLAG)
AC_SUBST(SHAREDFLAG)
AC_SUBST(LDPASSTHROUGH)

AC_SUBST(CXXFLAGS)
AC_SUBST(EXTRAFLAGS)
AC_SUBST(EXTRAINCLUDES)
AC_SUBST(FORTFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(MPICCFLAGS)
AC_SUBST(MPICXXFLAGS)
AC_SUBST(DEPFLAG)

AC_SUBST(PYTHONCFLAGS)
AC_SUBST(PYTHONCONFFLAGS)
AC_SUBST(NUMPYFLAGS)
AC_SUBST(NUMPYCFLAGS)
AC_SUBST(HDF5FLAGS)

AC_SUBST(FFTWFLAGS)

PYTHONCONFFLAGS=
LIBTARGETFLAGS=
JAMTOOLSETOPTS=
LDINSTALLNAME="-o"
LDRPATH=
FORTLINK=
NUMPYFLAGS=
EXTRAINCLUDES=
NUMPYCFLAGS=

# =======================================================================
# Selection of approved compiler sets for Spheral++.
# =======================================================================
AC_MSG_CHECKING(for compilers)
AC_ARG_WITH(compilers,
[  --with-compilers=ARG ..................... (gnu,clang,clang-ibm,vacpp,intel,pgi) choose a compiler suite],
[
   COMPILERS=$withval
],
[
   COMPILERS="gnu"
]
)
AC_MSG_RESULT($COMPILERS)

OSNAME=`uname -s`
AC_CHECK_PROG(GCCTEST, gcc, `which gcc`, nope)
AC_CHECK_PROG(MPICXXTEST, mpicxx, `which mpicxx`, nope)
case $COMPILERS in
   gnu)
      if test $OSNAME = "Linux" -a "$GCCTEST" != "nope"; then
         CC=gcc
         CXX=g++
         FORT=gfortran
         MPICC=mpicc
         #MPICCFLAGS="-cc=$CC"
         if test "$MPICXXTEST" != "nope"; then
            MPICXX=mpicxx
            #MPICXXFLAGS="-cxx=$CXX"
         else
            MPICXX=mpig++
            #MPICXXFLAGS="-cc=$CXX"
         fi
         CMAKECC=$CC
         CMAKECXX=$CXX
         GCCXMLCC=$CMAKECC
         GCCXMLCXX=$CMAKECXX
         PYTHONCC=$CC
         PYTHONCXX=$CXX
         PARMETISCC=$MPICC
         MPI4PYCC=$MPICC
         CXXFLAGS+=" -std=c++11 -march=native"

      else
         CC=gcc
         CXX=g++
         FORT=gfortran
         MPICC=mpicc
         MPICXX=mpic++
         CMAKECC=$CC
         CMAKECXX=$CXX
         GCCXMLCC=$CMAKECC
         GCCXMLCXX=$CMAKECXX
         PYTHONCC=$CC
         PYTHONCXX=$CXX
         PARMETISCC=$MPICC
         MPI4PYCC=$MPICC
         CXXFLAGS+=" -std=c++11 -march=native"
         if test $OSNAME = "Darwin"; then
           CXXFLAGS+=" -mmacosx-version-min=10.7 -stdlib=libc++"
         fi
      fi
      ;;

   clang)
      CC=clang
      CXX=clang++
      FORT=gfortran
      MPICC=mpicc
      MPICXX=mpicxx
      MPICCFLAGS=
      MPICXXFLAGS=
      CMAKECC=clang
      CMAKECXX=clang++
      GCCXMLCC=$CMAKECC
      GCCXMLCXX=$CMAKECXX
      PYTHONCC=$CC
      PYTHONCXX=$CXX
      PARMETISCC=$MPICC
      MPI4PYCC=$MPICC
      CXXFLAGS+=" -std=c++11 -Wno-undefined-var-template -march=native"
      if test $OSNAME = "Darwin"; then
        CXXFLAGS+=" -mmacosx-version-min=10.7 -stdlib=libc++"
      fi
      ;;

   clang-ibm)
      CC=clang
      CXX=clang++
      FORT=gfortran
      MPICC=mpiclang-gpu
      MPICXX=mpiclang++-gpu
      MPICCFLAGS=
      MPICXXFLAGS=
      CMAKECC=clang
      CMAKECXX=clang++
      GCCXMLCC=$CMAKECC
      GCCXMLCXX=$CMAKECXX
      PYTHONCC=gcc
      PYTHONCXX=g++
      PARMETISCC=$MPICC
      MPI4PYCC=$MPICC
      CXXFLAGS+=" -std=c++11 -DEIGEN_DONT_VECTORIZE"
      ;;

   gcc-bg)
      CC=mpigcc
      CXX=mpig++
      FORT=mpigfortran
      MPICC=mpigcc
      MPICXX=mpig++
      MPICCFLAGS=
      MPICXXFLAGS=
      CMAKECC=gcc
      CMAKECXX=g++
      GCCXMLCC=$CMAKECC
      GCCXMLCXX=$CMAKECXX
      PYTHONCC=bggcc
      PYTHONCXX=bgg++
      PARMETISCC=$MPICC
      MPI4PYCC=mpixlc_r
      CXXFLAGS+=" -std=c++11 -DEIGEN_DONT_VECTORIZE"
      #LDFLAGS+=" -dynamic"
      HDF5FLAGS+=" --enable-shared=no --enable-static=yes --enable-static-exec=yes"
      ;;

   clang-bg)
      CC=bgclang
      CXX=bgclang++11
      FORT=gfortran-4.7.2-fastmpi
      MPICC=mpiclang
      MPICXX=mpiclang++11
      MPICCFLAGS=
      MPICXXFLAGS=
      CMAKECC=gcc
      CMAKECXX=g++
      GCCXMLCC=$CMAKECC
      GCCXMLCXX=$CMAKECXX
      PYTHONCC=gcc
      PYTHONCXX=g++
      PARMETISCC=$MPICC
      MPI4PYCC=mpixlc_r
      CXXFLAGS+=" -std=c++11 -DEIGEN_DONT_VECTORIZE"
      HDF5FLAGS+=" --enable-shared=no --enable-static=yes --enable-static-exec=yes"
      ;;

   vacpp-bg)
      CC=xlc_r
      CXX=xlC_r
      FORT=mpixlf-fastmpi
      MPICC=mpixlc_r-fastmpi
      MPICXX=mpixlcxx_r-fastmpi
      MPICCFLAGS=
      MPICXXFLAGS=
      CMAKECC=gcc
      CMAKECXX=g++
      GCCXMLCC=$CMAKECC
      GCCXMLCXX=$CMAKECXX
      PYTHONCC=gcc
      PYTHONCXX=g++
      PARMETISCC=$MPICC
      MPI4PYCC=$MPICC
      CXXFLAGS+=" -qlanglvl=extended0x -DEIGEN_DONT_ALIGN -DEIGEN_DONT_VECTORIZE "
      HDF5FLAGS+=" --enable-shared=no --enable-static=yes --enable-static-exec=yes"
      ;;

   vacpp)
      CC=xlc_r
      CXX=xlC_r
      MPICC=mpicc
      MPICXX=mpic++ 
      CMAKECC=gcc
      CMAKECXX=g++
      GCCXMLCC=/usr/tcetmp/packages/gcc/gcc-4.9.3/bin/gcc
      GCCXMLCXX=/usr/tcetmp/packages/gcc/gcc-4.9.3/bin/g++
      PYTHONCC=/usr/tcetmp/packages/gcc/gcc-4.9.3/bin/gcc
      PYTHONCXX=/usr/tcetmp/packages/gcc/gcc-4.9.3/bin/g++
      PARMETISCC=$MPICC
      MPI4PYCC=$MPICC
      CFLAGS+=" "
      CXXFLAGS+=" -std=c++11 -qnoinline -qnoxlcompatmacros -qmaxmem=16384  -DEIGEN_DONT_ALIGN -DEIGEN_DONT_VECTORIZE "
      ;;

   intel)
      CC=icc
      CXX=icpc
      FORT=ifort
      MPICC=mpicc  # mpiicc
      MPICXX=mpic++ # mpiicpc
      PYTHONCC=$CC
      PYTHONCXX=$CXX
      CMAKECC=$CC
      CMAKECXX=$CXX
      GCCXMLCC=gcc
      GCCXMLCXX=g++
      CMAKECC=gcc
      CMAKECXX=g++
      PARMETISCC=$MPICC
      MPI4PYCC=$MPICC
      CXXFLAGS+=" -std=c++11 -flto"
      NUMPYFLAGS="--fcompiler=intelem"
      NUMPYCFLAGS="CFLAGS=-no-ip"
      ;;

   pgi)
      CC=pgcc
      CXX=pgCC
      FORT=mpif90
      MPICC=mpicc
      MPICXX=mpicxx
      PYTHONCC=$CC
      PYTHONCXX=$CXX
      CMAKECC=$CC
      CMAKECXX=$CXX
      GCCXMLCC=
      GCCXMLCXX=
      CMAKECC=pgcc
      CMAKECXX=pbCC
      PARMETISCC=$MPICC
      MPI4PYCC=$MPICC
      NUMPYFLAGS=
      # 111  - statement is unreachable
      # 186  - pointless comparison of unsigned integer with zero
      # 611  - overloaded virtual function
      CXXFLAGS+=" --display_error_number --diag_suppress 111 --diag_suppress 186 --diag_suppress 611"
      PYTHONCFLAGS+=" -fPIC"
      PYTHONCONFFLAGS+=" --enable-shared"
      ;;

   *)
      CC=gcc
      CXX=g++
      MPICC=mpicc
      MPICXX=mpiCC
      PYTHONCC=$CC
      PYTHONCXX=$CXX
      CMAKECC=$CC
      CMAKECXX=$CXX
      GCCXMLCC=$CC
      GCCXMLCXX=$CXX
      PARMETISCC=$MPICC
      MPI4PYCC=$MPICC
      #PYTHONCONFFLAGS="--with-gcc=$PYTHONCC"
      ;;

esac

# If we're not compiling with MPI, default the MPI compilers
# to CC and CXX.
if test $MPIENABLED = "no"; then
  MPICC=$CC
  MPICXX=$CXX
  MPICCFLAGS=
  MPICXXFLAGS=
fi

# On Mac with Xcode we need to change some link options
if test $OSNAME = "Darwin"; then
  LDFLAGS+=" -stdlib=libc++ -mmacosx-version-min=10.9"
fi

## On 64 bit Darwin we have to diddle the python configure 
#if test $OSNAME = "Darwin"; then
#   PYTHONCONFFLAGS="--enable-framework=$prefix DESTDIR=$prefix"
#   #PYTHONCONFFLAGS="'MACOSX_DEPLOYMENT_TARGET=10.5' --enable-framework=$prefix --enable-universalsdk" # --disable-toolbox-glue"
#fi

# =======================================================================
# MPI compilers
# =======================================================================
AC_MSG_CHECKING(for MPICC)
AC_ARG_WITH(MPICC,
[  --with-MPICC=ARG ......................... manually set the mpi C++ compiler to ARG],
[
   MPICC=$withval
   AC_MSG_RESULT($MPICC)
],
[
   AC_MSG_RESULT($MPICC)
]
)

AC_MSG_CHECKING(for MPICXX)
AC_ARG_WITH(MPICXX,
[  --with-MPICXX=ARG ........................ manually set the mpi C++ compiler to ARG],
[
   MPICXX=$withval
   AC_MSG_RESULT($MPICXX)
],
[
   AC_MSG_RESULT($MPICXX)
]
)

# =======================================================================
# Generic serial compilers
# =======================================================================
AC_MSG_CHECKING(for CC)
AC_ARG_WITH(CC,
[  --with-CC=ARG ............................ manually set C compiler to ARG],
[
   CC=$withval
   AC_MSG_RESULT($CC)
],
[
   AC_MSG_RESULT($CC)
]
)

AC_MSG_CHECKING(for CXX)
AC_ARG_WITH(CXX,
[  --with-CXX=ARG ........................... manually set C++ compiler to ARG],
[
   CXX=$withval
   AC_MSG_RESULT($CXX)
],
[
   AC_MSG_RESULT($CXX)
]
)

# =======================================================================
# Python compilers
# =======================================================================
AC_MSG_CHECKING(for python-CC)
AC_ARG_WITH(python-CC,
[  --with-python-CC=ARG ..................... manually set the CC compiler for python],
[
   PYTHONCC=$withval
   AC_MSG_RESULT($PYTHONCC)
],
[
   AC_MSG_RESULT($PYTHONCC)
]
)

AC_MSG_CHECKING(for python-CXX)
AC_ARG_WITH(python-CXX,
[  --with-python-CXX=ARG .................... manually set the CXX compiler for python],
[
   PYTHONCXX=$withval
   AC_MSG_RESULT($PYTHONCXX)
],
[
   AC_MSG_RESULT($PYTHONCXX)
]
)

# =======================================================================
# mpi4py compilers
# =======================================================================
AC_MSG_CHECKING(for MPI4PYCC)
AC_ARG_WITH(MPI4PYCC,
[  --with-MPI4PYCC=ARG ...................... manually set the CC compiler for mpi4py],
[
   MPI4PYCC=$withval
   AC_MSG_RESULT($MPI4PYCC)
],
[
   AC_MSG_RESULT($MPI4PYCC)
]
)

# =======================================================================
# GCCXML compilers
# =======================================================================
AC_MSG_CHECKING(for GCCXMLCC)
AC_ARG_WITH(GCCXMLCC,
[  --with-GCCXMLCC=ARG ...................... manually set the CC compiler for building gccxml],
[
   GCCXMLCC=$withval
   AC_MSG_RESULT($GCCXMLCC)
],
[
   AC_MSG_RESULT($GCCXMLCC)
]
)

AC_MSG_CHECKING(for GCCXMLCXX)
AC_ARG_WITH(GCCXMLCXX,
[  --with-GCCXMLCXX=ARG ..................... manually set the CXX compiler for building gccxml],
[
   GCCXMLCXX=$withval
   AC_MSG_RESULT($GCCXMLCXX)
],
[
   AC_MSG_RESULT($GCCXMLCXX)
]
)


# =======================================================================
# Cmake compilers
# =======================================================================
AC_MSG_CHECKING(for cmake-CC)
AC_ARG_WITH(cmake-CC,
[  --with-cmake-CC=ARG ...................... manually set the CC compiler for cmake],
[
   CMAKECC=$withval
   AC_MSG_RESULT($CMAKECC)
],
[
   AC_MSG_RESULT($CMAKECC)
]
)

AC_MSG_CHECKING(for cmake-CXX)
AC_ARG_WITH(cmake-CXX,
[  --with-cmake-CXX=ARG ..................... manually set the CXX compiler for cmake],
[
   CMAKECXX=$withval
   AC_MSG_RESULT($CMAKECXX)
],
[
   AC_MSG_RESULT($CMAKECXX)
]
)

# =======================================================================
# ParMETIS compilers
# =======================================================================
AC_MSG_CHECKING(for parmetis-CC)
AC_ARG_WITH(parmetis-CC,
[  --with-parmetis-CC=ARG ................... manually set the CC compiler for ParMETIS],
[
   PARMETISCC=$withval
   AC_MSG_RESULT(PARMETISCC)
],
[
   AC_MSG_RESULT($PARMETISCC)
]
)

# =======================================================================
# Fortran compiler
# =======================================================================
AC_MSG_CHECKING(for fortran)
AC_ARG_WITH(fortran,
[  --with-fortran=ARG ....................... manually set the fortran compiler],
[
   FORT=$withval
   AC_MSG_RESULT($FORT)
],
[
   AC_MSG_RESULT($FORT)
]
)

# ======================================================================
# Set up the linker to be used generating shared libs
# ======================================================================
AC_MSG_CHECKING(for LD)
AC_ARG_WITH(LD,
  [  --with-LD=ARG ............................ manually set linker to ARG],
    LD = $withval,
  [
  if test ! "$LD"; then 
    if test $MPIENABLED = "yes"; then
      LD="$MPICXX"
    else
      LD="$CXX"
    fi
  fi
  ]
)
AC_MSG_RESULT($LD)

# ======================================================================
# Try and detect what sort of C++ compiler we're using, and set the jam
# info appropriately.
# ======================================================================
AC_MSG_CHECKING(for Boost toolset)
CXXCOMPILERTYPE=
JAMTOOLSET=
BOOSTEXT=
SHAREDFLAG="-shared"
DEPFLAG="-E -M -w"

# Check for the GNU compiler and version
cat > .cxxtype.cc << EOF
#ifdef __GNUC__
yes;
#endif
EOF
$CXX -E .cxxtype.cc > .cxxtype.out
if test -n "`grep 'yes' .cxxtype.out`"; then
  CXXCOMPILERTYPE=GNU
  if test -n "`$CXX --version | grep '4\.0'`"; then
    COMPILERVERSION=40
  elif test -n "`$CXX --version | grep '4\.1'`"; then
    COMPILERVERSION=41
  elif test -n "`$CXX --version | grep '4\.2'`"; then
    COMPILERVERSION=42
  elif test -n "`$CXX --version | grep '4\.3'`"; then
    COMPILERVERSION=43
  fi
fi
rm -f .cxxtype.cc .cxxtype.out

# Check for the intel compiler and version
cat > .cxxtype.cc << EOF
#ifdef __INTEL_COMPILER
yes;
__INTEL_COMPILER
#endif
EOF
$CXX -E .cxxtype.cc > .cxxtype.out
if test -n "`grep 'yes' .cxxtype.out`"; then
  CXXCOMPILERTYPE=INTEL
  if test -n "`grep '800' .cxxtype.out`"; then
    COMPILERVERSION=800
  elif test -n "`grep '810' .cxxtype.out`"; then
    COMPILERVERSION=810
  elif test -n "`grep '910' .cxxtype.out`"; then
    COMPILERVERSION=910
  fi
fi
rm -f .cxxtype.cc .cxxtype.out

# Check for the PGI compiler and version
cat > .cxxtype.cc << EOF
#ifdef __PGI
yes;
#endif
EOF
$CXX -E .cxxtype.cc > .cxxtype.out
if test -n "`grep 'yes' .cxxtype.out`"; then
  CXXCOMPILERTYPE=PGI
fi
rm -f .cxxtype.cc .cxxtype.out

# Check for KCC
cat > .cxxtype.cc << EOF
#ifdef __KAI__
yes;
#endif
EOF
$CXX -E .cxxtype.cc > .cxxtype.out
if test -n "`grep 'yes' .cxxtype.out`"; then
  CXXCOMPILERTYPE=KAI
fi
rm -f .cxxtype.cc .cxxtype.out

# Check for xlC
cat > .cxxtype.cc << EOF
#ifdef __IBMCPP__
yes;
#endif
EOF
$CXX -E .cxxtype.cc > .cxxtype.out
if test -n "`grep 'yes' .cxxtype.out`"; then
  CXXCOMPILERTYPE=VACPP
fi
rm -f .cxxtype.cc .cxxtype.out

# Set the flag for passing arguments to the linker.
LDPASSTHROUGH=""
if test $CXXCOMPILERTYPE = "GNU" -o $CXXCOMPILERTYPE = "INTEL" -o $CXXCOMPILERTYPE = "VACPP"; then
  LDPASSTHROUGH="-Wl,"
fi

# A few AIX specializations.
if test "$OSNAME" = "AIX"; then
  CFLAGS="$CFLAGS -maix64"
  LDFLAGS="$LDFLAGS -Wl,-bmaxstack:0x10000000 -Wl,-bmaxdata:0x70000000"
  SHAREDFLAG="-shared -Wl,-G -Wl,-bexpall -Wl,-bexpfull"

elif test "$OSNAME" = "Linux"; then # -a "$CXXCOMPILERTYPE" != "INTEL"; then
  # On the gnu linker we can throw the rpath flag to avoid having to set the LD_LIBRARY_PATH
  # variable.
  LDRPATH="$LDRPATH ${LDPASSTHROUGH}-rpath=\$(libdir)"

elif test "$OSNAME" = "Darwin"; then
  LDRPATH="$LDRPATH ${LDPASSTHROUGH}-rpath \$(libdir)"
  LDINSTALLNAME="-install_name @rpath/\${@} -o"

#   # On Mac OS X Darwin, you install libraries with an "dylib_install_name" flag to avoid
#   # using rpath or DYLD_LIBRARY_PATH nonsense when linking to these libraries.
#    LIBTARGETFLAGS="$LIBTARGETFLAGS ${LDPASSTHROUGH} -dylib_install_name \$(LIBDIR)/\$(LIBTARGET)"

fi

# We have to make sure the LIBDIR directory exists, or the -rpath flag will cause the
# compiler tests to fail.
# if test ! -e $LIBDIR; then
#   mkdir -p $LIBDIR
# fi

# I'm not sure why we need this flag, but without this the boost::shared_ptr is 
# hanging in a thread lock.
CXXFLAGS="$CXXFLAGS -DBOOST_DISABLE_THREADS"

# Default the DYNLIBFLAG to be the SHAREDFLAG, which is true for everybody but
# Apple.  >:(
DYNLIBFLAG="$SHAREDFLAG"

case $CXXCOMPILERTYPE in 
GNU)
  if test "$OSNAME" = "Darwin"; then
     FORTLINK=" -L/opt/local/lib -lf2c"
     CFLAGS="$CFLAGS -fPIC"
     CXXFLAGS="$CXXFLAGS -fPIC -DHAVE_XCPT -DGNUCXX"
     FORTFLAGS="$FORTFLAGS -fPIC"
     DYNLIBFLAG="-dynamiclib -undefined suppress -flat_namespace"
     SHAREDFLAG="-bundle -undefined suppress -flat_namespace"
     JAMTOOLSET=darwin
     LDFLAGS="$LDFLAGS -undefined dynamic_lookup"
  else
     FORTLINK=" -lgfortran"
     CFLAGS="$CFLAGS -fpic -fexceptions"
     CXXFLAGS="$CXXFLAGS -fpic -fexceptions -DHAVE_XCPT -DGNUCXX"
     FORTFLAGS="$FORTFLAGS -fpic" #  -ffpe-trap=invalid,zero,overflow,underflow,denormal
     JAMTOOLSET=gcc
     BOOSTEXT="-$JAMTOOLSET$COMPILERVERSION"
     if test "$OSNAME" = "AIX"; then
        LDFLAGS="$LDFLAGS -shared-libgcc -lstdc++ -lsupc++"
     fi
  fi
  DEPFLAG="-MM"
  ;;
INTEL)
  # The -wd654 suppresses the "virtual methods partially overridden warning", which lots of Spheral++ code
  # emits by design.
  CFLAGS="$CFLAGS -fpic" #  -wd654"
  CXXFLAGS="$CXXFLAGS -fpic" #  -wd654"
  FORTFLAGS="$FORTFLAGS -fpic"
  #LIBS="$LIBS -lrt -lcxa -lirc"
  JAMTOOLSET="intel-linux"
  BOOSTEXT="-il"
#  if test "$COMPILERVERSION" = "800"; then
#    JAMTOOLSET="intel-8.0"
#  elif test "$COMPILERVERSION" = "810"; then
#    JAMTOOLSET="intel-8.1"
#  elif test "$COMPILERVERSION" = "910"; then
#    JAMTOOLSET="intel-9.1"
#  fi
  ;;
PGI)
  CFLAGS="$CFLAGS -fPIC"
  CXXFLAGS="$CXXFLAGS -fPIC"
  FORTFLAGS="$FORTFLAGS -fPIC"
  ;;
KAI)
  CXXFLAGS="$CXXFLAGS --restrict -DHAVE_XCPT"
  FORTFLAGS="$FORTFLAGS -fpic"
  if test "$OSNAME" = "AIX"; then
    CXXFLAGS="$CXXFLAGS -qstaticinline -qnofullpath"
  fi
  SHAREDFLAG=""
  DEPFLAG="-M --no_code_gen"
  CFLAGS="$CFLAGS -g"
  JAMTOOLSET=kcc
  BOOSTEXT="-$JAMTOOLSET"
  ;;
VACPP)
  FORTFLAGS="$FORTFLAGS " 
  SHAREDFLAG="$SHAREDFLAG -G -qmkshrobj"
  DEPFLAG="-M"
  #DEPFLAG="-M -E"
  #DEPENDRULES="dependrules.aix"
  CFLAGS="$CFLAGS -g"
  JAMTOOLSET=vacpp 
  BOOSTEXT="-xlc"
  ;;
esac

# We seem to be always getting the mulit-threaded thingy now with boost.
BOOSTEXT="$BOOSTEXT-mt"

# These boost libraries wreak havoc on AIX and Darwin, so turn 'em off.
if test "$OSNAME" = "AIX" -o "$OSNAME" = "Darwin"; then
  JAMOPTS="$JAMOPTS --without-test --without-serialization --without-wave --without-regex --without-iostreams"
fi

# On AIX we have to modify the toolset options passed to jam.
if test "$OSNAME" = "AIX"; then
  JAMTOOLSETOPTS=" : 3.4 : gcc-3.4.3 "
fi

# On Darwin we build with frameworks, which means passing more flags to the
# jam build.
# if test "$OSNAME" = "Darwin"; then
#   JAMTOOLSETOPTS="$JAMTOOLSETOPTS : : $prefix/Python.framework "
# fi

# Shared library extensions.
DYLIBEXT=so
SHLIBEXT=so

# If we're on Mac OS X, we need to use the 'Darwin' toolset, and we need to 
# use Apple's screwy shared library extension.
if test "$OSNAME" = "Darwin"; then
  JAMTOOLSET=darwin
  DYLIBEXT=dylib
fi

AC_MSG_RESULT($JAMTOOLSET)

# =======================================================================
# openmp or not
# =======================================================================
AC_MSG_CHECKING(for --without-openmp)
AC_ARG_WITH(openmp,
[  --without-openmp ......................... build without OpenMP],
[
   AC_MSG_RESULT(yes)
   PYTHONPKGS+=" OpenMP"
],
[
   AC_MSG_RESULT(no)
   PYTHONPKGS+=" OpenMP"
   FFTWFLAGS+=" --enable-openmp"
   if test $CXXCOMPILERTYPE = "VACPP"; then
      CXXFLAGS+=" "
      EXTRAFLAGS+="-qsmp=omp -qoffload -I/usr/tcetmp/packages/cuda-9.0.176/include    "
   else
      CXXFLAGS+=" -fopenmp"
    #  CXXFLAGS+=" "
    #  EXTRAFLAGS+=" -qsmp=omp -qoffload -I/usr/tcetmp/packages/cuda-9.0.176/include    "
   fi
]
)

# =======================================================================
# UVM (unified memory for GPUs and such)
# =======================================================================
AC_MSG_CHECKING(for uvm)
AC_ARG_WITH(uvm,
[  --with-uvm ............................... enable unified memory (only for use with OpenMP)],
[
   AC_MSG_RESULT(yes)
   EXTRAFLAGS+=" -DUSE_UVM"
   EXTRAFLAGS+=" -I/usr/tcetmp/packages/cuda-9.0.184/include -fopenmp-targets=nvptx64-nvidia-cuda -fopenmp-implicit-declare-target"
],
[
   AC_MSG_RESULT(no)
]
)

# =======================================================================
# gprof
# =======================================================================
AC_MSG_CHECKING(for --with-gprof)
AC_ARG_WITH(gprof,
[  --with-gprof ............................. compile with gprof stuff turned on],
[
  AC_MSG_RESULT(yes)
  if test "$CXXCOMPILERTYPE" = "GNU"; then
    PYTHONCFLAGS+=" -pg"
    CFLAGS+=" -pg"
    CXXFLAGS+=" -pg"
    LDFLAGS+=" -pg"
  elif test "$CXXCOMPILERTYPE" = "INTEL"; then
    PYTHONCFLAGS+=" -p -g"
    CFLAGS+=" -p -g"
    CXXFLAGS+=" -p -g"
    LDFLAGS+=" -pg"
  fi
],
[
  AC_MSG_RESULT(no)
])

# =======================================================================
# gperftools
# =======================================================================
AC_MSG_CHECKING(for --with-gperftools)
AC_ARG_WITH(gperftools,
[  --with-gperftools ........................ compile linking with gperftools (optionally specify link path to libprofiler)],
[
  PYTHONPKGS+=" Gperftools"
  if test $withval = "yes"; then
    LDFLAGS+=" -L$withval -lprofiler"
    AC_MSG_RESULT($withval)
  else
    LDFLAGS+=" -lprofiler"
    AC_MSG_RESULT(yes)
  fi
],
[
  AC_MSG_RESULT(no)
])

])
