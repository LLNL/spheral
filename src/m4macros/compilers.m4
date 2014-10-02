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
AC_SUBST(FORTFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(MPICCFLAGS)
AC_SUBST(MPICXXFLAGS)
AC_SUBST(DEPFLAG)

AC_SUBST(PYTHONCFLAGS)
AC_SUBST(PYTHONCONFFLAGS)

PYTHONCONFFLAGS=
LIBTARGETFLAGS=
JAMTOOLSETOPTS=
LDINSTALLNAME="-o"
LDRPATH=
FORTLINK=

# =======================================================================
# Selection of approved compiler sets for Spheral++.
# =======================================================================
AC_MSG_CHECKING(for compilers)
AC_ARG_WITH(compilers,
[  --with-compilers=ARG ..................... (gnu,vacpp,intel,emsolve) choose a compiler suite],
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
AC_CHECK_PROG(GCC446TEST, gcc-4.4.6, `which gcc-4.4.6`, nope)
case $COMPILERS in
   gnu)
      if test $OSNAME = "Linux" -a "$GCC446TEST" != "nope"; then
         CC=gcc-4.4.6
         CXX=g++-4.4.6
         FORT=gfortran
         MPICC=mpicc
         MPICXX=mpig++
         MPICCFLAGS="-cc=$CC"
         MPICXXFLAGS="-cc=$CXX"
         CMAKECC=$CC
         CMAKECXX=$CXX
         GCCXMLCC=$CMAKECC
         GCCXMLCXX=$CMAKECXX
         PYTHONCC=$CC
         PYTHONCXX=$CXX
         PARMETISCC=$MPICC
         JAMTOOLSETOPTS=" : 4.4 : gcc-4.4.6 "

      elif test $OSNAME = "Linux" -a "$GCCTEST" != "nope"; then
         CC=gcc
         CXX=g++
         FORT=gfortran
         MPICC=mpicc
         MPICXX=mpiCC
         if test "$GCC333TEST" != "nope"; then
            CMAKECC=gcc-3.3.3
            CMAKECXX=g++-3.3.3
         else
            CMAKECC=$CC
            CMAKECXX=$CXX
         fi
         GCCXMLCC=$CMAKECC
         GCCXMLCXX=$CMAKECXX
         PYTHONCC=$CC
         PYTHONCXX=$CXX
         PARMETISCC=$MPICC

      elif test $OSNAME = "Darwin"; then
         CC=clang
         CXX=clang++
         FORT=gfortran
         MPICC=mpicc
         MPICXX=mpicxx
         MPICCFLAGS="-cc=clang"
         MPICXXFLAGS="-cxx=clang++"
         CMAKECC=clang
         CMAKECXX=clang++
         GCCXMLCC=$CMAKECC
         GCCXMLCXX=$CMAKECXX
         PYTHONCC=$CC
         PYTHONCXX=$CXX
         PARMETISCC=$MPICC

      else
         CC=gcc
         CXX=g++
         FORT=gfortran
         MPICC=mpicc # $SRCDIR/helpers/mpicc
         MPICXX=mpic++ # $SRCDIR/helpers/mpic++
         CMAKECC=$CC
         CMAKECXX=$CXX
         GCCXMLCC=$CMAKECC
         GCCXMLCXX=$CMAKECXX
         PYTHONCC=$CC
         PYTHONCXX=$CXX
         PARMETISCC=$MPICC

      fi
      ;;

   vacpp)
      CC=/usr/local/tools/compilers/ibm/xlc-8.0.0.12a
      CXX=/usr/local/tools/compilers/ibm/xlC-8.0.0.12a
      MPICC=/usr/local/tools/compilers/ibm/mpxlc-8.0.0.12a
      MPICXX=/usr/local/tools/compilers/ibm/mpxlC-8.0.0.12a 
      CMAKECC=$CC
      CMAKECXX=$CXX
      GCCXMLCC=gcc-3.2.3
      GCCXMLCXX=g++-3.2.3
      PYTHONCC=$CC
      PYTHONCXX=$CXX
      PARMETISCC=$MPICC
      ;;

   intel)
      CC=icc
      CXX=icpc
      FORT=ifort
      MPICC=mpiicc
      MPICXX=mpiicpc
      PYTHONCC=icc
      PYTHONCXX=icpc
      CMAKECC=$CC
      CMAKECXX=$CXX
      GCCXMLCC=gcc
      GCCXMLCXX=g++
      CMAKECC=gcc
      CMAKECXX=g++
      PARMETISCC=$MPICC
      ;;

   emsolve)
      CC=icc
      CXX=icpc
      MPICC=$SRCDIR/helpers/mpiicc.emsolve
      MPICXX=$SRCDIR/helpers/mpiicpc.emsolve
      PYTHONCC=icc
      PYTHONCXX=icpc
      GCCXMLCC=gcc-3.2.1
      GCCXMLCXX=g++-3.2.1
      CMAKECC=gcc-3.2.1
      CMAKECXX=g++-3.2.1
      PARMETISCC=$MPICC
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

## On 64 bit Darwin we have to diddle the python configure 
#if test $OSNAME = "Darwin"; then
#   PYTHONCONFFLAGS="--enable-framework=$SPHERALDIR DESTDIR=$SPHERALDIR"
#   #PYTHONCONFFLAGS="'MACOSX_DEPLOYMENT_TARGET=10.5' --enable-framework=$SPHERALDIR --enable-universalsdk" # --disable-toolbox-glue"
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
   PYTHONCC=$CC
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
   PYTHONCXX=$CXX
   AC_MSG_RESULT($PYTHONCXX)
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
   GCCXMLCC=$CC
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
   GCCXMLCXX=$CXX
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
   CMAKECC=$CC
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
   CMAKECXX=$CXX
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
   PARMETISCC=$MPICC
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
if test $CXXCOMPILERTYPE = "GNU" -o $CXXCOMPILERTYPE = "INTEL"; then
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
  LDRPATH="$LDRPATH ${LDPASSTHROUGH}-rpath=$TOPLIBDIR ${LDPASSTHROUGH}-rpath=$LIBDIR"

elif test "$OSNAME" = "Darwin"; then
  LDRPATH="$LDRPATH ${LDPASSTHROUGH}-rpath $TOPLIBDIR ${LDPASSTHROUGH}-rpath $LIBDIR"
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
  FORTLINK=" -lgfortran"
  if test "$OSNAME" = "Darwin"; then
     CFLAGS="$CFLAGS -fPIC"
     CXXFLAGS="$CXXFLAGS -fPIC -DHAVE_XCPT -DGNUCXX"
     FORTFLAGS="$FORTFLAGS -fPIC"
     DYNLIBFLAG="-dynamiclib -undefined suppress -flat_namespace"
     SHAREDFLAG="-bundle -undefined suppress -flat_namespace"
     JAMTOOLSET=darwin
     LDFLAGS="$LDFLAGS -undefined dynamic_lookup"
  else
     CFLAGS="$CFLAGS -fpic -fexceptions"
     CXXFLAGS="$CXXFLAGS -fpic -fexceptions -DHAVE_XCPT -DGNUCXX"
     FORTFLAGS="$FORTFLAGS -fpic"
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
  CFLAGS="$CFLAGS -fpic -wd654"
  CXXFLAGS="$CXXFLAGS -fpic -wd654"
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
  CXXFLAGS="$CXXFLAGS -I/usr/include -DHAVE_XCPT -qstaticinline -qtempinc -qrtti=dynamiccast"
  FORTFLAGS="$FORTFLAGS -fpic"
  SHAREDFLAG="$SHAREDFLAG -G -qmkshrobj"
  DEPFLAG="-M -E"
  DEPENDRULES="dependrules.aix"
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
#   JAMTOOLSETOPTS="$JAMTOOLSETOPTS : : $SPHERALDIR/Python.framework "
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

])
