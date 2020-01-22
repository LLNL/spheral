# =======================================================================
# Spheral specific autoconf macros
# =======================================================================

# AC_DEFUN([CHOOSE_PKGS],[
# AC_SUBST(CXXPKGS)

# AC_MSG_CHECKING(for --with-geometry-only)
# AC_ARG_WITH(geometry-only,
# [  --with-geometry-only ..................... compile just the Geometry package],
# [
#    AC_MSG_RESULT(yes)
#    CXXPKGS="Geometry Utilities $CXXPKGS"
# ],[
#    AC_MSG_RESULT(no)
#    CXXPKGS="Geometry NodeList Field FieldOperations Kernel Material Neighbor DataBase Boundary ArtificialViscosity Physics Hydro ExternalForce Gravity Integrator FileIO DataOutput Utilities $CXXPKGS"
# ])
# ])
   
AC_DEFUN([SETUP_SPHERAL_ENV],[

AC_SUBST(SPHERALBUILDDIR)
AC_SUBST(HEADERDIR)
AC_SUBST(CXXFLAGS)
AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYFFLEENTRY)
AC_SUBST(LDFLAGS)
AC_SUBST(LDRPATH)
AC_SUBST(LIBS)
AC_SUBST(GEOMETRY_ONLY)
AC_SUBST(ALL)
AC_SUBST(ALLTOP)
AC_SUBST(CXXONLY)
AC_SUBST(USE_R3D)
AC_SUBST(TPINCS)
AC_SUBST(TPLIBS)
AC_SUBST(EXTRATHIRDPARTYTARGETS)
AC_SUBST(PIPTARGETS)
AC_SUBST(CMAKEEXE)
AC_SUBST(BOOSTTARGET)
AC_SUBST(SILOTARGET)

LDRPATH=
HEADERDIR=
TPLIBS=
BOOSTTARGET=
SILOTARGET=

PIPTARGETS+=" PYB11Generator sphinx sphinx_rtd_theme pytest-runner numpy-stl twine h5py cython"

AC_MSG_CHECKING(for spheral build directory)
#SPHERALBUILDDIR=`echo $PWD | sed -e "s/\/spheral\/src$//g;"`
SPHERALBUILDDIR=`echo $PWD`
AC_MSG_RESULT($SPHERALBUILDDIR)

# We default prefix to a system subdirectory of the build tree.
HOST="`uname -s`"
#AC_PREFIX_DEFAULT($SPHERALBUILDDIR)

# Choose the packages we're building.
AC_MSG_CHECKING(for --with-geometry-only)
AC_ARG_WITH(geometry-only,
[  --with-geometry-only ..................... compile just the Geometry package],
[
   AC_MSG_RESULT(yes)
   CXXPKGS="Geometry Utilities $CXXPKGS"
   CXXPKGLIBS="$CXXPKGS"
   GEOMETRY_ONLY=1
],[
   AC_MSG_RESULT(no)
   CXXPKGS="Geometry NodeList Field FieldOperations Kernel Material Neighbor DataBase Boundary Physics ArtificialViscosity Hydro ExternalForce Gravity Integrator FileIO DataOutput Utilities NodeGenerators SPH RK CRKSPH FVCRKH SVPH Mesh Damage SolidMaterial Strength ArtificialConduction SimulationControl $CXXPKGS"
   CXXPKGLIBS="Geometry NodeList Field FieldOperations Kernel Material Neighbor DataBase Boundary Physics ArtificialViscosity Hydro ExternalForce Gravity Integrator FileIO DataOutput Utilities NodeGenerators SPH RK CRKSPH FVCRKH SVPH Mesh Damage SolidMaterial Strength ArtificialConduction $CXXPKGLIBS"
   GEOMETRY_ONLY=0
])

echo "prefix is $prefix"
echo "HOST is $HOST"
echo "SPHERALBUILDDIR is $SPHERALBUILDDIR"
echo "LIBDIR is $LIBDIR"
echo "CXXPKGS is $CXXPKGS"
echo "CXXPKGLIBS is $CXXPKGLIBS"

# =======================================================================
# IBMs are weird about shared objects.  Here's a pile of crap (POC) that
# we have to do to get things running on the IBMs.
# =======================================================================
AC_SUBST(IMPMODS)
AC_SUBST(AIXLIBS)
AC_SUBST(PYFFLEENTRY)
AC_SUBST(MAKEIMPORTFILE)
AC_SUBST(CHECKLIBS)
AC_SUBST(DEPENDRULES)
AC_SUBST(AIXSHELL)
AC_SUBST(CONFIG_SHELL)
IMPMODS=""
AIXLIBS=""
PYFFLEENTRY=""
MAKEIMPORTFILE="$srcdir/helpers/generateDummyImportFile"
CHECKLIBS="$srcdir/helpers/checkLibsForUndefined"
DEPENDRULES="dependrules.generic"
AIXSHELL=""
CONFIG_SHELL=$SHELL
AC_MSG_CHECKING(python.exp required for linking)

LIBS=
if test "`uname -s`" = "AIX"; then
  #IMPMODS="$CXXPKGS"
  #PYFFLEENTRY="-e initlibPyffle"
  #MAKEIMPORTFILE=$srcdir/helpers/generateImportFile"

  # 32 bit
  #AIXLIBS="/lib/crt0.o -lm"
  #LDFLAGS="$LDFLAGS -Wl,-bbigtoc -Wl,-brtl -Wl,-bdynamic"

  # 64 bit
  AIXLIBS="/lib/crt0_64.o -lm"
  LDFLAGS="$LDFLAGS -L/lib -Wl,-bbigtoc -Wl,-bexpall -Wl,-brtl -Wl,-bdynamic -Wl,-b64"

  # This is nuts, but building gcc with the stock AIX sh takes a *day*, so force
  # the use of the bash on this platform.  HACK!
  CONFIG_SHELL="/usr/local/bin/bash"
  AIXSHELL="SHELL=/usr/local/bin/bash CONFIG_SHELL=/usr/local/bin/bash"

fi
echo "SHAREDFLAG is $SHAREDFLAG"

# -----------------------------------------------------------------
# We must be on a 32 bit intel processor in order to use the psyco 
# python accelerator.
# -----------------------------------------------------------------
# if (test -n "`uname -a | grep i386`" -o -n "`uname -a | grep i486`" \
#       -o -n "`uname -a | grep i586`" -o -n "`uname -a | grep i686`"); then
#   EXTRATHIRDPARTYTARGETS+=" .psyco-1.3-src.date"
# fi

# -----------------------------------------------------------------
# If we're on AIX, set a few special third party lib options.
# -----------------------------------------------------------------
AC_SUBST(GCCXMLDIST)
if test "`uname -s`" = "AIX"; then
  GCCXMLDIST="gccxml-cvssnapshot-2008-02-04.tar.bz2"
  #GCCXMLDIST="gccxml-0.6.0.tar.bz2"
else
  GCCXMLDIST="gccxml-cvssnapshot-2008-02-04.tar.bz2"
fi

# -----------------------------------------------------------------
# On Darwin we have to tell silo to link with readline
# We also need to pass Boost a few specialized flags.
# -----------------------------------------------------------------
AC_SUBST(SILOFLAGS)
SILOFLAGS=""
if test "`uname -s`" = "Darwin"; then
  SILOFLAGS="LDFLAGS=-lreadline"
  CXXFLAGS="$CXXFLAGS -DBOOST_DATE_TIME_NO_LOCALE"
fi

# -----------------------------------------------------------------
# Select the script to use building the WildMagic third party 
# target.
# -----------------------------------------------------------------
AC_MSG_CHECKING(build WildMagic)
AC_SUBST(BUILDWILDMAGIC)
AC_SUBST(WILDMAGICTARGET)
AC_SUBST(WMLIBEXT)
if test "`uname -s`" = "Darwin"; then
  BUILDWILDMAGIC="MacBuildWm5.csh"
  WILDMAGICTARGET="Release"
  WMLIBEXT="a"
else
  BUILDWILDMAGIC="WildMagic5p4_make.csh"
  WILDMAGICTARGET="ReleaseDynamic"
  WMLIBEXT="so"
fi
AC_MSG_RESULT($BUILDWILDMAGIC)

# -----------------------------------------------------------------
# Optionally build an additional package of C++ testing functions.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for cxxtests)
AC_ARG_WITH(cxxtests,
[  --with-cxxtests .......................... optionally build the C++ testing methods],
[
   AC_MSG_RESULT(yes)
   CXXPKGS+=" CXXTests"
   CXXPKGLIBS+=" CXXTests"
   PYTHONPKGS+=" CXXTests"
],
[
   AC_MSG_RESULT(no)
]
)

# -----------------------------------------------------------------
# Optionally build install the GSL (Gnu Scientific Library) 
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-gsl)
AC_ARG_WITH(gsl,
[  --with-gsl ............................... optionally install the Gnu Scientific Library extensions],
[
   AC_MSG_RESULT(yes)
   EXTRATHIRDPARTYTARGETS+=" .gsl-1.14.date .pygsl-0.9.5.date"
],
[
   AC_MSG_RESULT(no)
]
)

# -----------------------------------------------------------------
# Optionally build install mpmath
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-mpmath)
AC_ARG_WITH(mpmath,
[  --with-mpmath ............................ optionally install the Gnu Scientific Library extensions],
[
   AC_MSG_RESULT(yes)
   PIPTARGETS+=" mpmath"
],
[
   AC_MSG_RESULT(no)
]
)

# -----------------------------------------------------------------
# Optionally do not build third party libs.
# -----------------------------------------------------------------
AC_SUBST(BUILDTHIRDPARTYTARGET)
AC_MSG_CHECKING(for --without-thirdPartyLibs)
AC_ARG_WITH(thirdPartyLibs,
[  --without-thirdPartyLibs ................. do not build the third party libraries],
[
    AC_MSG_RESULT(yes)
    BUILDTHIRDPARTYTARGET=""
],
[
    AC_MSG_RESULT(no)
    BUILDTHIRDPARTYTARGET="thirdPartyLibs"
])

# -----------------------------------------------------------------
# Optionally do not build numpy.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-numpy)
AC_ARG_WITH(numpy,
[  --without-numpy .......................... do not build the NumPy third party extension],
[
    AC_MSG_RESULT(yes)
],
[
    AC_MSG_RESULT(no)
    PIPTARGETS+=" numpy"
])

# -----------------------------------------------------------------
# Optionally do not build gnuplot python interface.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-gnuplot-py)
AC_ARG_WITH(numpy,
[  --without-gnuplot-py ..................... do not build the Gnuplot python extension],
[
    AC_MSG_RESULT(yes)
],
[
    AC_MSG_RESULT(no)
    EXTRATHIRDPARTYTARGETS+=" .gnuplot-py-1.8.date"
])

# -----------------------------------------------------------------
# Optionally do not install the python sobol package.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-sobol)
AC_ARG_WITH(sobol,
[  --without-sobol .......................... do not build the Sobol python third party extension],
[
    AC_MSG_RESULT(yes)
],
[
    AC_MSG_RESULT(no)
    PIPTARGETS+=" sobol"
])

# -----------------------------------------------------------------
# Allow the use of an existing cmake.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-cmake)
AC_ARG_WITH(cmake,
[  --with-cmake ............................. specify a cmake executable],
[
    CMAKEEXE=$withval
    AC_MSG_RESULT($CMAKEEXE)
],
[
    CMAKEEXE='$(prefix)/bin/cmake'
    AC_MSG_RESULT($CMAKEEXE)
])

# -----------------------------------------------------------------
# Allow the use of an existing hdf5.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-hdf5)
AC_ARG_WITH(hdf5,
[  --with-hdf5 .............................. specify a root path for an existing hdf5 install],
[
    TPINCS+=" -I $withval/include"
    TPLIBS+=" -L $withval/lib"
    AC_MSG_RESULT($withval)
],
[
    EXTRATHIRDPARTYTARGETS+=" .hdf5-1.10.4.date"
    AC_MSG_RESULT(no)
])

# -----------------------------------------------------------------
# Allow the use of an existing silo.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-silo)
AC_ARG_WITH(silo,
[  --with-silo .............................. specify a root path for an existing silo install],
[
    TPINCS+=" -I $withval/include"
    TPLIBS+=" -L $withval/lib"
    AC_MSG_RESULT($withval)
],
[
    EXTRATHIRDPARTYTARGETS+=" .silo-4.10.2-bsd.date"
    SILOTARGET=".silo-4.10.2-bsd.date"
    AC_MSG_RESULT(no)
])

# -----------------------------------------------------------------
# Allow the use of an existing boost.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-boost)
AC_ARG_WITH(boost,
[  --with-boost ............................. specify a root path for an existing boost install],
[
    TPINCS+=" -I $withval/include"
    AC_MSG_RESULT($withval)
],
[
    EXTRATHIRDPARTYTARGETS+=" .boost_1_69_0.date"
    BOOSTTARGET=".boost_1_69_0.date"
    AC_MSG_RESULT(no)
])

# -----------------------------------------------------------------
# Allow the use of an existing qhull.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-qhull)
AC_ARG_WITH(qhull,
[  --with-qhull ............................. specify a root path for an existing qhull install],
[
    TPINCS+=" -I $withval/include"
    TPLIBS+=" -L $withval/lib"
    AC_MSG_RESULT($withval)
],
[
    EXTRATHIRDPARTYTARGETS+=" .qhull-2015-src-7.2.0.date"
    AC_MSG_RESULT(no)
])

# -----------------------------------------------------------------
# Allow the use of an existing eigen.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-eigen)
AC_ARG_WITH(eigen,
[  --with-eigen ............................. specify a root path for an existing eigen install],
[
    TPINCS+=" -I $withval/include"
    AC_MSG_RESULT($withval)
],
[
    EXTRATHIRDPARTYTARGETS+=" .eigen.date"
    AC_MSG_RESULT(no)
])

# -----------------------------------------------------------------
# Optionally do a cxxonly build.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-cxxonly)
AC_ARG_WITH(cxxonly,
[  --with-cxxonly ........................... optionally do a cxxonly build (no third party, no python extensions)],
[
   AC_MSG_RESULT(yes)
   CXXONLY="yes"
   ALL="\$(LIBTARGET) \$(STATICLIBTARGET) \$(EXETARGETS) \$(THIRDPARTYLIBTARGET) \$(INCTARGETS) \$(INSTALLTARGETS) install_headers"
   ALLTOP="\$(CXXPKGS) install_headers"
   CXXFLAGS+=" -DCXXONLY"
   LIBDIR="$libdir"
],
[
   AC_MSG_RESULT(no)
   CXXONLY="no"
   ALL="\$(PYTHONTARGETS) \$(LIBTARGET) \$(STATICLIBTARGET) \$(MODTARGET) \$(BPLMODTARGET) \$(PBGMODTARGET) \$(PYB11TARGET) \$(STATICPBGMODTARGET) \$(EXETARGETS) \$(THIRDPARTYLIBTARGET) \$(INCTARGETS) \$(INSTALLTARGETS)"
   ALLTOP="$PYTHONPKGDIR \$(CXXPKGS) pycompileall"
]
)

])


