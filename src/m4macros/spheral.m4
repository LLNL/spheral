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

AC_SUBST(SPHERALDIR)
AC_SUBST(SRCDIR)
AC_SUBST(TOPLIBDIR)
AC_SUBST(LIBDIR)
AC_SUBST(CXXFLAGS)
AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYFFLEENTRY)
AC_SUBST(LDFLAGS)
AC_SUBST(LDRPATH)
AC_SUBST(LIBS)
AC_SUBST(GEOMETRY_ONLY)

LDRPATH=

# Prepare for user selected third party targets.
AC_SUBST(EXTRATHIRDPARTYTARGETS)
EXTRATHIRDPARTYTARGETS=""

AC_MSG_CHECKING(for spheral directory)
SRCDIR=`echo $PWD`
SPHERALDIR=`echo $PWD | sed -e "s/\/spheral\/src$//g;"`
AC_MSG_RESULT($SPHERALDIR)

TOPLIBDIR=$PYTHONROOT/lib
LIBDIR=$PYTHONSITEPKGDIR/Spheral

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
   CXXPKGS="Geometry NodeList Field FieldOperations Kernel Material Neighbor DataBase Boundary ArtificialViscosity Physics Hydro ExternalForce Gravity Integrator FileIO DataOutput Utilities NodeGenerators SimulationControl SPH CSPH Mesh Damage SolidMaterial SolidSPH Strength $CXXPKGS"
   CXXPKGLIBS="Geometry NodeList Field FieldOperations Kernel Material Neighbor DataBase Boundary ArtificialViscosity Physics Hydro ExternalForce Gravity Integrator FileIO DataOutput Utilities NodeGenerators SPH CSPH Mesh Damage SolidMaterial SolidSPH Strength $CXXPKGLIBS"
   GEOMETRY_ONLY=0
])

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
MAKEIMPORTFILE="$SRCDIR/helpers/generateDummyImportFile"
CHECKLIBS="$SRCDIR/helpers/checkLibsForUndefined"
DEPENDRULES="dependrules.generic"
AIXSHELL=""
CONFIG_SHELL=$SHELL
AC_MSG_CHECKING(python.exp required for linking)

if test "`uname -s`" = "AIX"; then
  #IMPMODS="$CXXPKGS"
  #PYFFLEENTRY="-e initlibPyffle"
  #MAKEIMPORTFILE=$SRCDIR/helpers/generateImportFile"
  LIBS=

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
  EXTRATHIRDPARTYTARGETS+=" .numpy-1.6.2.date .gnuplot-py-1.8.date"
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

# -----------------------------------------------------------------
# Optionally build an additional package of C++ testing functions.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for cxxtests)
AC_ARG_WITH(cxxtests,
[  --with-cxxtests .......................... optionally build the C++ testing methods],
[
   AC_MSG_RESULT(yes)
   CXXPKGS="$CXXPKGS CXXTests"
   CXXPKGLIBS="$CXXPKGLIBS CXXTests"
   PYTHONPKGS="$PYTHONPKGS CXXTests"
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

])


