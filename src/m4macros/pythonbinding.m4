# -----------------------------------------------------------------
# Determine what we're using to bind to python.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_PYTHONBINDING],[

AC_SUBST(PYTHONBINDING)
AC_SUBST(PYTHONPKGDIR)
AC_SUBST(PYTHONPKGS)
AC_SUBST(BOOSTROOT)
AC_SUBST(BOOSTLIBTARGETS)
AC_SUBST(BOOSTPYTHONTARGET)
AC_SUBST(GCCXMLTARGETS)
AC_SUBST(BPLPATH)
AC_SUBST(PYSTEPATH)
AC_SUBST(BPLINCS)
AC_SUBST(PYOPT)
AC_SUBST(MODULELINK)

PYOPT=""
BPLINCS=""
BOOSTLIBTARGETS="math"
BOOSTPYTHONTARGET=""

# -----------------------------------------------------------------
# Configure for using Boost.Python library
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-bpl)
AC_ARG_WITH(bpl,
[  --with-bpl ............................... use Boost.Python wrappings],
[
    AC_MSG_RESULT(yes)
    BOOSTPYTHONTARGET="python"
    GCCXMLTARGETS="$(GCCXMLDATE)"
    PYTHONBINDING="BPL"
    PYTHONPKGDIR="BPLWraps"
    if test "$GEOMETRY_ONLY" = "1"; then
      PYTHONPKGS="$PYTHONPKGS Geometry Utilities"
    else
      PYTHONPKGS="$PYTHONPKGS Geometry DataOutput NodeList NodeIterators Field FieldOperations Kernel SplineKernel Neighbor Material FileIO DataBase Boundary ArtificialViscosity Physics Hydro ExternalForce Gravity Integrator CXXTypes Utilities Python NodeGenerators"
    fi
    PYOPT="$PYOPT -w"
    BPLINCS="-DPYSTE -I\$(BOOSTROOT) -I\$prefix/include/python\$(PYTHONVERSION) \$(patsubst %, -I\$(SRCTOP)/%, \$(CXXPKGS))"
    MODULELINK="-L\$(LIBDIR) -lboost_python \$(PKGLIBS)"
    if test "`uname -s`" = "AIX"; then
       MODULELINK="$MODULELINK -e init\$(PKGNAME)"
       BPLINCS="$BPLINCS -D_WCHAR_T"
    fi
],
[
    AC_MSG_RESULT(no)
    GCCXMLTARGETS=""
])

AC_MSG_CHECKING(for --without-pybindgen)
AC_ARG_WITH(pybindgen,
[  --without-pybindgen ...................... do not use pybindgen wrappings],
[
    AC_MSG_RESULT(yes)
],
[
    AC_MSG_RESULT(no)
    PYTHONBINDING="PYBINDGEN"
    PYTHONPKGDIR="PBGWraps"
    PYTHONPKGS="Geometry CXXTypes Silo DataOutput NodeList Field Kernel Neighbor Material FileIO DataBase Boundary Physics ArtificialViscosity Hydro ExternalForce Gravity Integrator Utilities NodeGenerators FieldOperations SPH CRKSPH SVPH TaylorSPH Mesh Damage SolidMaterial SolidSPH Strength ArtificialConduction $PYTHONPKGS"
    PYOPT="$PYOPT"
    MODULELINK="-L\$(LIBDIR) \$(PKGLIBS)"
    if test "`uname -s`" = "AIX"; then
       MODULELINK="$MODULELINK -e init\$(PKGNAME)"
    fi
])

AC_MSG_CHECKING(for --with-boostroot)
AC_ARG_WITH(boostroot,
  [  --with-boostroot=ARG ..................... set the path to boost],
  [
    AC_MSG_RESULT($withval)
    BOOSTROOT="$withval"
  ],
  [
    AC_MSG_RESULT(none)
    BOOSTROOT="\$(prefix)/include/boost"
  ]
)

BPLPATH="$LIBDIR"

echo "PYTHONBINDING is $PYTHONBINDING"
echo "PYTHONPKGDIR is $PYTHONPKGDIR"
echo "PYTHONPKGS is $PYTHONPKGS"
echo "PYOPT is $PYOPT"
echo "BOOSTROOT is $BOOSTROOT"
echo "BPLPATH is $BPLPATH"
echo "PYSTEPATH is $PYSTEPATH"

])
