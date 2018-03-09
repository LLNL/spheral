# -----------------------------------------------------------------
# Determine what we're using to bind to python.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_PYTHONBINDING],[

AC_SUBST(EXTRATHIRDPARTYTARGETS)
AC_SUBST(PYTHONBINDING)
AC_SUBST(PYTHONPKGDIR)
AC_SUBST(PYTHONPKGS)
AC_SUBST(BOOSTROOT)
AC_SUBST(BOOSTLIBTARGETS)
AC_SUBST(INCS)
AC_SUBST(PYOPT)
AC_SUBST(MODULELINK)
AC_SUBST(EXTRATHIRDPARTYTARGETS)
AC_SUBST(PIPTARGETS)

PYOPT=""
BOOSTLIBTARGETS="math"
PIPTARGETS+=" pybindgen==0.17.0"

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
    PYTHONPKGS="Geometry CXXTypes PolyClipper Silo DataOutput NodeList Field Kernel Neighbor Material FileIO DataBase Boundary Physics ArtificialViscosity Hydro ExternalForce Gravity Integrator Utilities NodeGenerators FieldOperations SPH CRKSPH SVPH Mesh Damage SolidMaterial Strength ArtificialConduction $PYTHONPKGS"
    PYOPT="$PYOPT"
    MODULELINK="-L\$(LIBDIR) \$(PKGLIBS)"
    if test "`uname -s`" = "AIX"; then
       MODULELINK="$MODULELINK -e init\$(PKGNAME)"
    fi
    if test "$CXXCOMPILERTYPE" = "INTEL"; then
       PYOPT=" -O0 -no-ipo"
    fi
])

# -----------------------------------------------------------------
# Configure using pybind11 library for python bindings
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-pybind11)
AC_ARG_WITH(pybind11,
[  --with-pybind11 .......................... use pybind11 wrappings],
[
    AC_MSG_RESULT(yes)
    PYTHONBINDING="PYBIND11"
    PYTHONPKGDIR="Pybind11Wraps"
    PYTHONPKGS+=" CXXTypes Geometry Silo DataOutput NodeList Field Kernel Neighbor Material FileIO DataBase Boundary Physics ArtificialViscosity Hydro"
    INCS+="-I\$(prefix)/include -I\$prefix/include/python\$(PYTHONVERSION) \$(patsubst %, -I\$(SRCTOP)/%, \$(CXXPKGS))"
    MODULELINK="-L\$(LIBDIR) \$(PKGLIBS)"
    if test "`uname -s`" = "AIX"; then
       MODULELINK="$MODULELINK -e init\$(PKGNAME)"
    fi
],
[
    AC_MSG_RESULT(no)
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

echo "PYTHONBINDING is $PYTHONBINDING"
echo "PYTHONPKGDIR is $PYTHONPKGDIR"
echo "PYTHONPKGS is $PYTHONPKGS"
echo "PYOPT is $PYOPT"
echo "BOOSTROOT is $BOOSTROOT"

])
