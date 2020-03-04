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
PIPTARGETS+=" pybindgen==0.17.0"      # if nothing else, polytope currently requires this

AC_MSG_CHECKING(for --with-pybindgen)
AC_ARG_WITH(pybindgen,
[  --with-pybindgen ......................... use pybindgen wrappings],
[
    AC_MSG_RESULT(no)
    PYTHONBINDING="PYBINDGEN"
    PYTHONPKGDIR="PBGWraps"
    PYTHONPKGS="Geometry CXXTypes PolyClipper Silo DataOutput NodeList Field Kernel Neighbor Material FileIO RK DataBase Boundary Physics ArtificialViscosity Hydro ExternalForce Gravity Integrator Utilities NodeGenerators FieldOperations SPH CRKSPH Mesh Damage SolidMaterial Strength ArtificialConduction $PYTHONPKGS"
    PYOPT="$PYOPT"
    MODULELINK="-L\$(LIBDIR) \$(PKGLIBS)"
    if test "`uname -s`" = "AIX"; then
       MODULELINK="$MODULELINK -e init\$(PKGNAME)"
    fi
    if test "$CXXCOMPILERTYPE" = "INTEL"; then
       PYOPT=" -O0 -no-ipo"
    fi
],
[
    AC_MSG_RESULT(yes)
])

# -----------------------------------------------------------------
# Configure using pybind11 library for python bindings
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-pybind11)
AC_ARG_WITH(pybind11,
[  --without-pybind11 ....................... do not use pybind11 wrappings],
[
    AC_MSG_RESULT(no)
],
[
    AC_MSG_RESULT(yes)
    PYTHONBINDING="PYBIND11"
    PYTHONPKGDIR="Pybind11Wraps"
    PYTHONPKGS+=" CXXTypes Geometry PolyClipper Silo DataOutput NodeList Field FieldList Kernel Neighbor Material FileIO Utilities RK DataBase Boundary Physics Hydro ExternalForce Gravity Integrator NodeGenerators FieldOperations SPH CRKSPH SVPH ArtificialViscosity Mesh Damage SolidMaterial Strength ArtificialConduction"
    INCS+="-I\$(prefix)/include -I\$prefix/include/python\$(PYTHONVERSION) \$(patsubst %, -I\$(SRCTOP)/%, \$(CXXPKGS))"
    MODULELINK="-L\$(LIBDIR) \$(PKGLIBS)"
    PIPTARGETS+=" decorator"
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

echo "PYTHONBINDING is $PYTHONBINDING"
echo "PYTHONPKGDIR is $PYTHONPKGDIR"
echo "PYTHONPKGS is $PYTHONPKGS"
echo "PYOPT is $PYOPT"
echo "BOOSTROOT is $BOOSTROOT"

])
