# ----------------------------------------------------------------
# Configuration for FractalGravity package.
# ----------------------------------------------------------------
AC_DEFUN([SETUP_FRACTALGRAVITY],[

AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(EXTRATHIRDPARTYTARGETS)

AC_MSG_CHECKING(for FractalGravity)
AC_ARG_WITH(FractalGravity,
[  --with-FractalGravity .................... optionally build the FractalGravity package],
[
   AC_MSG_RESULT(yes)
   CXXPKGS="$CXXPKGS FractalStruct"
   CXXPKGLIBS="$CXXPKGLIBS FractalStruct"
   PYTHONPKGS="$PYTHONPKGS FractalStruct"
   EXTRATHIRDPARTYTARGETS+=" .fftw-3.3.7.date"
],
[
   AC_MSG_RESULT(no)
]
)

])
