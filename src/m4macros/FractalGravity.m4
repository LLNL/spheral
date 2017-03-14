# ----------------------------------------------------------------
# Configuration for FractalGravity package.
# ----------------------------------------------------------------
AC_DEFUN([SETUP_FRACTALGRAVITY],[

AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYTHONPKGS)

AC_MSG_CHECKING(for FractalGravity)
AC_ARG_WITH(FractalGravity,
[  --with-FractalGravity .................... optionally build the FractalGravity package],
[
   AC_MSG_RESULT(yes)
   CXXPKGS="$CXXPKGS FractalGravity"
   CXXPKGLIBS="$CXXPKGLIBS FractalGravity"
   PYTHONPKGS="$PYTHONPKGS FractalGravity"
],
[
   AC_MSG_RESULT(no)
]
)

])
