# ----------------------------------------------------------------
# Configuration for FractalGravity package.
# ----------------------------------------------------------------
AC_DEFUN([SETUP_FRACTALGRAVITY],[

AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(EXTRATHIRDPARTYTARGETS)
AC_SUBST(TPINCS)
AC_SUBST(TPLIBS)

AC_MSG_CHECKING(for FractalGravity)
AC_ARG_WITH(FractalGravity,
[  --with-FractalGravity .................... optionally build the FractalGravity package],
[
   AC_MSG_RESULT(yes)
   CXXPKGS+=" FractalStruct"
   CXXPKGLIBS+=" FractalStruct"
   PYTHONPKGS+=" FractalGravity"
   EXTRATHIRDPARTYTARGETS+=" .fftw-3.3.7.date .hypre-v2.14.0.date"
   TPINCS+=" -I\$(prefix)/HYPRE/include"
   TPLIBS+=" -lfftw3 -L\$(prefix)/HYPRE/lib -lHYPRE"
],
[
   AC_MSG_RESULT(no)
]
)

])
