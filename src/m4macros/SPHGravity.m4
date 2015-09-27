# ----------------------------------------------------------------
# Configuration for the SPH gravity package.
# ----------------------------------------------------------------
AC_DEFUN([SETUP_SPHGRAVITY],[

AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(USEPETSC)

AC_MSG_CHECKING(for SPHGravity)
AC_ARG_WITH(SPHGravity,
[  --with-SPHGravity ........................ optionally build the SPH Gravity physics package],
[
   AC_MSG_RESULT(yes)
   CXXPKGS="$CXXPKGS SPHGravity"
   CXXPKGLIBS="$CXXPKGLIBS SPHGravity"
   PYTHONPKGS="$PYTHONPKGS SPHGravity"
   USEPETSC="1"
],
[
   AC_MSG_RESULT(no)
]
)

])
