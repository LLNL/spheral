# ----------------------------------------------------------------
# Configuration for MHD package.
# ----------------------------------------------------------------
AC_DEFUN([SETUP_MHD],[

AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(USEPETSC)

AC_MSG_CHECKING(for MHD)
AC_ARG_WITH(MHD,
[  --with-MHD ............................... optionally build the MHD physics package],
[
   AC_MSG_RESULT(yes)
   CXXPKGS="$CXXPKGS MHD"
   CXXPKGLIBS="$CXXPKGLIBS MHD"
   PYTHONPKGS="$PYTHONPKGS MHD"
   USEPETSC="1"
],
[
   AC_MSG_RESULT(no)
]
)

])
