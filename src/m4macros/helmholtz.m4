# -----------------------------------------------------------------
# Optional HELMHOLTZ interface
# -----------------------------------------------------------------
AC_DEFUN([SETUP_HELMHOLTZ],[
AC_SUBST(HELMSRCS)
AC_SUBST(HELMFSRCS)
AC_SUBST(HELMLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(LDFLAGS)

# -----------------------------------------------------------------
# Optionally build the HELMHOLTZ package
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-helmholtz)
AC_ARG_WITH(helmholtz,
[  --with-helmholtz ......................... optionally build the interface to HELMHOLTZ (requires the external HELMHOLTZ library)],
[
   AC_MSG_RESULT(yes)
   HELMSRCS="HelmholtzEquationOfStateInst.cc.py"
   HELMFSRCS="public_helm.f90 invert_helm.f90"
   PYTHONPKGS+=" Helmholtz" #does something go here!??!
   LDFLAGS+=" $FORTLINK"
],
[
   AC_MSG_RESULT(no)
   HELMSRCS=""
   HELMFSRCS=""
]
)

# -----------------------------------------------------------------
# Optionally override the link arguments for HELMHOLTZ
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-helmholtz-link)
AC_ARG_WITH(helmholtz-link,
[  --with-helmholtz-link=ARG ................ change how to link with the external HELMHOLTZ library],
[
   AC_MSG_RESULT($withval)
   HELMLIBS="$withval"
],
[
   AC_MSG_RESULT(no)
   HELMLIBS=""
]
)

])

