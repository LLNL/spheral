# -----------------------------------------------------------------
# Optional GEODYN interface
# -----------------------------------------------------------------
AC_DEFUN([SETUP_GEODYN],[
AC_SUBST(GEODYNSRCS)
AC_SUBST(GEODYNLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(LDFLAGS)

# -----------------------------------------------------------------
# Optionally build the GEODYN package
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-geodyn)
AC_ARG_WITH(geodyn,
[  --with-geodyn ............................ optionally build the interface to GEODYN (requires the external GEODYN library)],
[
   AC_MSG_RESULT(yes)
   GEODYNSRCS="GeodynInst.cc.py"
   PYTHONPKGS+=" Geodyn"
],
[
   AC_MSG_RESULT(no)
   GEODYNSRCS=""
]
)

# -----------------------------------------------------------------
# Optionally override the link arguments for GEODYN
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-geodyn-link)
AC_ARG_WITH(geodyn-link,
[  --with-geodyn-link=ARG ................... change how to link with the external GEODYN library],
[
   AC_MSG_RESULT($withval)
   GEODYNLIBS="$withval"
],
[
   AC_MSG_RESULT(no)
   GEODYNLIBS=""
]
)

])

