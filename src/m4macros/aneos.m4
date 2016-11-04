# -----------------------------------------------------------------
# Optional ANEOS interface
# -----------------------------------------------------------------
AC_DEFUN([SETUP_ANEOS],[
AC_SUBST(ANEOSSRCS)
AC_SUBST(ANEOSFSRCS)
AC_SUBST(ANEOSLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(LDFLAGS)

# -----------------------------------------------------------------
# Optionally build the ANEOS package
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-aneos)
AC_ARG_WITH(aneos,
[  --with-aneos ............................. optionally build the interface to ANEOS (requires the external ANEOS library)],
[
   AC_MSG_RESULT(yes)
   ANEOSSRCS="ANEOSInst.cc.py"
   ANEOSFSRCS="ANEOS_initialize.f"
   PYTHONPKGS+=" ANEOS"
   LDFLAGS+=" $FORTLINK"
],
[
   AC_MSG_RESULT(no)
   ANEOSSRCS=""
   ANEOSFSRCS=""
]
)

# -----------------------------------------------------------------
# Optionally override the link arguments for ANEOS
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-aneos-link)
AC_ARG_WITH(aneos-link,
[  --with-aneos-link=ARG .................... change how to link with the external ANEOS library],
[
   AC_MSG_RESULT($withval)
   ANEOSLIBS="$withval"
],
[
   AC_MSG_RESULT(no)
   ANEOSLIBS=""
]
)

])

