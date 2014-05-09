# -----------------------------------------------------------------
# opensubdiv -- a collection of math & science oriented python extensions.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_OPENSUBDIV],[
AC_SUBST(OPENSUBDIVTARGETS)
AC_SUBST(OPENSUBDIVLIBS)
AC_SUBST(CXXFLAGS)

AC_MSG_CHECKING(for --without-opensubdiv)
AC_ARG_WITH(opensubdiv,
[  --without-opensubdiv ..................... do not build the opensubdiv package],
[
    AC_MSG_RESULT(yes)
    OPENSUBDIVTARGETS=""
    OPENSUBDIVLIBS=""
],
[
    AC_MSG_RESULT(no)
    OPENSUBDIVTARGETS=".OpenSubdiv-master.date"
    OPENSUBDIVLIBS="-losdCPU -losdutil"
    CXXFLAGS="$CXXFLAGS -DHAVE_OPENSUBDIV"
])

])

