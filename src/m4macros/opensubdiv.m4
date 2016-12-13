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
    OPENSUBDIVTARGETS=".OpenSubdiv-2_5_0.date"
    #OPENSUBDIVLIBS="\$(prefix)/lib/libosdCPU.a \$(prefix)/lib/libosdutil.a"
    if test "`uname -s`" = "Darwin"; then
        OPENSUBDIVLIBS="\$(prefix)/lib/libosdCPU.a \$(prefix)/lib/libosdutil.a"
    else
        OPENSUBDIVLIBS="-losdCPU -losdutil"
    fi
    CXXFLAGS="$CXXFLAGS -DHAVE_OPENSUBDIV"
])

])

