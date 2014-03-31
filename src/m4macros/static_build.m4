# -----------------------------------------------------------------
# static_build.  Optionally build the CXXLIBS as static library
# targets.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_STATIC_BUILD],[
AC_SUBST(DYLIBEXT)
AC_SUBST(LD)
AC_SUBST(LDINSTALLNAME)

AC_MSG_CHECKING(for --with-static-libs)
AC_ARG_WITH(static-libs,
[  --with-static-libs ....................... link C++ libs statically (only works for cxxonly builds)],
[
    AC_MSG_RESULT(yes)
    DYLIBEXT=a 
],
[
    AC_MSG_RESULT(no)
])

])

