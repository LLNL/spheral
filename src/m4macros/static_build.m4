# -----------------------------------------------------------------
# static_build.  Optionally build the CXXLIBS as static library
# targets.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_STATIC_BUILD],[
AC_SUBST(LIBTOOL)
AC_SUBST(LIBTOOLLINKFLAGS)
AC_SUBST(DYLIBEXT)

LIBTOOL="libtool"
LIBTOOLLINKFLAGS="--mode=link --tag=CXX"

AC_MSG_CHECKING(for --with-static-libs)
AC_ARG_WITH(static-libs,
[  --with-static-libs ....................... link C++ libs statically (only works for cxxonly builds)],
[
    AC_MSG_RESULT(yes)
    DYLIBEXT=la
],
[
    AC_MSG_RESULT(no)
])

])

