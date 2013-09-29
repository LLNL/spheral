# -----------------------------------------------------------------
# matplotlib -- a nifty plotting package for python.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_MATPLOTLIB],[
AC_SUBST(MATPLOTLIBTARGETS)

AC_MSG_CHECKING(for --with-matplotlib)
AC_ARG_WITH(matplotlib,
[  --with-matplotlib ........................ build the matplotlib python graphics package],
[
    AC_MSG_RESULT(yes)
    MATPLOTLIBTARGETS=".matplotlib-1.0.0.date"
],
[
    AC_MSG_RESULT(no)
    MATPLOTLIBTARGETS=""
])

])

