# -----------------------------------------------------------------
# matplotlib -- a nifty plotting package for python.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_MATPLOTLIB],[
AC_SUBST(EXTRATHIRDPARTYTARGETS)

AC_MSG_CHECKING(for --with-matplotlib)
AC_ARG_WITH(matplotlib,
[  --with-matplotlib ........................ build the matplotlib python graphics package],
[
    AC_MSG_RESULT(yes)
    EXTRATHIRDPARTYTARGETS+=" matplotlib_pip_install.date"
],
[
    AC_MSG_RESULT(no)
])

])

