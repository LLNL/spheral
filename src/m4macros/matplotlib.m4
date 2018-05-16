# -----------------------------------------------------------------
# matplotlib -- a nifty plotting package for python.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_MATPLOTLIB],[
AC_SUBST(EXTRATHIRDPARTYTARGETS)
AC_SUBST(PIPTARGETS)

AC_MSG_CHECKING(for --with-matplotlib)
AC_ARG_WITH(matplotlib,
[  --without-matplotlib ..................... do not build the matplotlib python graphics package],
[
    AC_MSG_RESULT(yes)
],
[
    AC_MSG_RESULT(no)
    PIPTARGETS+=" matplotlib"
])

])

