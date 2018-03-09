# -----------------------------------------------------------------
# scipy -- a collection of math & science oriented python extensions.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_SCIPY],[
AC_SUBST(PIPTARGETS)

AC_MSG_CHECKING(for --without-scipy)
AC_ARG_WITH(scipy,
[  --without-scipy .......................... do not build the scipy python package],
[
    AC_MSG_RESULT(yes)
],
[
    AC_MSG_RESULT(no)
    PIPTARGETS+=" scipy"
])

])

