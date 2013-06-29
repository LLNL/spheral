# -----------------------------------------------------------------
# scipy -- a collection of math & science oriented python extensions.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_SCIPY],[
AC_SUBST(SCIPYTARGETS)

AC_MSG_CHECKING(for --with-scipy)
AC_ARG_WITH(scipy,
[  --with-scipy ............................. build the scipy python package],
[
    AC_MSG_RESULT(yes)
    SCIPYTARGETS=".scipy-0.10.1.date"
],
[
    AC_MSG_RESULT(no)
    SCIPYTARGETS=""
])

])

