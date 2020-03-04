AC_DEFUN([SETUP_DEBUG_PRINT],[
AC_SUBST(CXXFLAGS)
AC_SUBST(PYOPT)

AC_MSG_CHECKING(for --with-debug)
AC_ARG_WITH(debug, 
[  --with-debug ............................. turn on debug printing],
[
    AC_MSG_RESULT(yes)
    CXXFLAGS="$CXXFLAGS -DDEBUG"
    echo "Debug printing is on."
],
[
    AC_MSG_RESULT(no)
])

AC_MSG_CHECKING(for --with-boundchecking)
AC_ARG_WITH(boundchecking, 
[  --with-boundchecking ..................... turn on debug printing],
[
    AC_MSG_RESULT(yes)
    CXXFLAGS="$CXXFLAGS -D_GLIBCXX_DEBUG"
],
[
    AC_MSG_RESULT(no)
])
])
