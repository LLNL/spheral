dnl -----------------------------------------------------------------
dnl dbc.m4 -- support for Design By Contract macros
dnl -----------------------------------------------------------------

AC_DEFUN([SETUP_DBC],
[
  AC_MSG_CHECKING(for --with-dbc)
  AC_ARG_WITH(dbc,
  [  --with-dbc=MODE .......................... choose DBC mode [all,pre,none]],
  [
    AC_MSG_RESULT($withval)
    if test $withval = "all"; then
      CXXFLAGS="$CXXFLAGS -DDBC_COMPILE_ALL"
    fi
    if test $withval = "pre"; then
      CXXFLAGS="$CXXFLAGS -DDBC_COMPILE_PRE"
    fi
  ],
  [
    AC_MSG_RESULT(none, defaulting to no contracts)
    CXXFLAGS+=" -DNDEBUG"
  ])
])

