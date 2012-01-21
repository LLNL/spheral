
dnl -----------------------------------------------------------------
dnl gadget.m4 -- support for Gadget, that neato gravity code.
dnl -----------------------------------------------------------------

dnl Specify Gadget's directory.
AC_SUBST(GADGET_PATH)
AC_DEFUN([CHECK_GADGET_PATH],
[
  AC_MSG_CHECKING(for specified Gadget path)
  AC_ARG_WITH(gadget-path,
  [  --with-gadget-path=DIR ................... location of Gadget directory ],
  [
    AC_MSG_RESULT($withval)
    GADGET_PATH="$withval"
  ],
  [
    AC_MSG_RESULT(none)
  ])
])

AC_DEFUN([SETUP_GADGET],
[
  AC_MSG_CHECKING(for --enable-gadget)
  AC_ARG_ENABLE(gadget,
  [  --enable-gadget .......................... turn on support for Gadget ],
  [
    AC_MSG_RESULT(yes)
    CHECK_GADGET_PATH
    CXXFLAGS="$CXXFLAGS -DUSE_GADGET -I$GADGET_PATH/S-Gadget"
    CFLAGS="$CFLAGS -DUSE_GADGET -I$GADGET_PATH/S-Gadget"
    CPPFLAGS="$CPPFLAGS -DUSE_GADGET -I$GADGET_PATH/S-Gadget"
    CXXPKGS="$CXXPKGS Gadget"
  ],
  [
    AC_MSG_RESULT(no)
  ])
])

