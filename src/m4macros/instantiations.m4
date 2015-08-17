# -----------------------------------------------------------------
# Choose which dimensionalities to instantiate.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_INSTANTIATIONS],[
AC_SUBST(INST1D)
AC_SUBST(INST2D)
AC_SUBST(INST3D)
AC_SUBST(DIMS)

# -----------------------------------------------------------------
# 1D
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-1d)
AC_ARG_WITH(1d,
[  --without-1d ............................. turn off instantiation of the 1D algorithms],
[
   AC_MSG_RESULT(yes)
   INST1D="no"
],
[
   AC_MSG_RESULT(no)
   INST1D="yes"
   DIMS+="1 "
]
)

# -----------------------------------------------------------------
# 2D
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-2d)
AC_ARG_WITH(2d,
[  --without-2d ............................. turn off instantiation of the 2D algorithms],
[
   AC_MSG_RESULT(yes)
   INST2D="no"
],
[
   AC_MSG_RESULT(no)
   INST2D="yes"
   DIMS+="2 "
]
)

# -----------------------------------------------------------------
# 3D
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-3d)
AC_ARG_WITH(3d,
[  --without-3d ............................. turn off instantiation of the 3D algorithms],
[
   AC_MSG_RESULT(yes)
   INST3D="no"
],
[
   AC_MSG_RESULT(no)
   INST3D="yes"
   DIMS+="3 "
]
)

])

