# -----------------------------------------------------------------
# Polytope
# -----------------------------------------------------------------
AC_DEFUN([SETUP_POLYTOPE],[
AC_SUBST(USE_POLYTOPE)
AC_SUBST(CXXFLAGS)
AC_SUBST(EXTRATHIRDPARTYTARGETS)
AC_SUBST(POLYTOPEFLAGS)
AC_SUBST(POLYTOPELIBS)
AC_SUBST(USE_TRIANGLE)
AC_SUBST(USE_TETGEN)

POLYTOPEFLAGS=
POLYTOPELIBS=

# -----------------------------------------------------------------
# Optionally build polytope without Triangle
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-triangle)
AC_ARG_WITH(triangle,
[  --without-triangle ....................... optionally build polytope without Triangle],
[
   AC_MSG_RESULT(yes)
   USE_TRIANGLE=0
],
[
   AC_MSG_RESULT(no)
   USE_TRIANGLE=1
   POLYTOPELIBS+=" -ltriangle"
]
)

# -----------------------------------------------------------------
# Optionally build polytope without Tetgen
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-tetgen)
AC_ARG_WITH(tetgen,
[  --without-tetgen ......................... optionally build polytope without Tetgen],
[
   AC_MSG_RESULT(yes)
   USE_TETGEN=0
],
[
   AC_MSG_RESULT(no)
   USE_TETGEN=1
   POLYTOPELIBS+=" -ltetgen"
]
)

# -----------------------------------------------------------------
# Optionally turn off polytope entirely.  This disables building 
# the tessellation extensions of Spheral.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --without-polytope)
AC_ARG_WITH(polytope,
[  --without-polytope ....................... optionally build without polytope (disables all tessellation extensions)],
[
   AC_MSG_RESULT(yes)
   USE_POLYTOPE=0
   POLYTOPEFLAGS=
   POLYTOPELIBS=
   CXXFLAGS+=" -DNOPOLYTOPE"
],
[
   AC_MSG_RESULT(no)
   USE_POLYTOPE=1
   EXTRATHIRDPARTYTARGETS+=" \$(POLYTOPEBUILDDATE)"
   POLYTOPEFLAGS+=" prefix=\$(prefix) boost_root=\$(prefix) use_python=1 build_tests=0 python_exe=$PYTHON python_version=$PYTHONVERSION"
   POLYTOPELIBS+=" -lpolytope -lvoro_2d -lvoro_3d"
]
)

])

