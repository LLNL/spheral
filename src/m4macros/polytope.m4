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
AC_SUBST(PYTHONPKGS)

POLYTOPEFLAGS=
POLYTOPELIBS=

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
   POLYTOPEFLAGS+=" -Wno-dev  -DCMAKE_INSTALL_PREFIX:PATH=\$(prefix) -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=\$(CC)  -DCMAKE_CXX_COMPILER=\$(CXX) -DBUILD_SHARED_LIBS=ON -DTESTING=OFF -DBOOST_ROOT=\$(prefix) -DHDF5_ROOT=\$(prefix) -DUSE_SILO=ON -DUSE_PYTHON=ON -DPYTHON_EXE=\$(PYTHONEXE) -G 'Unix Makefiles'"
   POLYTOPELIBS+=" -lpolytope"
   PYTHONPKGS+=" polytope"
]
)

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
# Optionally build polytope with Tetgen
#
# Defaulting off 'til the current polytope supports Tetgen again.
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-tetgen)
AC_ARG_WITH(tetgen,
[  --with-tetgen ............................ optionally build polytope with Tetgen],
[
   AC_MSG_RESULT(yes)
   USE_TETGEN=1
   POLYTOPELIBS+=" -ltetgen"
],
[
   AC_MSG_RESULT(no)
   USE_TETGEN=0
]
)

])

