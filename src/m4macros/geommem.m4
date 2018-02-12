# -----------------------------------------------------------------
# Select which version of GeomVector we're going to use.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_GEOMMEM],[
AC_SUBST(CXXFLAGS)

# -----------------------------------------------------------------
# Choose which version of the goemetry memory patterns we want to 
# use
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-geommem)
AC_ARG_WITH(geommem,
[  --with-geommem=ARG........................ (default,array,eigen) choose a variant of Geometry internal storage],
[
   GEOMMEM=$withval
],
[
   GEOMMEM="default"
]
)
AC_MSG_RESULT($GEOMMEM)

case $GEOMMEM in
     array)
        CXXFLAGS+=" -DGEOMMEM_ARRAY"
        ;;
     eigen)
        CXXFLAGS+=" -DGEOMMEM_EIGEN"
        ;;
esac

])
