# =======================================================================
# Debug flags
# =======================================================================
AC_DEFUN([AC_CHECK_OPTIMIZATION], [
AC_SUBST(OPT)
AC_SUBST(PYFFLE_OPT)
AC_SUBST(DISTRIBUTEDOPT)
AC_SUBST(BPLOPT)


# =======================================================================
# Check for the optimization.
# =======================================================================
AC_MSG_CHECKING(for --with-opt=)
AC_ARG_WITH(opt,
[  --with-opt=Val ........................... set optimization level (0,1,2,3) ],
[   
case $CXXCOMPILERTYPE in
GNU)
  if test "$withval" = "0";then
    OPT="-g"
  elif test "$withval" = "1";then
    OPT="-O1"
  elif test "$withval" = "10" -o "$withval" = "01";then
    OPT="-O1-g"
  elif test "$withval" = "2";then
    OPT="-O2"
  elif test "$withval" = "20" -o "$withval" = "02";then
    OPT="-O2-g"
  elif test "$withval" = "3";then
    OPT="-O3"
  elif test "$withval" = "4";then
    OPT="-O4-march=native"
  elif test "$withval" = "5";then
    OPT="-O5"
  elif test "$withval" = "6";then
    OPT="-O6"
  elif test "$withval" = "03" -o "$withval" = "30"; then
    OPT="-g -O3"
  else
    echo "Unknown optimization level, defaulting to -O"
    OPT="-O"
  fi
  PYFFLE_OPT="-g";;
KAI)
  if test "$withval" = "0";then
    OPT="+K0"
  elif test "$withval" = "1";then
    OPT="+K1 -O1"
    BPLOPT="$BPLOPT -O"
  elif test "$withval" = "2";then
    OPT="+K2 -O2"
    CXXFLAGS="$CXXFLAGS"
    CFLAGS="$CFLAGS"
    BPLOPT="$BPLOPT -O"
  elif test "$withval" = "3";then
    OPT="+K3 -O3"
    BPLOPT="$BPLOPT -O"
  else    
    echo "Unknown optimization level, defaulting to -O"
    OPT="-O"
  fi
  PYFFLE_OPT="+K0";;
INTEL)
  if test "$withval" = "0";then
    OPT="-g"
  elif test "$withval" = "1";then
    OPT="-O"
  elif test "$withval" = "2";then
    OPT="-O2"
  elif test "$withval" = "3";then
    OPT="-O3"
    BPLOPT="-O"
  elif test "$withval" = "4";then
    OPT="-O3 -ip -ansi-alias -no-prec-div"
  elif test "$withval" = "02" -o "$withval" = "20"; then
    OPT="-g -O2"
  elif test "$withval" = "03" -o "$withval" = "30"; then
    OPT="-g -O3"
  else    
    echo "Unknown optimization level, defaulting to -O"
    OPT="-O"
  fi
  BPLOPT="-O0"
  PYFFLE_OPT="-g";;
*)
  if test "$withval" = "0";then
    OPT="-g"
  elif test "$withval" = "1";then
    OPT="-O1"
  elif test "$withval" = "2";then
    OPT="-O2"
  elif test "$withval" = "3";then
    OPT="-O3"
  elif test "$withval" = "4";then
    OPT="-O4"
  elif test "$withval" = "5";then
    OPT="-O5"
  elif test "$withval" = "6";then
    OPT="-O6"
  else    
    echo "Unknown optimization level, defaulting to -O"
    OPT="-O"
  fi
  PYFFLE_OPT="-g";;
esac

if test "$withval" = "0"; then
     DISTRIBUTEDOPT=""
else
     OPT="$OPT -DNDEBUG"
     DISTRIBUTEDOPT=$OPT
fi

AC_MSG_RESULT($OPT)
echo "PYFFLE_OPT is $PYFFLE_OPT"
echo "DISTRIBUTEDOPT is $DISTRIBUTEDOPT"
echo "BPLOPT is $BPLOPT"
], [
AC_MSG_RESULT(no)
OPT="-O3"
PYFFLE_OPT="-g"
DISTRIBUTEDOPT=$OPT
echo "Defaulting optimization to $OPT"
])
echo "CXXFLAGS is $CXXFLAGS"
])
