# ------------------------------------------------------------------
# Python Macros
# ------------------------------------------------------------------


# -----------------------------------------------
# Look for python (override with --with-python or
# PYTHON in environment.  
# -----------------------------------------------

AC_DEFUN([AC_PROG_PYTHON],[

if test x$prefix != xNONE; then
  PATH="$prefix/bin:$PATH" 
fi
AC_SUBST(PYTHON)
AC_SUBST(PYTHONVERSION)
AC_SUBST(PYTHONROOT)
AC_SUBST(PYTHONSITEPKGDIR)
AC_SUBST(PYLIBDIR)

AC_MSG_CHECKING(for --with-python)
AC_ARG_WITH(python,[  --with-python=python1.6a2 ................ use non-standard python],[
  AC_MSG_RESULT($withval)
  PYTHON=$withval
],[
  AC_MSG_RESULT(no)
  # If we're using Mac OS X, we need to use the version that shipped with the 
  # operating system...
  PYTHONVERSION=2.7
#   if test "`uname -s`" = "Darwin"; then
#     PYTHONROOT=/System/Library/Frameworks/Python.framework/Versions/$PYTHONVERSION
#     # Otherwise, proceed normally.
#   else 
#     PYTHONROOT=`echo $PWD | sed -e "s/spheral\/src$//g;"`
#   fi
  PYTHONROOT=`echo $PWD | sed -e "s/spheral\/src$//g;"`
  PYTHON=$PYTHONROOT/bin/python
  PYTHONSITEPKGDIR=$PYTHONROOT/lib/python$PYTHONVERSION/site-packages
  PYLIBDIR=$PYTHONSITEPKGDIR
#   if test -z "$PYTHON"; then
#     AC_PATH_PROG(PYTHON,python)
#   fi
])

#AC_MSG_CHECKING(for execute access)
#if test ! -x "$PYTHON"; then
#  AC_MSG_RESULT(no)
#  AC_MSG_ERROR(Python in '$PYTHON' not executable)
#else
#  AC_MSG_RESULT(yes)
#fi
])


# AC_DEFUN([SETUP_PYTHON_ENV],[
# # -----------------------------------------------
# # Check Python version 
# # -----------------------------------------------
# AC_MSG_CHECKING(python version string)
# AC_SUBST(PYVERSION)
# PYVERSION=`$PYTHON -c "import sys,re; print re.compile(r'\\d+\\.\\d+').match(sys.version).group()" 2>/dev/null`
# if test -z "$PYVERSION"; then
#   AC_MSG_RESULT(1.4 or earlier)
#   AC_MSG_ERROR(Python must be 1.4 or higher)
# else
#   AC_MSG_RESULT($PYVERSION)
# fi

## -----------------------------------------------
## Disutils?
## -----------------------------------------------
#AC_MSG_CHECKING(distutils?)
#PYDISTUTILS=`$PYTHON -c "import distutils; print 'yes'" 2>/dev/null`
#AC_MSG_RESULT($PYDISTUTILS)

#if test -z "$PYDISTUTILS"; then 
#	AC_MSG_ERROR(Must have distutils to configure Spheral)	
#fi

# # -----------------------------------------------
# # Checks for python header files.
# # -----------------------------------------------

# AC_MSG_CHECKING(for python include location)
# AC_SUBST(PYINCLUDEDIR)
# PYINCLUDEDIR=`$PYTHON -c "from distutils.sysconfig import get_python_inc; print get_python_inc()" 2>/dev/null`
# AC_MSG_RESULT($PYINCLUDEDIR)

# AC_MSG_CHECKING(that include directory exists)
# if test -d "$PYINCLUDEDIR" -a -r "$PYINCLUDEDIR" -a -x "$PYINCLUDEDIR"; then
#   AC_MSG_RESULT(yes)
#   CPPFLAGS="$CPPFLAGS -I$PYINCLUDEDIR"
# else
#   AC_MSG_RESULT(no)
#   AC_MSG_ERROR(No installed headers for version $PYVERSION)
# fi

# # -----------------------------------------------
# # Checks for python lib files...
# # -----------------------------------------------
# AC_MSG_CHECKING(for python library location)
# AC_SUBST(PYLIBDIR)
# PYLIBDIR=`$PYTHON -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" 2>/dev/null`
# AC_MSG_RESULT($PYLIBDIR)

# AC_MSG_CHECKING(that lib directory is accessable)
# if test -d "$PYLIBDIR" -a -r "$PYLIBDIR" -a -x "$PYLIBDIR"; then
#   AC_MSG_RESULT(yes)
# else
#   AC_MSG_RESULT(no)
#   AC_MSG_WARN(No installed lib for version $PYVERSION)
# fi

## -----------------------------------------------
## Checks for python lib/config files...
## -----------------------------------------------
#AC_MSG_CHECKING(for python lib/config location)
#AC_SUBST(PYLIBCONFIGDIR)
#if test -z "$PYDISTUTILS"; then
#  PYLIBCONFIGDIR=$PYLIBDIR/config
#else
#  changequote(X,Y)
#  PYLIBCONFIGDIR=`$PYTHON -c "from distutils.sysconfig import get_makefile_filename; import os; print os.path.split(get_makefile_filename())[0]" 2>/dev/null`
#  changequote([,])
#fi
#AC_MSG_RESULT($PYLIBCONFIGDIR)

#AC_MSG_CHECKING(that lib/config directory is accessible)
#if test -d "$PYLIBCONFIGDIR" -a -r "$PYLIBCONFIGDIR" -a -x "$PYLIBCONFIGDIR"; then
#  AC_MSG_RESULT(yes)
#  LIBS="$LIBS -L$PYLIBCONFIGDIR"
#else
#  AC_MSG_RESULT(no)
#  AC_MSG_WARN(No installed lib/config for version $PYVERSION)
#fi

## -----------------------------------------------
## Checks for python directory for libpythonx.x.a
## -----------------------------------------------
#AC_MSG_CHECKING(libpython$PYVERSION is there)
#if test -r "$PYLIBCONFIGDIR/libpython$PYVERSION.a"; then
#  AC_MSG_RESULT(yes)
#  LIBS="$LIBS -lpython$PYVERSION"
#else
#  AC_MSG_RESULT(not there)
#  AC_MSG_WARN(libpython$PYVERSION.a wasn't where I thought it was)
#  AC_MSG_WARN(You may need to add a --with-libs=-L/somewhere/else or --with-distribution=/path/to/pythonstuff)
#fi

## -----------------------------------------------
## See if there exists a python.exp file
## if so, put it in the right place (src/)
## for the final link
## -----------------------------------------------
#AC_MSG_CHECKING(python.exp file)
#if test -r $PYLIBCONFIGDIR/python.exp; then
#  rm -f python.exp
#  # copy to current directory to see if link succedes 
#  cp $PYLIBCONFIGDIR/python.exp . 
#  AC_MSG_RESULT($PYLIBCONFIGDIR/python.exp)
#else
#  AC_MSG_RESULT(no)
#fi

## -----------------------------------------------
## See if we can steal the Python libraries...
## -----------------------------------------------
#AC_MSG_CHECKING(PYLDFLAGS)
#AC_SUBST(PYLDFLAGS)
#changequote(Z,W)
#PYLDFLAGS=`$PYTHON -c "from distutils.sysconfig import parse_makefile, get_makefile_filename; print parse_makefile(get_makefile_filename())['LDFLAGS']" `
#changequote([,])
#AC_MSG_RESULT($PYLDFLAGS)

#AC_MSG_CHECKING(PYLINKFORSHARED)
#AC_SUBST(PYLINKFORSHARED)
#changequote(Z,W)
#PYLINKFORSHARED=`$PYTHON -c "from distutils.sysconfig import parse_makefile, get_makefile_filename; print parse_makefile(get_makefile_filename())['LINKFORSHARED']" `
#changequote([,])
#AC_MSG_RESULT($PYLINKFORSHARED)

#AC_MSG_CHECKING(PYLOCALMODLIBS)
#AC_SUBST(PYLOCALMODLIBS)
#changequote(Z,W)
#PYLOCALMODLIBS=`$PYTHON -c "from distutils.sysconfig import parse_makefile, get_makefile_filename; print parse_makefile(get_makefile_filename())['LOCALMODLIBS']" `
#changequote([,])
#AC_MSG_RESULT($PYLOCALMODLIBS)

#AC_MSG_CHECKING(PYBASEMODLIBS)
#AC_SUBST(PYBASEMODLIBS)
#changequote(Z,W)
#BASEMODLIBS=`$PYTHON -c "from distutils.sysconfig import parse_makefile, get_makefile_filename; print parse_makefile(get_makefile_filename())['BASEMODLIBS']" `
#changequote([,])
#AC_MSG_RESULT($PYBASEMODLIBS)

# AC_MSG_CHECKING(PYLIBS)
# AC_SUBST(PYLIBS)
# changequote(Z,W)
# PYLIBS=`$PYTHON -c "from distutils.sysconfig import parse_makefile, get_makefile_filename; print parse_makefile(get_makefile_filename())['LIBS']" `
# changequote([,])
# AC_MSG_RESULT($PYLIBS)

#AC_MSG_CHECKING(PYLDLAST)
#AC_SUBST(PYDLAST)
#changequote(Z,W)
#PYLDLAST=`$PYTHON -c "from distutils.sysconfig import parse_makefile, get_makefile_filename; print parse_makefile(get_makefile_filename())['LDLAST']" `
#changequote([,])
#AC_MSG_RESULT($PYLDLAST)

#AC_MSG_CHECKING(Python library options)
#AC_SUBST(LIBS)
#LIBS="$LIBS $PYLDFLAGS $PYLOCALMODLIBS $PYBASEMODLIBS $pyLIBS $PYLDLAST"
## LIBS="$LIBS $PYLDFLAGS $PYLINKFORSHARED $PYLOCALMODLIBS $PYBASEMODLIBS $pyLIBS $PYLDLAST"
#echo "LIBS is $LIBS"
#AC_MSG_RESULT($LIBS)


## -----------------------------------------------
## Checks for required libraries
## -----------------------------------------------
#AC_CHECK_LIB(m,pow)

## need library stuck at the end, so we define what happens on success
#AC_CHECK_LIB(python$PYVERSION,Py_GetBuildInfo,
#[LIBS="$LIBS -lpython$PYVERSION"]
#)

## -----------------------------------------------
## See if we can intuit what extra Python link libs are 
## 
## Fri Jan 26 16:24:40 PST 2001 - if link
## fails, will try to link with libgcc and lm (since
## there is a chance we are running into the
## __eprintf undefined symbol bug 
## -----------------------------------------------
#AC_MSG_CHECKING(Python links as is)
#AC_TRY_LINK(extern int Py_Main(int,char**);, Py_Main(0,0),[
#  AC_MSG_RESULT(yes)
#],[
#  AC_MSG_RESULT(no)

#  # libm and libgcc.a 
#  AC_MSG_CHECKING(link with libgcc.a and -lm)
#  LIBGCC=`/usr/local/bin/gcc --print-libgcc-file-name`
#  LIBS="$LIBS $LIBGCC -lm"
#AC_TRY_LINK(extern int Py_Main(int,char**);, Py_Main(0,0),[
#  AC_MSG_RESULT(yes)
#],[
#  AC_MSG_RESULT(no)
#  AC_MSG_WARN(Python doesn't seem to link.  Look at config.log.  You may need to add --with-libs or --with-includes or just hack the Makefile)
#])
#  # - end check
#])

#AC_SUBST(PYLIBS)
#PYLIBS=$LIBS


])
