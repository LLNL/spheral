# ------------------------------------------------------------------
# Python Macros
# ------------------------------------------------------------------


# -----------------------------------------------
# Look for python (override with --with-python or
# PYTHON in environment.  
# -----------------------------------------------

AC_DEFUN([AC_PROG_PYTHON],[

AC_SUBST(PYTHONTARGET)
AC_SUBST(PYTHONROOT)
AC_SUBST(PYTHON)
AC_SUBST(PYTHONVERSION)
AC_SUBST(TOPLIBDIR)
AC_SUBST(LIBDIR)
AC_SUBST(PYLIBDIR)
AC_SUBST(TPINCS)

PYTHONROOT="\$(prefix)"
TPINCS=

AC_MSG_CHECKING(for --with-python)
AC_ARG_WITH(python,[  --with-python=/usr/local/2.6/bin/python .. use non-standard python],[
  AC_MSG_RESULT($withval)
  PYTHON=$withval
  PYTHONVERSION=`${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_version()"`
  PYTHONTARGET=${PYTHON}
  TPINCS+=" -I `${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()"`"
],[
  AC_MSG_RESULT(no)
  PYTHON="\$(exec_prefix)/bin/python"
  PYTHONVERSION=2.7
  PYTHONTARGET=".Python-2.7.15.date"
])

PYLIBDIR=$PYTHONROOT/lib/python$PYTHONVERSION/site-packages
TOPLIBDIR=$PYTHONROOT/lib
LIBDIR=$TOPLIBDIR

])
