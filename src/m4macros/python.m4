# ------------------------------------------------------------------
# Python Macros
# ------------------------------------------------------------------


# -----------------------------------------------
# Look for python (override with --with-python or
# PYTHON in environment.  
# -----------------------------------------------

AC_DEFUN([AC_PROG_PYTHON],[

AC_SUBST(PYTHONROOT)
AC_SUBST(PYTHON)
AC_SUBST(PYTHONVERSION)
AC_SUBST(PYLIBDIR)
AC_SUBST(TOPLIBDIR)
AC_SUBST(LIBDIR)

AC_MSG_CHECKING(for --with-python)
AC_ARG_WITH(python,[  --with-python=/usr/local/2.6/bin/python .. use non-standard python],[
  AC_MSG_RESULT($withval)
  PYTHON=$withval
  PYTHONVERSION=`${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_version()"`
  PYLIBDIR=`${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_lib()"`
  TOPLIBDIR=`${PYTHON} -c "import distutils.sysconfig, os.path; print os.path.split(os.path.split(distutils.sysconfig.get_python_lib())[[0]])[[0]]"`
  LIBDIR=$PYLIBDIR/Spheral
],[
  AC_MSG_RESULT(no)
  PYTHONROOT="\$(prefix)"
  PYTHON="\$(exec_prefix)/bin/python"
  PYTHONVERSION=2.7
  PYLIBDIR=$PYTHONROOT/lib/python$PYTHONVERSION/site-packages
  TOPLIBDIR=$PYTHONROOT/lib
  LIBDIR=$PYLIBDIR/Spheral
])

])
