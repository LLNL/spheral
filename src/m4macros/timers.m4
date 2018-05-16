dnl -----------------------------------------------------------------
dnl Timers
dnl -----------------------------------------------------------------

AC_DEFUN([SETUP_TIMERS],[
AC_SUBST(TIMERTARGETS)

AC_MSG_CHECKING(for --with-timers)
AC_ARG_WITH(timers,
[  --with-timers ............................ turn on the Timer class profiling],
[
  AC_MSG_RESULT(yes)
  CXXFLAGS="$CXXFLAGS -DTIMER"
  CFLAGS="$CFLAGS -DTIMER"
  CPPFLAGS="$CPPFLAGS -DTIMER"
  TIMERTARGETS="Timer.cc SpheralTimers.cc"
],
[
  AC_MSG_RESULT(no)
  TIMERTARGETS="Timer.cc SpheralTimers.cc"
])
])

AC_SUBST(PAPIROOT)
AC_SUBST(PAPILIBS)
AC_DEFUN([SETUP_PAPI],[

AC_MSG_CHECKING(for --with-papiroot)
AC_ARG_WITH(papiroot,
[  --with-papiroot=ARG ...................... set the base path for PAPI libraries and headers],
[
  AC_MSG_RESULT($withval)
  PAPIROOT=$withval
],
[
  AC_MSG_RESULT(no)
  PAPIROOT=/usr/local/papi
])

AC_MSG_CHECKING(for --with-papi)
AC_ARG_WITH(papi,
[  --with-papi .............................. activate PAPI for class profiling],
[
  AC_MSG_RESULT(yes)
  CXXFLAGS="$CXXFLAGS -DPAPI -I$PAPIROOT/include"
  CFLAGS="$CFLAGS -DPAPI"
  CPPFLAGS="$CPPFLAGS -DPAPI"
  EXTRAFLAGS="$EXTRAFLAGS -I$PAPIROOT/include"
  PAPILIBS="$PAPIROOT/lib/libpapi.a -lpmapi"
  #PAPILIBS="-L/usr/local/papi/lib -lpapi -L/usr/local/lib -lperfctr"
],
[
  AC_MSG_RESULT(no)
  PAPILIBS=
])

])
