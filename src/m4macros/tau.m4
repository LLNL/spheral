# -----------------------------------------------------------------
# Tau performance profiling library.
# -----------------------------------------------------------------
AC_SUBST(TAUFLAGS)
AC_SUBST(TAULIBS)
AC_SUBST(TAUVERSION)
AC_SUBST(TAUTARGET)
AC_SUBST(TAUCONFIGUREFLAGS)
AC_SUBST(TAUMAKEFILE)
AC_SUBST(TAUARCH)
AC_SUBST(TAUCC)
AC_SUBST(TAUCXX)
AC_SUBST(PAPIROOT)
AC_SUBST(TAUPAPIFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(OPT)
AC_SUBST(DISTRIBUTEDOPT)

AC_DEFUN([SETUP_TAU],[
TAUVERSION=2.21.4
AC_MSG_CHECKING(for --with-tau)
AC_ARG_WITH(tau,
[  --with-tau ............................... turn on Tau class profiling],
[
  AC_MSG_RESULT(yes)
  TAUFLAGS="\$(TAU_INCLUDE) \$(TAU_MPI_INCLUDE) \$(TAU_DEFS) -DPROFILING_ON"
  TAULIBS="\$(TAU_MPI_LIBS) \$(TAU_LDFLAGS) \$(TAU_SHLIBS)"
  TAUTARGET=".tau-\$(TAUVERSION).date"

  OSNAME=`uname -s`
  OSMACH=`uname -m`
  case $OSNAME in
    AIX)
      TAUARCH="rs6000"
      if test $MPIENABLED = "yes"; then
        TAUCONFIGUREFLAGS="-mpiinc=/usr/lpp/ppe.poe/include -mpilib=/usr/lpp/ppe.poe/lib"
      fi
      ;;
   Linux)
     if test $MPIENABLED = "yes"; then
       TAUCONFIGUREFLAGS="-mpiinc=/usr/lib/mpi/include -mpilib=/usr/lib/mpi/lib -useropt=-O1"
     fi
     case $OSMACH in
       ia64)
         TAUARCH="ia64"
         ;;
       x86_64)
         TAUARCH="x86_64"
         TAUCONFIGUREFLAGS="-mpiinc=/usr/local/tools/mvapich-gnu/include -mpilib=/usr/local/tools/mvapich-gnu/lib/shared -useropt=-O1"
         ;;
       *)
         TAUARCH="i386_linux"
         ;;
     esac
     #TAUCONFIGUREFLAGS="$TAUCONFIGUREFLAGS -LINUXTIMERS"  # Currently does not work on LC opterons  :(
     ;;
   Darwin)
     TAUARCH="apple"
     ;;
   *)
     echo "WARNING: unable to determine TAU architecture."
     ;;
  esac

  if test $MPIENABLED = "yes"; then
    TAUCONFIGUREFLAGS="$TAUCONFIGUREFLAGS -mpi"
  else
    TAUCONFIGUREFLAGS="$TAUCONFIGUREFLAGS"
  fi
  TAUMAKEFILE="include \$(SPHERALTOP)/thirdPartyLibs/tau-$(TAUVERSION)/Makefile.tau"
],
[
  AC_MSG_RESULT(no)
  TAUTARGET=".tau-\$(TAUVERSION).dummydate"
  TAUFLAGS="-I \$(SPHERALTOP)/thirdPartyLibs/tau-\$(TAUVERSION)/include"
  TAUMAKEFILE="TAU_INCLUDE = \$(TAUFLAGS)"
])

AC_MSG_CHECKING(for --with-tauCC)
AC_ARG_WITH(tauCC,
[  --with-tauCC=ARG ......................... set the CC compiler for Tau],
[
  AC_MSG_RESULT($withval)
  TAUCC=$withval
],
[
  TAUCC=$CC
  AC_MSG_RESULT($TAUCC)
])

AC_MSG_CHECKING(for --with-tauCXX)
AC_ARG_WITH(tauCXX,
[  --with-tauCXX=ARG ........................ set the C++ compiler for Tau],
[
  AC_MSG_RESULT($withval)
  TAUCXX=$withval
],
[
  TAUCXX=$CXX
  AC_MSG_RESULT($TAUCXX)
])

AC_MSG_CHECKING(for --with-papiroot)
AC_ARG_WITH(papiroot,
[  --with-papiroot=ARG ...................... set the base path for PAPI libraries and headers],
[
  AC_MSG_RESULT($withval)
  PAPIROOT=$withval
],
[
  PAPIROOT="/usr/local/tools/papi"
  AC_MSG_RESULT($PAPIROOT)
])

AC_MSG_CHECKING(for --with-papi)
AC_ARG_WITH(papi,
[  --with-papi .............................. activate PAPI hooks for Tau profiling],
[
  AC_MSG_RESULT(yes)
  TAUCONFIGUREFLAGS="$TAUCONFIGUREFLAGS -papi=${PAPIROOT} -PAPIWALLCLOCK -MULTIPLECOUNTERS"
],
[
  AC_MSG_RESULT(no)
])

# Set the compiler choices to TAU.
TAUCONFIGUREFLAGS="$TAUCONFIGUREFLAGS -cc=\$(TAUCC) -c++=\$(TAUCXX)"

AC_MSG_CHECKING(for --with-tauflags)
AC_ARG_WITH(tauflags,
[  --with-tauflags=ARG ...................... send additional configure parameters to Tau],
[
  TAUCONFIGUREFLAGS="$TAUCONFIGUREFLAGS $withval"
  AC_MSG_RESULT($TAUCONFIGUREFLAGS)
],
[
  TAUCONFIGUREFLAGS="$TAUCONFIGUREFLAGS -PROFILE -PROFILECALLPATH"
  AC_MSG_RESULT(Defaulting to $TAUCONFIGUREFLAGS)
])



])
