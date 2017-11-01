# ------------------------------------------------------------------------------
# Set the system dependent distributed make options.
# ------------------------------------------------------------------------------
AC_DEFUN([SETUP_DISTRIBUTEDMAKE], [

AC_SUBST(SUBMITDISTRIBUTEDMAKE)
AC_SUBST(DISTRIBUTEDMAKEOPTS)

# ------------------------------------------------------------------------------
# Set up defaults for various systems I know about.
# ------------------------------------------------------------------------------
SUBMITDISTRIBUTEDMAKE="mpirun -np 5"
DISTRIBUTEDMAKEOPTS="numLocalTasks=1"
case "`uname -n`" in
pengra*)
  SUBMITDISTRIBUTEDMAKE="srun -l -p pdebug -n 4 -N 2"
  ;;
mcr*)
  SUBMITDISTRIBUTEDMAKE="srun -l -p pdebug -n 10 -N 5"
  ;;
bigdev*)
  SUBMITDISTRIBUTEDMAKE="srun -l -p pkull -n 12 -N 6"
  ;;
esac

# ------------------------------------------------------------------------------
# Check for SUBMITDISTRIBUTEDMAKE.
# ------------------------------------------------------------------------------
AC_MSG_CHECKING(for submitDistributedMake)
AC_ARG_WITH(submitDistributedMake,
  [  --with-submitDistributedMake=ARG ......... choose a method of invoking the distribute make],
  [
    SUBMITDISTRIBUTEDMAKE=$withval
  ],
)
AC_MSG_RESULT($SUBMITDISTRIBUTEDMAKE)

# ------------------------------------------------------------------------------
# Check for any optional DISTRIBUTEDMAKEOPTS.
# ------------------------------------------------------------------------------
AC_MSG_CHECKING(for distributedMakeOpts)
AC_ARG_WITH(distributedMakeOpts,
  [  --with-distributedMakeOpts=ARG ........... set optional arguements to DistributedMake (see spheral/src/helpers/DistributedMake)],
  [
    DISTRIBUTEDMAKEOPTS=$withval
  ],
)
AC_MSG_RESULT($DISTRIBUTEDMAKEOPTS)


])
