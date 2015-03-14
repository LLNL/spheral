# -----------------------------------------------------------------
# PETSc -- A library providing (among other things) linear solvers.
# -----------------------------------------------------------------
AC_DEFUN([SETUP_PETSC],[
AC_SUBST(PETSCTARGETS)
AC_SUBST(PETSCOPTS)
AC_SUBST(MPIENABLED)
AC_SUBST(USEPETSC)

# Set up a string containing build options for PETSc.
PETSCOPTS="--with-shared=1 --with-dynamic=0 --with-fc=0 --with-python=1 --download-c-blas-lapack=1"
if test "${MPIENABLED}" = "yes"; then
   PETSCOPTS="${PETSCOPTS} --with-cc=${MPICC} --with-mpi=1 --download-parmetis=1"
else
   PETSCOPTS="${PETSCOPTS} --with-cc=${CC} --with-mpi=0"
fi

AC_MSG_CHECKING(for --with-petsc)
AC_ARG_WITH(petsc,
[  --with-petsc ............................. build the PETSc libraries],
[
    AC_MSG_RESULT(yes)
    PETSCTARGETS=".petsc-2.3.3-p8.date .petsc4py-0.7.5.date"
],
[
    if test -n "${USEPETSC}"; then
    AC_MSG_RESULT(yes, needed by something)
    PETSCTARGETS=".petsc-2.3.3-p8.date .petsc4py-0.7.5.date"
    else
    AC_MSG_RESULT(no)
    PETSCTARGETS=""
    fi
])

])

