# -----------------------------------------------------------------
# MPI
# -----------------------------------------------------------------
AC_DEFUN([SETUP_MPI],[
AC_SUBST(CXXFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(CXXPKGS)
AC_SUBST(CXXPKGLIBS)
AC_SUBST(PYTHONPKGS)
AC_SUBST(MPIPYTHONINTERFACETARGET)
AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(HDF5FLAGS)
AC_SUBST(POLYTOPEFLAGS)
AC_SUBST(EXTRATHIRDPARTYTARGETS)

HDF5FLAGS=

AC_MSG_CHECKING(for --without-mpi)
AC_ARG_WITH(mpi,
[  --without-mpi ............................ turn off mpi],
[
    AC_MSG_RESULT(yes)
    MPIPYTHONINTERFACETARGET="fakempi.py"
    MPIENABLED="no"
    MPICC=$CC
    MPICXX=$CXX
    CXXPKGS="$CXXPKGS PythonMPIInterfaces"
    POLYTOPEFLAGS="$POLYTOPEFLAGS CC='\$(CC)' CXX='\$(CXX)'"
],
[
    AC_MSG_RESULT(no)
    EXTRATHIRDPARTYTARGETS+=" .mpi4py-1.3.date"
    MPIPYTHONINTERFACETARGET="mpi_mpi4py.py"
    #HDF5FLAGS="$HDF5FLAGS --enable-parallel"
    CXXFLAGS="$CXXFLAGS -DUSE_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DOMPI_SKIP_MPICXX"
    CFLAGS="$CFLAGS -DUSE_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DOMPI_SKIP_MPICXX"
    CPPFLAGS="$CPPFLAGS -DUSE_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DOMPI_SKIP_MPICXX"
    CXXPKGS="$CXXPKGS Distributed PythonMPIInterfaces"
    CXXPKGLIBS="$CXXPKGLIBS Distributed"
    PYTHONPKGS="$PYTHONPKGS Distributed"
    MPICC=$CC
    MPICXX=$CXX
    MPIENABLED="yes"
    POLYTOPEFLAGS="$POLYTOPEFLAGS CC='\$(MPICC)' CXX='\$(MPICXX)' CFLAGS='\$(MPICCFLAGS)' CXXFLAGS='\$(MPICXXFLAGS)' MPI=1"

    # # On Apple we will exclude the C++ bindings.
    # if test "`uname -s`" = "Darwin"; then
    #   CFLAGS="$CFLAGS -DOMPI_SKIP_MPICXX"
    #   CPPFLAGS="$CPPFLAGS -DOMPI_SKIP_MPICXX"
    # fi
])

# -----------------------------------------------------------------
# Deadlock detection
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-deadlockDetection)
AC_ARG_WITH(deadlockDetection,
[  --with-deadlockDetection ................. enable local MPI deadlock detection],
[
    CXXFLAGS="$CXXFLAGS -DUSE_MPI_DEADLOCK_DETECTION"
    CFLAGS="$CFLAGS -DUSE_MPI_DEADLOCK_DETECTION"
    CPPFLAGS="$CPPFLAGS -DUSE_MPI_DEADLOCK_DETECTION"
    AC_MSG_RESULT(yes)
],
[
    AC_MSG_RESULT(no)
])

])

# -----------------------------------------------------------------
# ParMETIS
# -----------------------------------------------------------------
AC_DEFUN([SETUP_PARMETIS],[
AC_SUBST(PARMETISPATH)
AC_SUBST(PARMETISTARGET)
AC_SUBST(PARMETISBPLTARGET)
AC_SUBST(PARMETISINCS)
AC_SUBST(PARMETISLIBS)
PARMETISPATH=$SPHERALDIR
AC_MSG_CHECKING(for --with-parmetis)
AC_ARG_WITH(parmetis,
[  --with-parmetis .......................... compile with ParMETIS],
[
    EXTRATHIRDPARTYTARGETS+=" .parmetis-4.0.3.date"
    PARMETISPATH=$SPHERALDIR
    PARMETISTARGET="ParmetisRedistributeNodesInst.cc"
    PARMETISBPLTARGET="ParmetisRedistributeNodes.pyste"
    PARMETISINCS="-I$PARMETISPATH/include"
    PARMETISLIBS="-L$PARMETISPATH/lib -lparmetis"
    PYTHONPKGS="$PYTHONPKGS Parmetis"
    AC_MSG_RESULT(yes)
],
[
    PARMETISPATH="$SPHERALDIR"
    PARMETISTARGET=
    PARMETISBPLTARGET=
    PARMETISINCS=
    PARMETISLIBS=
    AC_MSG_RESULT(no)
])
])

