# ===========================================================================
# File:   configure.ac
# Author: Martin Casado
#
# Copyright (C) 2000 University of California Regents               
# ===========================================================================

# Process this file with autoconf to produce a configure script.
AC_INIT(Physics/Physics.hh)

# All mcaro files are stored in macros/
AC_CONFIG_AUX_DIR(m4macros)

# We want the INSTALL script.
# AC_PROG_INSTALL <-- this doesn't work.
# FIXME: The method below is non-portable... we'll fix this sometime. :-)
AC_SUBST(INSTALL)
INSTALL=`pwd`/m4macros/install-sh

# -------------------------------------------------
# Setup python
# macro definitions are in m4macros/python.m4
# -------------------------------------------------
AC_PROG_PYTHON

# -------------------------------------------------
# Setup the Spheral Environment 
# (m4macros/spheral.m4)
# -------------------------------------------------
SETUP_SPHERAL_ENV

# -------------------------------------------------
# Setup the DistributedMake options
# (m4macros/distributedmake.m4)
# -------------------------------------------------
SETUP_DISTRIBUTEDMAKE

# -------------------------------------------------
# setup the mympKCC script (m4macros/kcc.m4)
# -------------------------------------------------
# SETUP_KCC

# -------------------------------------------------
# setup geommem (m4macros/geommem.m4)
# -------------------------------------------------
SETUP_GEOMMEM

# -------------------------------------------------
# setup polytope (m4macros/polytope.m4)
# -------------------------------------------------
SETUP_POLYTOPE

# -------------------------------------------------
# setup mpi (m4macros/mpi.m4)
# -------------------------------------------------
SETUP_MPI

# -----------------------------------------------
# allow user to define compiler at the command line, plus look
# for KCC first, defined in m4macros/compilers.m4
# -----------------------------------------------
SETUP_COMPILERS_THE_WAY_I_WANT

# ----------------------------------------------
# Setup up compiler flags for spheral depending
# on the compiler
# ----------------------------------------------
AC_CHECK_OPTIMIZATION
SETUP_DEBUG_PRINT
SETUP_DBC

# -------------------------------------------------
# setup ParMETIS (m4macros/mpi.m4)
# -------------------------------------------------
SETUP_PARMETIS

# -------------------------------------------------
# setup postgres support (m4macros/postgres.m4)
# -------------------------------------------------
SETUP_POSTGRES

# -------------------------------------------------
# setup SPH Gravity stuff (m4macros/SPHGravity.m4)
# -------------------------------------------------
SETUP_SPHGRAVITY

# -------------------------------------------------
# setup MHD stuff (m4macros/mhd.m4)
# -------------------------------------------------
SETUP_MHD

# -------------------------------------------------
# setup FractalGravity (m4macros/FractalGravity.m4)
# -------------------------------------------------
SETUP_FRACTALGRAVITY

# -------------------------------------------------
# setup PETSc (m4macros/petsc.m4)
# -------------------------------------------------
SETUP_PETSC

# -------------------------------------------------
# setup python binding (m4macros/pythonbinding.m4)
# -------------------------------------------------
SETUP_PYTHONBINDING

# -------------------------------------------------
# setup matplotlib (m4macros/matplotlib.m4)
# -------------------------------------------------
SETUP_MATPLOTLIB

# -------------------------------------------------
# setup scipy (m4macros/scipy.m4)
# -------------------------------------------------
SETUP_SCIPY

# -------------------------------------------------
# setup ANEOS (m4macros/aneos.m4)
# -------------------------------------------------
SETUP_ANEOS

# -----------------------------------------------
# Output the makefile
# -----------------------------------------------
AC_OUTPUT(helpers/makefile_master \
          makefile \
          PBGWraps/makefile \
          PBGWraps/SpheralModulesBindings.py)
