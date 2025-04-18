###############################################
Building Spheral on Livermore Computing systems
###############################################

This guide explains the process for setting up and running Spheral on LC (Livermore Computing) systems.

.. warning::

   Do not run the provided terminal commands on login nodes. Make sure to grab an allocation first, using ``salloc -N 2`` for SLURM machines or ``flux alloc -xN 2`` for Flux machines. Do not use the ``srun`` or ``flux run`` commands when running ``tpl-manager.py`` or ``spheral-ats``.

Cloning/Updating Spheral
========================

.. include:: ../include/cloning.rst.inc


Third Party Libraries (TPLs)
============================

Spheral uses Spack to build and install TPLs used by Spheral. Spheral includes Spack configurations and environments pre-configured for LC systems in ``scripts/spack/configs/$SYS_TYPE`` and ``scripts/spack/environments/$SYS_TYPE``. This greatly simplifies building TPLs for Spheral.

.. include:: ../include/tpls.rst.inc

Configuring
===========

.. include:: ../include/configure.rst.inc

Build and Install
=================

.. include:: ../include/building.rst.inc

Running Tests
=============

.. include:: ../include/tests.rst.inc
