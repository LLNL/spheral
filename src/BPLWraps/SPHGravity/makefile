# -*-makefile-*-
# Spheral++ SPHGravity BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = SPHGravity
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	SPHGravity.pyste
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I. -I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) \
	-I$(SPHERALTOP)/../../include/pyMPI2.4 -I$(SPHERALTOP)/../../include/python2.5
