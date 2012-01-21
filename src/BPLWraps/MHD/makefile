# -*-makefile-*-
# Spheral++ MHD BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = MHD
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	ConductingFluidNodeList.pyste \
	PriceMonaghanDissipation.pyste \
	MHD.pyste 

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I. -I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) \
	-I$(SPHERALTOP)/../../include/pyMPI2.4 -I$(SPHERALTOP)/../../include/python2.5

