# -*-makefile-*-
# Spheral++ DataOutput BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = DataOutput
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	RestartableObject.pyste \
	RestartRegistrar.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
