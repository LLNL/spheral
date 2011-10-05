# -*-makefile-*-
# Spheral++ Physics BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Physics
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	Physics.pyste \
	GenericHydro.pyste \
	GenericBodyForce.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
