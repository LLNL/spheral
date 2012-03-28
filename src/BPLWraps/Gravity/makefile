# -*-makefile-*-
# Spheral++ Gravity BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Gravity
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	NBodyGravity.pyste
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
