# -*-makefile-*-
# Spheral++ ExternalForce BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = ExternalForce
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	PointPotential.pyste \
	ConstantAcceleration.pyste \
	LinearAcceleration.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
