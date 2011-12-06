# -*-makefile-*-
# Spheral++ Integrator BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Integrator
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	Integrator.pyste \
	PredictorCorrector.pyste \
	SynchronousRK1.pyste \
	SynchronousRK2.pyste \
	CheapSynchronousRK2.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
