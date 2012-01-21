# -*-makefile-*-
# Spheral++ Material BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Material
PKGDIR = $(PKGNAME)
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	PhysicalConstants.pyste \
	EquationOfState.pyste \
	GammaLawGas.pyste \
	PolytropicEquationOfState.pyste \
	IsothermalEquationOfState.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
