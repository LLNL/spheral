# -*-makefile-*-
# Spheral++ Material package makefile
#-------------------------------------------------------------------------------
PKGNAME = Material
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	EquationOfStateInst.cc \
	MKSUnits.cc \
	CGSUnits.cc \
	CosmologicalUnits.cc \
	SolarUnits.cc \
	GammaLawGasInst.cc \
	PolytropicEquationOfStateInst.cc \
	IsothermalEquationOfStateInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
