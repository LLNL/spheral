# -*-makefile-*-
# Spheral++ Material package makefile
#-------------------------------------------------------------------------------
PKGNAME = Material
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	EquationOfState.cc \
	MKSUnits.cc \
	CGSUnits.cc \
	CosmologicalUnits.cc \
	PhysicalConstants.cc \
	GammaLawGas.cc \
	PolytropicEquationOfState.cc \
	IsothermalEquationOfState.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
