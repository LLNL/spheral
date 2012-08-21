# -*-makefile-*-
# Spheral++ Material package makefile
#-------------------------------------------------------------------------------
PKGNAME = Material
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	PhysicalConstants.cc \
	EquationOfStateInst.cc \
	GammaLawGasInst.cc \
	PolytropicEquationOfStateInst.cc \
	IsothermalEquationOfStateInst.cc

PYTHONTARGETS = \
	MaterialUnits.py \
	MaterialEquationsOfState.py

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
