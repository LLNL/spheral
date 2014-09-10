# -*-makefile-*-
# Spheral++ Material package makefile
#-------------------------------------------------------------------------------
PKGNAME = Material
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	PhysicalConstants.cc \
	EquationOfStateInst.cc \
	GammaLawGasInst.cc \
	PolytropicEquationOfStateInst.cc \
	IsothermalEquationOfStateInst.cc \
	$(HELMSRCS)

F90SRCTARGETS = \
	$(HELMFSRCS)

PYTHONTARGETS = \
	MaterialUnits.py \
	MaterialEquationsOfState.py

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)
EXTRALIBS = $(HELMLIBS)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
