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
	HelmholtzEquationOfStateInst.cc

F90SRCTARGETS = \
	public_helm.f90 \
	invert_helm.f90

PYTHONTARGETS = \
	MaterialUnits.py \
	MaterialEquationsOfState.py

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
