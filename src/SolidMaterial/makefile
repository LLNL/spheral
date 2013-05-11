# -*-makefile-*-
# SolidSpheral++ SolidMaterial package makefile
#-------------------------------------------------------------------------------
PKGNAME = SolidMaterial
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	LinearPolynomialEquationOfStateInst.cc \
	GruneisenEquationOfStateInst.cc \
	MurnahanEquationOfStateInst.cc \
	TillotsonEquationOfStateInst.cc \
	StrengthModelInst.cc \
	SteinbergGuinanStrengthInst.cc \
	SteinbergGuinanLundStrengthInst.cc \
	PorousEquationOfStateInst.cc \
	PorousStrengthModelInst.cc \
	StrainPorosityInst.cc \
	$(ANEOSSRCS)

FSRCTARGETS = \
	$(ANEOSFSRCS)

PYTHONTARGETS = \
	SolidMaterialUnits.py \
	SolidMaterialEquationsOfState.py

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

EXTRALIBS = $(ANEOSLIBS)

#-------------------------------------------------------------------------------
include $(SOLIDTOP)/helpers/makefile_master
