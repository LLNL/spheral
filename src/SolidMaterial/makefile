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
	StrainPorosityInst.cc

PYTHONTARGETS = \
	SolidMaterialUnits.py \
	SolidMaterialEquationsOfState.py

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SOLIDTOP)/helpers/makefile_master
