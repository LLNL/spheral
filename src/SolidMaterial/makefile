# -*-makefile-*-
# SolidSpheral++ SolidMaterial package makefile
#-------------------------------------------------------------------------------
PKGNAME = SolidMaterial
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	LinearPolynomialEquationOfStateInst.cc \
	GruneisenEquationOfStateInst.cc \
	MurnahanEquationOfStateInst.cc \
	TillotsonEquationOfStateInst.cc \
	SteinbergGuinanStrengthInst.cc \
	SteinbergGuinanLundStrengthInst.cc \
	PorousEquationOfStateInst.cc \
	StrainPorosityInst.cc

PYTHONTARGETS = \
	SolidMaterialUnits.py \
	SolidMaterialEquationsOfState.py

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SOLIDTOP)/helpers/makefile_master
