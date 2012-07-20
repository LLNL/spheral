# -*-makefile-*-
# SolidSpheral++ SolidMaterial package makefile
#-------------------------------------------------------------------------------
PKGNAME = SolidMaterial
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	LinearPolynomialEquationOfState.cc \
	GruneisenEquationOfState.cc \
	MurnahanEquationOfState.cc \
	TillotsonEquationOfState.cc \
	SteinbergGuinanStrength.cc \
	SteinbergGuinanLundStrength.cc \
	PorousEquationOfState.cc \
	StrainPorosity.cc

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SOLIDTOP)/helpers/makefile_master
