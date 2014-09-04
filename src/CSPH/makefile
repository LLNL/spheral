# -*-makefile-*-
# Spheral++ CSPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = CSPH
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeCSPHSumMassDensityInst.cc \
	computeCSPHCorrectionsInst.cc \
	computeHullVolumesInst.cc \
	computeHVolumesInst.cc \
	centerOfMass.cc \
	CSPHHydroBaseInst.cc \
	interpolateCSPHInst.cc \
	gradientCSPHInst.cc \
	HVolumePolicyInst.cc \
#	CSPHSpecificThermalEnergyPolicyInst.cc

PYTHONTARGETS = \
	CSPHHydros.py

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
