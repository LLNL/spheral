# -*-makefile-*-
# Spheral++ CSPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = CSPH
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeCSPHSumMassDensityInst.cc \
	computeHullSumMassDensityInst.cc \
	computeCSPHCorrectionsInst.cc \
	computeCSPHEvaluationInst.cc \
	computeHullVolumesInst.cc \
	computeHVolumesInst.cc \
	computeNeighborHullInst.cc \
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
