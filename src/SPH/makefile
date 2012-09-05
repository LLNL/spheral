# -*-makefile-*-
# Spheral++ SPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = SPH
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeSPHSumMassDensityInst.cc \
	computeSPHOmegaGradhCorrectionInst.cc \
	computeSumVoronoiCellMassDensityInst.cc \
	computeHydrostaticEquilibriumPressure.cc \
	SPHHydroBaseInst.cc

PYTHONTARGETS = \
	SPHHydros.py \
	SPHUtilities.py

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
