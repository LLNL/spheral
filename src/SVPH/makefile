# -*-makefile-*-
# Spheral++ SVPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = SVPH
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeSVPHCorrectionsInst.cc \
	computeSVPHCorrectionsOnFacesInst.cc \
	computeSumVoronoiCellMassDensityFromFacesInst.cc \
	sampleFieldListSVPHInst.cc \
	gradientFieldListSVPHInst.cc \
	SVPHHydroBaseInst.cc \
	SVPHFacetedHydroBaseInst.cc \
	SVPHCorrectionsPolicyInst.cc \
	SVPHMassDensityPolicyInst.cc \
	CompatibleFaceSpecificThermalEnergyPolicyInst.cc \
	MeshIdealHPolicyInst.cc \
	SpecificThermalEnergyVolumePolicyInst.cc

PYTHONTARGETS = \
	SVPHHydros.py

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
