# -*-makefile-*-
# Spheral++ NodeGenerators makefile
#-------------------------------------------------------------------------------
PKGNAME = NodeGenerators
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)

SRCTARGETS = \
	generateCylDistributionFromRZ.cc \
	fillFacetedVolume.cc

PYTHONTARGETS = \
	NodeGeneratorBase.py \
	GenerateNodeDistribution2d.py \
	GenerateNodeDistribution3d.py \
	CompositeNodeDistribution.py \
	distributeNodesGeneric.py \
	DistributeNodes.py \
	ParMETISDistributeNodes.py \
	MortonOrderDistributeNodes.py \
	PeanoHilbertDistributeNodes.py \
	VoronoiDistributeNodes.py \
	NestedGridDistributeNodes.py \
	GzipFileNodeGenerator.py \
	CubicNodeGenerator.py \
	GenerateSphericalShellSection.py \
	AbaqusNodeGenerator.py \
	CentroidalVoronoiRelaxation.py

SPHERALTOP = ..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

