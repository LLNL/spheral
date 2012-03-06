# -*-makefile-*-
# Spheral++ NodeGenerators BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = NodeGenerators
PKGDIR = $(PKGNAME)/
BPLMODTARGET = $(PKGNAME).$(SHLIBEXT)

BPLTARGETS = \
	NodeGenerators.cpp

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
	NestedGridDistributeNodes.py \
	GzipFileNodeGenerator.py \
	CubicNodeGenerator.py \
	GenerateSphericalShellSection.py

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
