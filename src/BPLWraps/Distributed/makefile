# -*-makefile-*-
# Spheral++ Distributed BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Distributed
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	DistributedBoundary.pyste \
	NestedGridDistributedBoundary.pyste \
	DomainNode.pyste \
	RedistributeNodes.pyste \
	DistributeByXPosition.pyste \
	NestedGridRedistributeNodes.pyste \
	SpaceFillingCurveRedistributeNodes.pyste \
	MortonOrderRedistributeNodes.pyste \
	PeanoHilbertOrderRedistributeNodes.pyste \
	SortAndDivideRedistributeNodes.pyste \
	$(PARMETISBPLTARGET)

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
EXTRAFLAGS = $(PARMETISINCS)
EXTRALIBS = $(PARMETISLIBS)
