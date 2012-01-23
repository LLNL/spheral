# -*-makefile-*-
# Spheral++ Neighbor BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Neighbor
PKGDIR = $(PKGNAME)
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	Neighbor.pyste \
	GridCellIndex.pyste \
	GridCellPlane.pyste \
	NestedGridNeighbor.pyste \
	ConnectivityMap.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
