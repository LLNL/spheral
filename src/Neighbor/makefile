# -*-makefile-*-
# Spheral++ Neighbor package makefile
#-------------------------------------------------------------------------------
PKGNAME = Neighbor
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	GridCellIndex.cc \
	GridCellPlane.cc \
	Neighbor.cc \
	NestedGridNeighbor.cc \
	ConnectivityMap.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#EXTRAFLAGS = -DDEBUG

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
