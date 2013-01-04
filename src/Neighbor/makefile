# -*-makefile-*-
# Spheral++ Neighbor package makefile
#-------------------------------------------------------------------------------
PKGNAME = Neighbor
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	GridCellIndexInst.cc \
	GridCellPlaneInst.cc \
	NeighborInst.cc \
	NestedGridNeighborInst.cc \
	TreeNeighborInst.cc \
	ConnectivityMapInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#EXTRAFLAGS = -DDEBUG

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
