# -*-makefile-*-
# Spheral++ Utilities package makefile
#-------------------------------------------------------------------------------
PKGNAME = Utilities
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	DBC.cc \
	Process.cc \
	iterateIdealHInst.cc \
	KeyTraits.cc \
	mortonOrderIndiciesInst.cc \
	peanoHilbertOrderIndiciesInst.cc \
	PeanoHilbertTransform.cc \
	nodeOrderingInst.cc \
	RegisterMPIDataTypes.cc \
	RedistributionRegistrar.cc \
	nodeBoundingBoxesInst.cc \
	boundingBoxInst.cc \
	globalBoundingVolumesInst.cc \
	lineSegmentIntersections.cc \
	integrateThroughMeshAlongSegmentInst.cc \
	numberDensityInst.cc \
	medianPositionInst.cc \
	coarsenBinnedValuesInst.cc \
	segmentIntersectEdges.cc \
	pointOnPolygon.cc \
	pointInPolygon.cc \
	pointOnPolyhedron.cc \
	pointInPolyhedron.cc

PYTHONTARGETS = \
	fitspline.py

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#EXTRAFLAGS = 
#EXTRALIBS = $(PAPILIBS) 
#EXTRALIBS = -lgsl -lgslcblas

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
