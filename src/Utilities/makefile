# -*-makefile-*-
# Spheral++ Utilities package makefile
#-------------------------------------------------------------------------------
PKGNAME = Utilities
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	DBC.cc \
	Process.cc \
	iterateIdealH.cc \
	KeyTraits.cc \
	mortonOrderIndicies.cc \
	peanoHilbertOrderIndicies.cc \
	PeanoHilbertTransform.cc \
	nodeOrdering.cc \
	RegisterMPIDataTypes.cc \
	RedistributionRegistrar.cc \
	spheralWildMagicConverters.cc \
	orientedBoundingBox.cc \
	nodeBoundingBoxes.cc \
	boundingBox.cc \
	globalBoundingVolumes.cc \
	lineSegmentIntersections.cc \
	pointDistances.cc \
	integrateThroughMeshAlongSegment.cc \
	numberDensity.cc \
	medianPosition.cc \
	coarsenBinnedValues.cc \
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
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
