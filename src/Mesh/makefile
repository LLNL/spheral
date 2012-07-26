# -*-makefile-*-
# Spheral++ Mesh package makefile
#-------------------------------------------------------------------------------
PKGNAME = Mesh
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	Cell.cc \
	VoroPP2d.cc \
	VoroPP.cc \
	MeshConstructionUtilities.cc \
	findMatchingVertex.cc \
	LineMesh.cc \
	LineEdge.cc \
	LineFace.cc \
	LineZone.cc \
	PolygonalMesh.cc \
	PolygonalEdge.cc \
	PolygonalFace.cc \
	PolygonalZone.cc \
	PolyhedralMesh.cc \
	PolyhedralEdge.cc \
	PolyhedralFace.cc \
	PolyhedralZone.cc \
	computeGenerators.cc \
	generateMesh.cc \
	MeshPolicy.cc

PYTHONTARGETS = \
	generateMesh.py \
	siloMeshDump.py

EXTRALIBS = -lvoro++_2d

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
