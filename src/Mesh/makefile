# -*-makefile-*-
# Spheral++ Mesh package makefile
#-------------------------------------------------------------------------------
PKGNAME = Mesh
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	MeshConstructionUtilities.cc \
	findMatchingVertexInst.cc \
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
	computeGeneratorsInst.cc \
	generateMeshInst.cc \
	MeshPolicyInst.cc \
	computeElementOwner.cc

PYTHONTARGETS = \
	generateMesh.py \
	siloMeshDump.py

#EXTRAINCLUDES = -I $(SPHERALDIR)/include/polytope
#EXTRALIBS = -l polytope

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
