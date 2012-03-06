# -*-makefile-*-
# Spheral++ Geometry BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Geometry
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	GeomVector.pyste \
	Geom3Vector.pyste \
	GeomTensor.pyste \
	GeomSymmetricTensor.pyste \
	EigenStruct.pyste \
	GeomPlane.pyste \
	GeomThirdRankTensor.pyste
#	GeomTensorFunctions.pyste

BPLTARGETS = \
	ExtendGeomVector.cpp \
	ExtendGeom3Vector.cpp \
	ExtendGeomTensor.cpp \
	ExtendGeomThirdRankTensor.cpp \
	ExtendJacobi.cpp

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
