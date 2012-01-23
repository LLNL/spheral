# -*-makefile-*-
# Spheral++ CXXTypes BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = CXXTypes
PKGDIR = $(PKGNAME)/
BPLMODTARGET = $(PKGNAME).$(SHLIBEXT)
BPLTARGETS = \
	STDStringWrap.cpp \
	STLVectorPrimitivesWrap.cpp \
	STLVectorFieldWrap.cpp \
	STLVectorFieldListWrap.cpp \
	STLVectorBoundaryWrap.cpp \
	STLVectorPhysicsWrap.cpp \
	STLVectorDistributedWrap.cpp \
	PairWraps.cpp \
	CXXTypes.cpp
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
