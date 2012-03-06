# -*-makefile-*-
# Spheral++ FieldOperations package makefile
#-------------------------------------------------------------------------------
PKGNAME = FieldOperations
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	smoothFields.cc \
	gradient.cc \
	divergence.cc \
	smoothFieldsMash.cc \
	smoothFieldsMash2.cc \
	sampleFieldsMash.cc \
	splatFieldsMash.cc \
	gradientMash.cc \
	divergenceMash.cc \
	gradientPairWise.cc \
	divergencePairWise.cc \
	sampleMultipleFieldsMash.cc \
	splatMultipleFieldsMash.cc \
	gradDivVectorFieldList.cc \
	gradDivVectorFieldListSimple.cc \
	gradDivVectorFieldListGolden.cc \
	gradDivVectorFieldListGolden2.cc \
	gradDivVectorFieldListMash.cc \
	gradDivVectorFieldListPairWise.cc \
	sampleMultipleFields2Lattice.cc \
	sampleMultipleFields2LatticeMash.cc \
	binFieldList2Lattice.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
