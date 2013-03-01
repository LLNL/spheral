# -*-makefile-*-
# Spheral++ FieldOperations package makefile
#-------------------------------------------------------------------------------
PKGNAME = FieldOperations
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	smoothFieldsInst.cc \
	gradientInst.cc \
	divergenceInst.cc \
	smoothFieldsMashInst.cc \
	smoothFieldsMash2Inst.cc \
	sampleFieldsMashInst.cc \
	splatFieldsMashInst.cc \
	gradientMashInst.cc \
	divergenceMashInst.cc \
	gradientPairWiseInst.cc \
	divergencePairWiseInst.cc \
	sampleMultipleFieldsMashInst.cc \
	splatMultipleFieldsMashInst.cc \
	gradDivVectorFieldListInst.cc \
	gradDivVectorFieldListSimpleInst.cc \
	gradDivVectorFieldListGoldenInst.cc \
	gradDivVectorFieldListGolden2Inst.cc \
	gradDivVectorFieldListMashInst.cc \
	gradDivVectorFieldListPairWiseInst.cc \
	sampleMultipleFields2LatticeInst.cc \
	sampleMultipleFields2LatticeMashInst.cc \
	binFieldList2LatticeInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
