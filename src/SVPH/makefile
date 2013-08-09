# -*-makefile-*-
# Spheral++ SVPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = SVPH
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeSVPHCorrectionsInst.cc \
	sampleFieldListSVPHInst.cc \
	gradientFieldListSVPHInst.cc \
	SVPHHydroBaseInst.cc \
	SVPHCorrectionsPolicyInst.cc

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
