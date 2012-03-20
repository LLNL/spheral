# -*-makefile-*-
# Spheral++ Utilities BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Utilities
PKGDIR = $(PKGNAME)/
BPLMODTARGET = $(PKGNAME).$(SHLIBEXT)
BPLTARGETS = \
	initializeTau.cpp \
	erff.cpp \
	newtonRaphson.cpp \
	globalNodeIDs.cpp \
	rotationMatrix.cpp \
	iterateIdealH.cpp \
	CSPHUtilities.cpp \
	visitWriter.cpp \
	mortonOrderIndicies.cpp \
	peanoHilbertOrderIndicies.cpp \
	nodeOrdering.cpp \
	packElement.cpp \
	Utilities.cpp

# PYSTETARGETS = \
# 	initializeTau.pyste \
# 	erff.pyste \
# 	newtonRaphson.pyste \
# 	globalNodeIDs.pyste \
# 	rotationMatrix.pyste \
# 	iterateIdealH.pyste \
# 	CSPHUtilities.pyste \
# 	visitWriter.pyste \
# 	mortonOrderIndicies.pyste

SPHERALTOP = ../..
EXTRALIBS = -lvisit_writer

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

EXTRAINCLUDES = $(TAU_INCLUDE)

