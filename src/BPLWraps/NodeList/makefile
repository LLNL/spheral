# -*-makefile-*-
# Spheral++ NodeList BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = NodeList
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	NodeListRegistrar.pyste \
	NodeList.pyste \
	FluidNodeList.pyste \
	SmoothingScales.pyste \
	FluidDerivatives.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I. -DBOOST_PYTHON_MAX_ARITY=46
