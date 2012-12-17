# -*-makefile-*-
# Spheral++ CSPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = CSPH
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeCSPHSumMassDensityInst.cc \
	computeCSPHCorrectionsInst.cc \
	CSPHHydroBaseInst.cc

PYTHONTARGETS = \
	CSPHHydros.py

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
