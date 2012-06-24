# -*-makefile-*-
# Spheral++ CSPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = CSPH
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeCSPHSumMassDensity.cc \
	computeCSPHCorrections.cc \
	CSPHHydroBase.cc

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
