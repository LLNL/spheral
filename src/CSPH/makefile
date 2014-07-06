# -*-makefile-*-
# Spheral++ CSPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = CSPH
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeCSPHSumMassDensityInst.cc \
	computeCSPHCorrectionsInst.cc \
	computeHullVolumesInst.cc \
	CSPHHydroBaseInst.cc

PYTHONTARGETS = \
	CSPHHydros.py

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
