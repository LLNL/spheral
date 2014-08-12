# -*-makefile-*-
# Spheral++ TaylorSPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = TaylorSPH
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	computeTaylorSPHCorrectionsInst.cc \
	TaylorSPHHydroBaseInst.cc

PYTHONTARGETS = \
	TaylorSPHHydros.py

LINKLIST = linklist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
