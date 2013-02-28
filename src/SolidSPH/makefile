# -*-makefile-*-
# SolidSpheral++ SolidSPH package makefile
#-------------------------------------------------------------------------------
PKGNAME = SolidSPH
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	SolidSPHHydroBaseInst.cc

PYTHONTARGETS = \
	SolidSPHHydros.py

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SOLIDTOP)/helpers/makefile_master
