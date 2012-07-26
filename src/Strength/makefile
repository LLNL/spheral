# -*-makefile-*-
# SolidSpheral++ Strength package makefile
#-------------------------------------------------------------------------------
PKGNAME = Strength
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	SolidFieldNames.cc \
	SolidNodeList.cc \
	DeviatoricStressPolicy.cc \
	PlasticStrainPolicy.cc \
	BulkModulusPolicy.cc \
	ShearModulusPolicy.cc \
	YieldStrengthPolicy.cc \
	StrengthSoundSpeedPolicy.cc

PYTHONTARGETS = \
	SolidNodeLists.py

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SOLIDTOP)/helpers/makefile_master
