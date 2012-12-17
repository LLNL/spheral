# -*-makefile-*-
# SolidSpheral++ Strength package makefile
#-------------------------------------------------------------------------------
PKGNAME = Strength
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	SolidFieldNames.cc \
	SolidNodeListInst.cc \
	DeviatoricStressPolicyInst.cc \
	PlasticStrainPolicyInst.cc \
	BulkModulusPolicyInst.cc \
	ShearModulusPolicyInst.cc \
	YieldStrengthPolicyInst.cc \
	StrengthSoundSpeedPolicyInst.cc

PYTHONTARGETS = \
	SolidNodeLists.py

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SOLIDTOP)/helpers/makefile_master
