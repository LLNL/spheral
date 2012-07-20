# -*-makefile-*-
# SolidSpheral++ Damage package makefile
#-------------------------------------------------------------------------------
PKGNAME = Damage
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	DamageModel.cc \
	TensorDamageModel.cc \
	StrainPolicy.cc \
	TensorStrainPolicy.cc \
	TensorDamagePolicy.cc \
	YoungsModulusPolicy.cc \
	LongitudinalSoundSpeedPolicy.cc \
	weibullFlawDistribution.cc \
	computeFragmentField.cc

PYTHONTARGETS = \
	identifyFragments.py \
	WeibullDamage.py \
	GradyKippTensorDamage.py

SPHERALTOP = ../../../spheral/src
SOLIDTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SOLIDTOP)/helpers/makefile_master
