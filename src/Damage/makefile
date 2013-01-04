# -*-makefile-*-
# SolidSpheral++ Damage package makefile
#-------------------------------------------------------------------------------
PKGNAME = Damage
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	DamageModelInst.cc \
	TensorDamageModelInst.cc \
	StrainPolicyInst.cc \
	TensorStrainPolicyInst.cc \
	TensorDamagePolicyInst.cc \
	YoungsModulusPolicyInst.cc \
	LongitudinalSoundSpeedPolicyInst.cc \
	weibullFlawDistributionInst.cc \
	computeFragmentFieldInst.cc

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
