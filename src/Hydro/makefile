# -*-makefile-*-
# Spheral++ Hydro package makefile
#-------------------------------------------------------------------------------
PKGNAME = Hydro
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	HydroFieldNames.cc \
	PressurePolicyInst.cc \
	SoundSpeedPolicyInst.cc \
	SpecificThermalEnergyPolicyInst.cc \
	NonSymmetricSpecificThermalEnergyPolicyInst.cc \
	PositionPolicyInst.cc \
	VolumePolicyInst.cc \
	VoronoiMassDensityPolicyInst.cc \
	SumVoronoiMassDensityPolicyInst.cc \
	SecondMomentHourglassControlInst.cc \
	ThirdMomentHourglassControlInst.cc \
	VoronoiHourglassControlInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
