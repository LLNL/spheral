# -*-makefile-*-
# Spheral++ ArtificialViscosity package makefile
#-------------------------------------------------------------------------------
PKGNAME = ArtificialViscosity
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	ArtificialViscosityInst.cc \
	MonaghanGingoldViscosityInst.cc \
	MorrisMonaghanReducingViscosityInst.cc \
	TensorMonaghanGingoldViscosityInst.cc \
	FiniteVolumeViscosityInst.cc \
	TensorSVPHViscosityInst.cc \
	VonNeumanViscosityInst.cc
#	MonaghanGingoldViscosityRZ.cc \
#	GradPressureViscosity.cc \
#	ArtificialViscosityList.cc \
#	MASHVonNeumanViscosity.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
