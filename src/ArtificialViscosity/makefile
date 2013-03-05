# -*-makefile-*-
# Spheral++ ArtificialViscosity package makefile
#-------------------------------------------------------------------------------
PKGNAME = ArtificialViscosity
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	ArtificialViscosityInst.cc \
	MonaghanGingoldViscosityInst.cc \
	TensorMonaghanGingoldViscosityInst.cc \
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
