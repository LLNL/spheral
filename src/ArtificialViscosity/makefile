# -*-makefile-*-
# Spheral++ ArtificialViscosity package makefile
#-------------------------------------------------------------------------------
PKGNAME = ArtificialViscosity
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	ArtificialViscosity.cc \
	MonaghanGingoldViscosity.cc \
	TMGViscosity.cc \
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
include $(SPHERALTOP)/helpers/makefile_taumaster
