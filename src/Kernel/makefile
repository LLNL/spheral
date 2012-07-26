# -*-makefile-*-
# Spheral++ Kernel package makefile
#-------------------------------------------------------------------------------
PKGNAME = Kernel
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	BSplineKernel.cc \
	W4SplineKernel.cc \
	GaussianKernel.cc \
	SuperGaussianKernel.cc \
	PiGaussianKernel.cc \
	HatKernel.cc \
	SincKernel.cc \
	NSincPolynomialKernel.cc \
	NBSplineKernel.cc \
	QuarticSplineKernel.cc \
	QuinticSplineKernel.cc \
	TableKernel.cc \
	VolumeIntegrationFunctions.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

EXTRAFLAGS =

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
