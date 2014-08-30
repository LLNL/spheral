# -*-makefile-*-
# Spheral++ Kernel package makefile
#-------------------------------------------------------------------------------
PKGNAME = Kernel
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	BSplineKernel.cc \
	W4SplineKernel.cc \
	GaussianKernel.cc \
	SuperGaussianKernel.cc \
	PiGaussianKernel.cc \
	HatKernel.cc \
	SincKernel.cc \
	NSincPolynomialKernelInst.cc \
	NBSplineKernelInst.cc \
	QuarticSplineKernel.cc \
	QuinticSplineKernel.cc \
	TableKernelInst.cc \
	VolumeIntegrationFunctionsInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

EXTRAFLAGS =

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
