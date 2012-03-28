# -*-makefile-*-
# Spheral++ Kernel BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Kernel
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	GaussianKernel.pyste \
	SuperGaussianKernel.pyste \
	PiGaussianKernel.pyste \
	SincKernel.pyste \
	NSincPolynomialKernel.pyste \
	TableKernel.pyste
BPLTARGETS = \
	ExtendTableKernel.cpp
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
