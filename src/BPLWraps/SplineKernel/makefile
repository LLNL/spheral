# -*-makefile-*-
# Spheral++ SplineKernel BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = SplineKernel
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	HatKernel.pyste \
	BSplineKernel.pyste \
	W4SplineKernel.pyste \
	QuarticSplineKernel.pyste \
	QuinticSplineKernel.pyste \
	NBSplineKernel.pyste
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
