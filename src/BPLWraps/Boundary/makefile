# -*-makefile-*-
# Spheral++ Boundary BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Boundary
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)

PYSTETARGETS = \
	Boundary.pyste \
	PlanarBoundary.pyste \
	ReflectingBoundary.pyste \
	RigidBoundary.pyste \
	PeriodicBoundary.pyste \
	ConstantVelocityBoundary.pyste \
	ConstantXVelocityBoundary.pyste \
	ConstantYVelocityBoundary.pyste \
	ConstantZVelocityBoundary.pyste \
	ConstantRVelocityBoundary.pyste \
	ConstantBoundary.pyste \
	SphericalBoundary.pyste \
	CylindricalBoundary.pyste \
	AxialSymmetryBoundary.pyste

BPLTARGETS = \
	ExtendBoundary.cpp

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
