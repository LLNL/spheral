# -*-makefile-*-
# Spheral++ ArtificialViscosity BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = ArtificialViscosity
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	ArtificialViscosity.pyste \
	MonaghanGingoldViscosity.pyste \
	MonaghanGingoldViscosityRZ.pyste \
	TensorMonaghanGingoldViscosity.pyste \

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I. -DBOOST_PYTHON_MAX_ARITY=25

