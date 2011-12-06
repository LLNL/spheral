# -*-makefile-*-
# Spheral++ Hydro BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Hydro
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	Hydro.pyste \
	SecondMomentHourglassControl.pyste \
	ThirdMomentHourglassControl.pyste
#	IdealHUpdate.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I. -DBOOST_PYTHON_MAX_ARITY=46
