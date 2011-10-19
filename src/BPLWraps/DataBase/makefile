# -*-makefile-*-
# Spheral++ DataBase BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = DataBase
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	StateBase.pyste \
	State.pyste \
	StateDerivatives.pyste \
	DataBase.pyste
BPLTARGETS = \
	ExtendDataBase.cpp
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
