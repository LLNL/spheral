# -*-makefile-*-
# Spheral++ FieldOperations BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = FieldOperations
PKGDIR = $(PKGNAME)/
BPLMODTARGET = $(PKGNAME).$(SHLIBEXT)
BPLTARGETS = \
	FieldOperations.cpp
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I. -I$(SPHERALTOP)/BPLWraps
