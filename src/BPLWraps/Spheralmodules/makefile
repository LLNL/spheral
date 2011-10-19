# -*-makefile-*-
# Spheral++ Spheralmodules package makefile
#-------------------------------------------------------------------------------
PKGNAME = Spheralmodules
PKGDIR = $(PKGNAME)
BPLMODTARGET = $(PKGNAME).$(SHLIBEXT)
BPLTARGETS = \
	Spheralmodules.cpp
LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
EXTRALIBS = $(patsubst %, $(LIBDIR)/%.so, $(FILTEREDBPLPKGS))
