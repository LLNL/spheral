# -*-makefile-*-
# Spheral++ FileIO BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = FileIO
PKGDIR = $(PKGNAME)
MODTARGET = $(PKGNAME).$(SHLIBEXT)

PYSTETARGETS = \
	FileIOEnums.pyste \
	FileIO.pyste \
	FlatFileIO.pyste \
	PyFileIO.pyste \
	vectorstringUtilities.pyste \
	fieldstringUtilities.pyste

BPLTARGETS = \
	ExtendFileIO.cpp

SRCTARGETS = \
	vectorstringUtilitiesSpecializations.cc

PYTHONTARGETS = \
	ExtendFlatFileIO.py

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
