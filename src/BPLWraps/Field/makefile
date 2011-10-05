# -*-makefile-*-
# Spheral++ Field BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = Field
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	FieldBase.pyste \
	Field1d.pyste \
	Field2d.pyste \
	Field3d.pyste \
	FieldList.pyste \
	FieldListSet.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
