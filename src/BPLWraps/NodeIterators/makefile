# -*-makefile-*-
# Spheral++ NodeIterators BPL wrap package makefile
#-------------------------------------------------------------------------------
PKGNAME = NodeIterators
PKGDIR = $(PKGNAME)/
MODTARGET = $(PKGNAME).$(SHLIBEXT)
PYSTETARGETS = \
	AllNodeIterator.pyste \
	InternalNodeIterator.pyste \
	GhostNodeIterator.pyste \
	MasterNodeIterator.pyste \
	CoarseNodeIterator.pyste \
	RefineNodeIterator.pyste

LINKLIST = linklist
SPHERALTOP = ../..

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I.
