# -*-makefile-*-
# Spheral++ Field package makefile
#-------------------------------------------------------------------------------
PKGNAME = Field
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	NodeIteratorBase.cc \
	AllNodeIterator.cc \
	InternalNodeIterator.cc \
	GhostNodeIterator.cc \
	MasterNodeIterator.cc \
	CoarseNodeIterator.cc \
	RefineNodeIterator.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
