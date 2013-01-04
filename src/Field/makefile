# -*-makefile-*-
# Spheral++ Field package makefile
#-------------------------------------------------------------------------------
PKGNAME = Field
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	NodeIteratorBaseInst.cc \
	AllNodeIteratorInst.cc \
	InternalNodeIteratorInst.cc \
	GhostNodeIteratorInst.cc \
	MasterNodeIteratorInst.cc \
	CoarseNodeIteratorInst.cc \
	RefineNodeIteratorInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
