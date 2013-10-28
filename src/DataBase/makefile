# -*-makefile-*-
# Spheral++ DataBase package makefile
#-------------------------------------------------------------------------------
PKGNAME = DataBase
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	StateBaseInst.cc \
	StateInst.cc \
	StateDerivativesInst.cc \
	DataBaseInst.cc \
	IncrementStateInst.cc \
	IncrementBoundedStateInst.cc \
	ReplaceStateInst.cc \
	ReplaceBoundedStateInst.cc \
	CopyStateInst.cc \
	IncrementFieldListInst.cc \
	IncrementBoundedFieldListInst.cc \
	ReplaceFieldListInst.cc \
	ReplaceBoundedFieldListInst.cc \
	CopyFieldListInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
##EXTRALIBS = -L$(SRCTOP)/lib -lField -lNodeList -lGeometry -lField -lNodeList -lGeometry

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
