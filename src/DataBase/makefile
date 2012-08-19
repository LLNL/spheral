# -*-makefile-*-
# Spheral++ DataBase package makefile
#-------------------------------------------------------------------------------
PKGNAME = DataBase
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	StateBaseInst.cc \
	StateInst.cc \
	StateDerivativesInst.cc \
	DataBaseInst.cc \
	IncrementStateInst.cc \
	IncrementBoundedStateInst.cc \
	ReplaceStateInst.cc \
	ReplaceBoundedStateInst.cc \
	CopyStateInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
##EXTRALIBS = -L$(SRCTOP)/lib -lField -lNodeList -lGeometry -lField -lNodeList -lGeometry

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
