# -*-makefile-*-
# Spheral++ DataBase package makefile
#-------------------------------------------------------------------------------
PKGNAME = DataBase
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	StateBase.cc \
	State.cc \
	StateDerivatives.cc \
	DataBase.cc \
	IncrementState.cc \
	IncrementBoundedState.cc \
	ReplaceState.cc \
	ReplaceBoundedState.cc \
	CopyState.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
##EXTRALIBS = -L$(SRCTOP)/lib -lField -lNodeList -lGeometry -lField -lNodeList -lGeometry

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
