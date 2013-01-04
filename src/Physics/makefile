# -*-makefile-*-
# Spheral++ Physics package makefile
#-------------------------------------------------------------------------------
PKGNAME = Physics
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	PhysicsInst.cc \
	GenericHydroInst.cc \
	GenericBodyForceInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
