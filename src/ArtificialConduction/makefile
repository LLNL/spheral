# -*-makefile-*-
# Spheral++ Artificial Conduction package makefile
#-------------------------------------------------------------------------------
PKGNAME = ArtificialConduction
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	ArtificialConductionInst.cc \
	ArtificialConductionPolicyInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# EXTRAFLAGS = -DDEBUG
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
