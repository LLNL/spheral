# -*-makefile-*-
# Spheral++ ExternalForce package makefile
#-------------------------------------------------------------------------------
PKGNAME = ExternalForce
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	PointPotentialInst.cc \
	NFWPotentialInst.cc \
	ConstantAccelerationInst.cc \
	LinearAccelerationInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
