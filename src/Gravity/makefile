# -*-makefile-*-
# Spheral++ Gadget package makefile
#-------------------------------------------------------------------------------
PKGNAME = Gravity
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	NBodyGravity.cc \
#	SPHGravity.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

EXTRAFLAGS = -I../Physics -I../Utilities -I../DataBase -I../Field 
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
