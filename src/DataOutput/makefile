# -*-makefile-*-
# Spheral++ DataOutput package makefile
#-------------------------------------------------------------------------------
PKGNAME = DataOutput
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	RestartRegistrar.cc \
	RestartableObject.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAINCLUDES =
##EXTRALIBS = 
##EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
