# -*-makefile-*-
# Spheral++ Integrator package makefile
#-------------------------------------------------------------------------------
PKGNAME = Integrator
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	Integrator.cc \
	PredictorCorrector.cc \
	SynchronousRK1.cc \
	SynchronousRK2.cc \
	CheapSynchronousRK2.cc \

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
