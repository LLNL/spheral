# -*-makefile-*-
# Spheral++ Integrator package makefile
#-------------------------------------------------------------------------------
PKGNAME = Integrator
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	IntegratorInst.cc \
	PredictorCorrectorInst.cc \
	SynchronousRK1Inst.cc \
	SynchronousRK2Inst.cc \
	SynchronousRK4Inst.cc \
	CheapSynchronousRK2Inst.cc \

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_master
