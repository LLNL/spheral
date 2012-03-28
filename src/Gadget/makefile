# -*-makefile-*-
# Spheral++ Gadget package makefile
#-------------------------------------------------------------------------------
PKGNAME = Gadget
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
  FakeGadgetICFile.cc \
	GadgetGravityForce.cc 

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

include $(SPHERALTOP)/helpers/makefile_master

EXTRAFLAGS = -I../Physics -I../Utilities -I../DataBase -I../Field -I$(GADGET_PATH)/S-Gadget 
GADGET_OBJECT_PATH = $(GADGET_PATH)/S-Gadget
EXTRALIBS = $(GADGET_OBJECT_PATH)/accel.o $(GADGET_OBJECT_PATH)/allocate.o $(GADGET_OBJECT_PATH)/allvars.o $(GADGET_OBJECT_PATH)/begrun.o $(GADGET_OBJECT_PATH)/cooling.o $(GADGET_OBJECT_PATH)/density.o $(GADGET_OBJECT_PATH)/endrun.o $(GADGET_OBJECT_PATH)/ewald.o $(GADGET_OBJECT_PATH)/forcetree.o $(GADGET_OBJECT_PATH)/global.o $(GADGET_OBJECT_PATH)/gravgrape.o $(GADGET_OBJECT_PATH)/gravtree.o $(GADGET_OBJECT_PATH)/hydra.o $(GADGET_OBJECT_PATH)/init.o $(GADGET_OBJECT_PATH)/io.o $(GADGET_OBJECT_PATH)/main.o $(GADGET_OBJECT_PATH)/potential.o $(GADGET_OBJECT_PATH)/potgrape.o $(GADGET_OBJECT_PATH)/predict.o $(GADGET_OBJECT_PATH)/read_ic.o $(GADGET_OBJECT_PATH)/read_ic_cluster.o $(GADGET_OBJECT_PATH)/restart.o $(GADGET_OBJECT_PATH)/run.o $(GADGET_OBJECT_PATH)/system.o $(GADGET_OBJECT_PATH)/timeline.o $(GADGET_OBJECT_PATH)/timestep.o $(GADGET_OBJECT_PATH)/veldisp.o 
# EXTRAFLAGS = -DVERBOSEDT
#-------------------------------------------------------------------------------
