# -*-makefile-*-
# Spheral++ MHD package makefile
#-------------------------------------------------------------------------------
PKGNAME = MHD
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	MHDFieldNames.cc \
	ConductingFluidNodeList.cc \
	CurrentDensityUpdatePolicy.cc \
	SpitzerResistivityUpdatePolicy.cc \
	MagnetosonicSpeedUpdatePolicy.cc \
	PriceMonaghanDissipation.cc \
	MHD.cc 

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

EXTRALIBS = -L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscksp -lpetscmat -lpetscvec -lpetsc

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
EXTRAINCLUDES = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) \
	-I$(SPHERALTOP)/../../include/pyMPI2.4 -I$(SPHERALTOP)/../../include/python$(PYTHONVERSION)

#TestSod:
#	$(CXX) $(PURIFY) $(CXXFLAGS) $(EXTRAFLAGS) -Bstatic -o tests/TestSod tests/TestSod.cc -L$(SRCTOP)/lib $(PKGLIBS) $(PKGLIBS) $(PKGLIBS) -lm
