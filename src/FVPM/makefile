# -*-makefile-*-
# Spheral++ MHD package makefile
#-------------------------------------------------------------------------------
PKGNAME = FVPM
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	TestFunction.cc \
	QuadRule.cc \
	CircularQuadRule.cc \
	GaussianCircularQuadRule.cc \
	gauss_legendre.cc \
	FVPMFluidDerivatives.cc 

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

EXTRALIBS =

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

#TestSod:
#	$(CXX) $(PURIFY) $(CXXFLAGS) $(EXTRAFLAGS) -Bstatic -o tests/TestSod tests/TestSod.cc -L$(SRCTOP)/lib $(PKGLIBS) $(PKGLIBS) $(PKGLIBS) -lm
