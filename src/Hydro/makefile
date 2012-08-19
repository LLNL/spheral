# -*-makefile-*-
# Spheral++ Hydro package makefile
#-------------------------------------------------------------------------------
PKGNAME = Hydro
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	HydroFieldNames.cc \
	PressurePolicyInst.cc \
	SoundSpeedPolicyInst.cc \
	SpecificThermalEnergyPolicyInst.cc \
	PositionPolicyInst.cc \
	VolumePolicyInst.cc \
	VoronoiMassDensityPolicyInst.cc \
	SumVoronoiMassDensityPolicyInst.cc \
	SecondMomentHourglassControlInst.cc \
	ThirdMomentHourglassControlInst.cc \
	VoronoiHourglassControlInst.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

##-------------------------------------------------------------------------------
#CXXPKGS = \
#	Geometry \
#	NodeList \
#	Field \
#	Kernel \
#	Material \
#	Neighbor \
#	DataBase \
#	Boundary \
#	ArtificialViscosity \
#	Physics \
#	Hydro \
#	Integrator
#PKGLIBS := $(patsubst %, -l%, $(CXXPKGS))
#TestSod:
#	$(CXX) $(PURIFY) $(CXXFLAGS) $(EXTRAFLAGS) -Bstatic -o tests/TestSod tests/TestSodInst.cc -L$(SRCTOP)/lib $(PKGLIBS) $(PKGLIBS) $(PKGLIBS) -lm
