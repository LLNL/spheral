# -*-makefile-*-
# Spheral++ Hydro package makefile
#-------------------------------------------------------------------------------
PKGNAME = Hydro
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	HydroFieldNames.cc \
	PressurePolicy.cc \
	SoundSpeedPolicy.cc \
	SpecificThermalEnergyPolicy.cc \
	PositionPolicy.cc \
	VolumePolicy.cc \
	VoronoiMassDensityPolicy.cc \
	SumVoronoiMassDensityPolicy.cc \
	SecondMomentHourglassControl.cc \
	ThirdMomentHourglassControl.cc \
	VoronoiHourglassControl.cc

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
#	$(CXX) $(PURIFY) $(CXXFLAGS) $(EXTRAFLAGS) -Bstatic -o tests/TestSod tests/TestSod.cc -L$(SRCTOP)/lib $(PKGLIBS) $(PKGLIBS) $(PKGLIBS) -lm
