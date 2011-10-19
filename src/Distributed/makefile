# -*-makefile-*-
# Spheral++ Distributed package makefile
#-------------------------------------------------------------------------------
PKGNAME = Distributed
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	DistributedBoundary.cc \
	NestedGridDistributedBoundary.cc \
	BoundingVolumeDistributedBoundary.cc \
	RedistributeNodes.cc \
	DistributeByXPosition.cc \
	NestedGridRedistributeNodes.cc \
	SortAndDivideRedistributeNodes.cc \
	SortAndDivideRedistributeNodes1d.cc \
	SortAndDivideRedistributeNodes2d.cc \
	SortAndDivideRedistributeNodes3d.cc \
	SpaceFillingCurveRedistributeNodes.cc \
	MortonOrderRedistributeNodes.cc \
	PeanoHilbertOrderRedistributeNodes.cc \
	VoronoiRedistributeNodes.cc \
	waitAllWithDeadlockDetection.cc \
	$(PARMETISTARGET)

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

EXTRAFLAGS += $(PARMETISINCS)
EXTRALIBS += $(PARMETISLIBS)

#-------------------------------------------------------------------------------
# Override the default PYFFLEOPT flag, because using g++ with -g on the OSF Dec
# machines causes an internal compiler error!
OPTFLAGS = ${DISTRIBUTEDOPT}
PYFFLEOPT = ${DISTRIBUTEDOPT}

