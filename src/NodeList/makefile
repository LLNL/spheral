# -*-makefile-*-
# Spheral++ NodeList package makefile
#-------------------------------------------------------------------------------
PKGNAME = NodeList
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	NodeListRegistrar.cc \
	NodeList.cc \
	FluidNodeList.cc \
	SmoothingScaleBase.cc \
	FixedSmoothingScale.cc \
	SPHSmoothingScale.cc \
	ASPHSmoothingScale.cc \
	nthNodalMoment.cc \
	generateVoidNodes.cc

PYTHONTARGETS = \
	FluidNodeLists.py \
	VoidNodeLists.py

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#EXTRAFLAGS = -DDEBUG

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

#-------------------------------------------------------------------------------
TestSphNodeList:  tests/TestSphNodeList.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -I/home/owen/Spheral++ -o tests/TestSphNodeList tests/TestSphNodeList.cc ${SRCTARGETS:.cc=.o} /home/owen/Spheral++/Field/FieldBase.o /home/owen/Spheral++/Field/Field.o $(LIBS)
