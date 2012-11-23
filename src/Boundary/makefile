# -*-makefile-*-
# Spheral++ Boundary package makefile
#-------------------------------------------------------------------------------
PKGNAME = Boundary
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	BoundaryInst.cc \
	PlanarBoundaryInst.cc \
	ReflectingBoundaryInst.cc \
	RigidBoundaryInst.cc \
	PeriodicBoundaryInst.cc \
	ConstantVelocityBoundaryInst.cc \
	ConstantXVelocityBoundaryInst.cc \
	ConstantYVelocityBoundaryInst.cc \
	ConstantZVelocityBoundaryInst.cc \
	ConstantRVelocityBoundaryInst.cc \
	ConstantBoundaryInst.cc \
	SphericalBoundary.cc \
	CylindricalBoundary.cc \
	AxialSymmetryBoundary.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

#EXTRAFLAGS = -DDEBUG

#------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

#-------------------------------------------------------------------------------
# PYTHONLIBS = -L/usr/local/lib/python1.5/config -L/usr/local/lib/python1.5/lib-dynload `find /usr/local/lib/python1.5/lib-dynload -name '*.so'` -lpython1.5
# SPHERALLIBS = -L$(SRCTOP)/lib -lPyffle -lBoundary -lField -lGeometry -lKernel -lMaterial -lNeighbor -lNodeList

# TestReflectingBoundary:  tests/TestReflectingBoundary.cc $(HDRS)
# 	$(CXX) $(PURIFY) $(CXXFLAGS) -I/home/owen/Spheral++ -o tests/TestReflectingBoundary tests/TestReflectingBoundary.cc $(LIBS) -lm

# 	#$(CXX) $(CXXFLAGS) -I/home/owen/Spheral++ -o tests/TestReflectingBoundary tests/TestReflectingBoundary.cc $(LIBS) $(PYTHONLIBS) $(SPHERALLIBS) $(SPHERALLIBS) -L$(BLITZDIR) -lblitz -lm -lreadline -ldl

# 	#$(CXX) $(CXXFLAGS) -I/home/owen/Spheral++ -o tests/TestReflectingBoundary tests/TestReflectingBoundary.cc $(LIBS) `cat $(SRCTOP)/linklist | tr " " "\n" | grep -v "Methods" | grep -v "Pyffle" | tr "\n" " "`

# 	#$(CXX) $(CXXFLAGS) -I/home/owen/Spheral++ -o tests/TestReflectingBoundary tests/TestReflectingBoundary.cc ${SRCTARGETS:.cc=.o} $(LIBS) $(SRCTOP)/lib/Spheralmodule.so /usr/local/apps/kull/src/Python-1.5.2/libpython1.5.a -L/usr/local/apps/kull/lib_v1.0 -L/usr/local/lib -L/usr/lib -lreadline -lnsl -lsocket -ldl -lmalloc -ltk4.2 -ltcl7.6 -lX11
