# -*-makefile-*-
# Spheral++ Geometry package makefile
#-------------------------------------------------------------------------------
PKGNAME = Geometry
PKGDIR = $(PKGNAME)/
LIBTARGET = libSpheral_$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	Dimension.cc \
	GeomVector.cc \
	GeomTensor.cc \
	GeomSymmetricTensor.cc \
	Geom3Vector.cc \
	GeomThirdRankTensor.cc \
	GeomPlane.cc \
	Jacobi.cc \
	Jacobi2.cc \
	buildEigenVector.cc \
	EigenStruct.cc \
	computeEigenValues.cc \
	fakemain.cc \
	GeomFacet2d.cc \
	GeomFacet3d.cc \
	GeomPolygon.cc \
	GeomPolyhedron.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAINCLUDES = -I $(STLPORTDIR)/stlport -I $(BOOSTDIR) -I $(BOOSTDIR)/boost
##EXTRALIBS = -L$(STLPORTDIR) -lstlport_gcc $(SRCTOP)/boost/libboost_python.a
##EXTRAINCLUDES = -I $(BOOSTDIR) # -I $(BOOSTDIR)/boost
##EXTRALIBS = $(SRCTOP)/boost/libboost_python.a
##EXTRALIBS = -bmap:libGeometry.map /lib/crt0.o

EXTRALIBS = -lqhullstatic

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

#-------------------------------------------------------------------------------
TestGeometry:	tests/TestGeometry.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -I/home/owen/Spheral++/Geometry -o tests/TestGeometry tests/TestGeometry.cc
TestMTL:	TestMTL.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -I.. -o TestMTL TestMTL.cc
