# -*-makefile-*-
# Spheral++ FileIO package makefile
#-------------------------------------------------------------------------------
PKGNAME = FileIO
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	FileIO.cc \
	FlatFileIO.cc \
	PyFileIO.cc \
	vectorstringUtilities.cc
#	DbFileIO.cc 
#	HDF5IO.cc

PYTHONTARGETS = \
	GzipFileIO.py

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

# Include the python path.
EXTRAINCLUDES = -I$(PYTHONROOT)/include/python$(PYTHONVERSION)

# Postgres includes and libs -- should not be necessary.
#EXTRAINCLUDES = -I/usr/include/postgresql/libpq++ -I/usr/include/postgresql
#EXTRALIBS = -lpq++ -lpq

#EXTRAINCLUDES = -I${HDFINCDIR} 
#EXTRALIBS = -L${HDFLIBDIR} -lhdf5_cpp -lhdf5
#EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master

# For some reason, optimizing FlatFileIO on AIX with g++ 3.2.1 causes a segv 
# when we try and new up the fstream in open()?!  For now just override to build
# debug in this package, since this package is going to be disk access limited
# anyway.
# OPTFLAGS = -g
