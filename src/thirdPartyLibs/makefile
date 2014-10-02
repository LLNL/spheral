# -*-makefile-*-
# Spheral++ thirdPartyLibs package makefile
#-------------------------------------------------------------------------------
SPHERALTOP = ..

#-------------------------------------------------------------------------------
# Include the master file to get compilers and such.
include $(SPHERALTOP)/helpers/makefile_master

GZIP = gzip
BZIP2 = bzip2
UNZIP = unzip
TAR = tar
CP = cp
MV = mv
PATCH = patch
MKDIR = mkdir -p
CURL = curl -k -L
AR = ar

#-------------------------------------------------------------------------------
# Targets.
#-------------------------------------------------------------------------------
# Python
PYTHONDIST = Python-2.7.6.tgz
PYTHONURL = http://www.python.org/ftp/python/2.7.6/$(PYTHONDIST)
PYTHONDIR = $(PYTHONDIST:.tgz=)
PYTHONBUILDDATE = $(patsubst %.tgz, .%.date, $(PYTHONDIST))

# Numpy
NUMPYDIST = numpy-1.7.2.tar.gz
NUMPYURL = http://sourceforge.net/projects/numpy/files/NumPy/1.7.2/$(NUMPYDIST)/download
NUMPYDIR = $(NUMPYDIST:.tar.gz=)
NUMPYBUILDDATE = $(patsubst %.tar.gz, .%.date, $(NUMPYDIST))

# PyGnuplot
GNUPLOTDIST = gnuplot-py-1.8.tar.gz
GNUPLOTURL = http://downloads.sourceforge.net/gnuplot-py/gnuplot-py-1.8.tar.gz?modtime=1209767882&big_mirror=0
GNUPLOTDIR = $(GNUPLOTDIST:.tar.gz=)
GNUPLOTBUILDDATE = $(patsubst %.tar.gz, .%.date, $(GNUPLOTDIST))

# Psyco python accelerator.
PSYCOURL = http://wyvern.cs.uni-duesseldorf.de/psyco/psyco-snapshot.tar.gz
PSYCODIR = codespeak.net
PSYCOBUILDDIR = $(PYSCODIR)/svn/psyco/dist
PSYCOBUILDDATE = .psyco.date

# # Element tree
# ELEMENTTREEDIST = elementtree-1.2.7-20070827-preview.zip
# ELEMENTTREEURL = http://effbot.org/media/downloads/$(ELEMENTTREEDIST)
# ELEMENTTREEDIR = $(ELEMENTTREEDIST:.zip=)
# ELEMENTTREEBUILDDATE = $(patsubst %.zip, .%.date, $(ELEMENTTREEDIST))

# Boost
BOOSTDIST = boost_1_55_0.tar.bz2
BOOSTURL = http://sourceforge.net/projects/boost/files/boost/1.55.0/$(BOOSTDIST)/download
BOOSTDIR = $(BOOSTDIST:.tar.bz2=)
BOOSTDATE = $(patsubst %.tar.bz2, .%.date, $(BOOSTDIST))
BOOSTPYTHONLIB0 = libboost_python$(BOOSTEXT).$(DYLIBEXT)
BOOSTPYTHONLIB = libboost_python.$(DYLIBEXT)

# Cmake
CMAKEEXE = $(SPHERALDIR)/bin/cmake
CMAKEDIST = cmake-3.0.2.tar.gz
CMAKEURL = http://www.cmake.org/files/v3.0/$(CMAKEDIST)
CMAKEDIR = $(CMAKEDIST:.tar.gz=)

# Gccxml
GCCXMLDIR = $(GCCXMLDIST:.tar.bz2=)
GCCXMLDATE = $(patsubst %.tar.bz2, .%.date, $(GCCXMLDIST))

# Parmetis
PARMETISDIST = parmetis-4.0.3.tar.gz
PARMETISURL = http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/$(PARMETISDIST)
PARMETISDIR = $(PARMETISDIST:.tar.gz=)
PARMETISBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PARMETISDIST))
#PARMETISMAKEFILE = Makefile.ParMETIS

# Tau
TAUDIST = tau-$(TAUVERSION).tar.gz
TAUURL=http://www.cs.uoregon.edu/Research/tau/tau_releases/$(TAUDIST)
#TAUURL=http://tau.uoregon.edu/tau.tgz
TAUDIR = $(TAUDIST:.tar.gz=)
TAUBUILDDATE = $(patsubst %.tar.gz, .%.date, $(TAUDIST))
TAUDUMMYDATE = $(patsubst %.tar.gz, .%.dummydate, $(TAUDIST))

# Pygtk
PYGTKDIST = pygtk-2.12.1.tar.gz
PYGTKURL = http://saimei.acc.umu.se/pub/GNOME/sources/pygtk/2.12/$(PYGTKDIST)
PYGTKDIR = $(PYGTKDIST:.tar.gz=)
PYGTKBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PYGTKDIST))

# Matplotlib
MATPLOTLIBDIST = matplotlib-1.3.1.tar.gz
MATPLOTLIBURL = http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.3.1/matplotlib-1.3.1.tar.gz/download
MATPLOTLIBDIR = $(MATPLOTLIBDIST:.tar.gz=)
MATPLOTLIBBUILDDATE = $(patsubst %.tar.gz, .%.date, $(MATPLOTLIBDIST))

# Scipy
SCIPYDIST = scipy-0.13.2.tar.gz
SCIPYURL = http://sourceforge.net/projects/scipy/files/scipy/0.13.2/$(SCIPYDIST)/download
SCIPYDIR = $(SCIPYDIST:.tar.gz=)
SCIPYBUILDDATE = $(patsubst %.tar.gz, .%.date, $(SCIPYDIST))

# F2PY
# F2PYDIST = F2PY-2.45.241_1926.tar.gz
# F2PYDIR = $(F2PYDIST:.tar.gz=)
# F2PYBUILDDATE = $(patsubst %.tar.gz, .%.date, $(F2PYDIST))

# GNU scientific library (for special functions).
GSLDIST = gsl-1.14.tar.gz
GSLURL = http://ftp.gnu.org/gnu/gsl/$(GSLDIST)
GSLDIR = $(GSLDIST:.tar.gz=)
GSLBUILDDATE = $(patsubst %.tar.gz, .%.date, $(GSLDIST))
GSLLIBSO = libgsl.$(DYLIBEXT)
GSLLIBA = libgsl.a

# pygsl -- python interface to GNU scientific library
PYGSLDIST = pygsl-0.9.5.tar.gz
PYGSLURL = http://sourceforge.net/projects/pygsl/files/pygsl/pygsl-0.9.5/$(PYGSLDIST)/download
PYGSLDIR = $(PYGSLDIST:.tar.gz=)
PYGSLBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PYGSLDIST))

# PETSc linear solver library.
PETSCURL = ftp://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-2.3.3-p8.tar.gz
PETSCDIST = petsc-2.3.3-p8.tar.gz
PETSCBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PETSCDIST))

# petsc4py
PETSC4PYDIST = petsc4py-0.7.5.tar.gz
PETSC4PYURL = http://petsc4py.googlecode.com/files/$(PETSC4PYDIST)
PETSC4PYDIR = $(PETSC4PYDIST:.tar.gz=)
PETSC4PYBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PETSC4PYDIST))

# Google hash_map.
HASHMAPDIST = sparsehash-1.1.tar.gz
HASHMAPURL = http://google-sparsehash.googlecode.com/files/$(HASHMAPDIST)
HASHMAPDIR = $(HASHMAPDIST:.tar.gz=)
HASHMAPBUILDDATE = $(patsubst %.tar.gz, .%.date, $(HASHMAPDIST))

# Visit writer software from Hank.
VISITWRITERDIST = visit_writer.tar.gz
VISITWRITERURL = http://sourceforge.net/projects/spheral/files/ThirdPartySources/$(VISITWRITERDIST)/download
VISITWRITERDIR = $(VISITWRITERDIST:.tar.gz=)
VISITWRITERBUILDDATE = $(patsubst %.tar.gz, .%.date, $(VISITWRITERDIST))

# ATS testing system.
ATSDIST = ats-5.2.tar.gz
ATSURL = http://ats.googlecode.com/files/$(ATSDIST)
ATSDIR = $(ATSDIST:.tar.gz=)
ATSBUILDDATE = $(patsubst %.tar.gz, .%.date, $(ATSDIST))

# FFTW
FFTWDIST = fftw-3.3.tar.gz
FFTWURL = http://www.fftw.org/$(FFTWDIST)
FFTWDIR = $(FFTWDIST:.tar.gz=)
FFTWBUILDDATE = $(patsubst %.tar.gz, .%.date, $(FFTWDIST))

# FFTW++
FFTWPPDIST = fftw++-1.09.tar.gz
FFTWPPURL = http://sourceforge.net/projects/fftwpp/files/$(FFTWPPDIST)/download
FFTWPPDIR = $(FFTWPPDIST:.tar.gz=)
FFTWPPBUILDDATE = $(patsubst %.tar.gz, .%.date, $(FFTWPPDIST))

# Array
ARRAYDIST = Array.h
ARRAYURL = http://www.math.ualberta.ca/~bowman/$(ARRAYDIST)
ARRAYBUILDDATE = $(patsubst %.h, .%.date, $(ARRAYDIST))

# ZLIB
ZLIBDIST = zlib-1.2.8.tar.gz
ZLIBURL = http://zlib.net/$(ZLIBDIST)
ZLIBDIR = $(ZLIBDIST:.tar.gz=)
ZLIBBUILDDATE = $(patsubst %.tar.gz, .%.date, $(ZLIBDIST))

# JPEG
JPEGDIST = jpegsrc.v8d.tar.gz
JPEGURL = http://www.ijg.org/files/$(JPEGDIST)
JPEGDIR = jpeg-8d
JPEGBUILDDATE = $(patsubst %.tar.gz, .%.date, $(JPEGDIST))

# PIL
PILDIST = Imaging-1.1.6.tar.gz
PILURL = http://effbot.org/downloads/$(PILDIST)
PILDIR = $(PILDIST:.tar.gz=)
PILBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PILDIST))

# pybindgen
PYBINDGENDIST = PyBindGen-0.17.0.tar.gz
PYBINDGENURL = https://pypi.python.org/packages/source/P/PyBindGen/$(PYBINDGENDIST)
PYBINDGENDIR = $(PYBINDGENDIST:.tar.gz=)
PYBINDGENBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PYBINDGENDIST))

# WildMagic
WILDMAGICDIST = WildMagic5p9.zip
WILDMAGICURL = http://sourceforge.net/projects/spheral/files/ThirdPartySources/$(WILDMAGICDIST)/download
WILDMAGICDIR = GeometricTools
WILDMAGICBUILDDATE = $(patsubst %.zip, .%.date, $(WILDMAGICDIST))

# Qhull
QHULLDIST = qhull-2012.1-src.tgz
QHULLURL = http://www.qhull.org/download/$(QHULLDIST)
QHULLDIR = $(QHULLDIST:-src.tgz=)
QHULLINCDIR = $(SPHERALDIR)/include/qhull
QHULLDATE = $(patsubst %-src.tgz, .%.date, $(QHULLDIST))

# mpi4py
MPI4PYDIST = mpi4py-1.3.tar.gz
MPI4PYURL = http://mpi4py.googlecode.com/files/$(MPI4PYDIST)
MPI4PYDIR = $(MPI4PYDIST:.tar.gz=)
MPI4PYBUILDDATE = $(patsubst %.tar.gz, .%.date, $(MPI4PYDIST))

# HDF5
HDF5DIST = hdf5-1.8.13.tar.bz2
HDF5URL = http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/$(HDF5DIST)
HDF5DIR = $(HDF5DIST:.tar.bz2=)
HDF5BUILDDATE = $(patsubst %.tar.bz2, .%.date, $(HDF5DIST))

# SILO
SILODIST = silo-4.9.1-bsd.tar.gz
SILOURL = https://wci.llnl.gov/codes/silo/silo-4.9.1/$(SILODIST)
SILODIR = $(SILODIST:.tar.gz=)
SILOBUILDDATE = $(patsubst %.tar.gz, .%.date, $(SILODIST))

# Eigen
EIGENDIST = Eigen-3.2.0.tar.bz2
EIGENURL = http://bitbucket.org/eigen/eigen/get/3.2.0.tar.bz2
EIGENDIR = eigen-eigen-ffa86ffb5570
EIGENBUILDDATE = .eigen.date

# # Voro++
# VOROPPDIST = voro++_0.3.1.tar.gz
# VOROPPURL = http://math.lbl.gov/voro++/download/dir/$(VOROPPDIST)
# VOROPPDIR = voro++
# VOROPPBUILDDATE = $(patsubst %.tar.gz, .%.date, $(VOROPPDIST))

# # Voro++ 2d specialization
# VOROPP2DDIST = voro++_2d.tbz
# VOROPP2DDIR = voro++_2d
# VOROPP2DINC = voro++_2d
# VOROPP2DBUILDDATE = $(patsubst %.tbz, .%.date, $(VOROPP2DDIST))

# GNU multiprecision library.
GMPDIST = gmp-5.0.2.tar.bz2
GMPURL = ftp://ftp.gmplib.org/pub/gmp-5.0.2/$(GMPDIST)
GMPDIR = $(GMPDIST:.tar.bz2=)
GMPBUILDDATE = $(patsubst %.tar.bz2, .%.date, $(GMPDIST))

# mpmath -- a python multiprecision math library
MPMATHDIST = mpmath-0.17.tar.gz
MPMATHURL = http://mpmath.googlecode.com/files/$(MPMATHDIST)
MPMATHDIR = $(MPMATHDIST:.tar.gz=)
MPMATHBUILDDATE = $(patsubst %.tar.gz, .%.date, $(MPMATHDIST))

# # OCCTREE 
# OCTTREEDIST = octree_code.zip
# OCTTREEURL = http://nomis80.org/code/$(OCTTREEDIST)
# OCTTREEDIR = $(OCTTREEDIST:.zip=)
# OCTTREEINCDIR = $(SPHERALDIR)/include/octtree
# OCTTREEBUILDDATE = $(patsubst %.zip, .%.date, $(OCTTREEDIST))

# The FLANN (Fast Library for Approximate Nearest Neighbors) library.
FLANNDIST = flann-1.7.1-src.zip
FLANNURL = http://people.cs.ubc.ca/~mariusm/uploads/FLANN/$(FLANNDIST)
FLANNDIR = $(FLANNDIST:.zip=)
FLANNBUILDDATE = $(patsubst %.zip, .%.date, $(FLANNDIST))

# The Point Cloud Library (PCL)
PCLDIST = PCL-1.5.1-Source.tar.bz2
PCLURL = http://dev.pointclouds.org/attachments/download/771/$(PCLDIST)
PCLDIR = $(PCLDIST:.tar.bz2=)
PCLBUILDDATE = $(patsubst %.tar.bz2, .%.date, $(PCLDIST))

# HYPRE parallel linear solver
HYPREDIST = hypre-2.8.0b.tar.gz
HYPREURL = https://computation.llnl.gov/casc/hypre/download/$(HYPREDIST)
HYPREDIR = $(HYPREDIST:.tar.gz=)
HYPREBUILDDATE = $(patsubst %.tar.gz, .%.date, $(HYPREDIST))

# Triangle (needed by polytope)
TRIANGLEDIST = triangle.zip
TRIANGLEURL = http://www.netlib.org/voronoi/$(TRIANGLEDIST)
TRIANGLEDIR = $(TRIANGLEDIST:.zip=)

# Tetgen (needed by polytope)
TETGENDIST = tetgen1.4.3.tar.gz
TETGENURL = http://ftp.mcs.anl.gov/pub/petsc/externalpackages/$(TETGENDIST)
TETGENDIR = $(TETGENDIST:.tar.gz=)

# Polytope
POLYTOPEVER = 0.5.21
POLYTOPEDIST = polytope-$(POLYTOPEVER).tbz
POLYTOPEURL = https://bitbucket.org/jjphatt/polytope/downloads/$(POLYTOPEDIST)
POLYTOPEDIR = $(POLYTOPEDIST:.tbz=)
POLYTOPEBUILDDATE = $(patsubst %.tbz, .%.date, $(POLYTOPEDIST))

# OpenSubDiv
OPENSUBDIVDIST = OpenSubdiv-2_5_0.zip
OPENSUBDIVURL = https://github.com/PixarAnimationStudios/OpenSubdiv/archive/v2_5_0.zip
OPENSUBDIVDIR = $(OPENSUBDIVDIST:.zip=)
OPENSUBDIVBUILDDATE = $(patsubst %.zip, .%.date, $(OPENSUBDIVDIST))

#-------------------------------------------------------------------------------
all:	$(PYTHONBUILDDATE) \
	$(BOOSTDATE) \
	$(GCCXMLTARGETS) \
	$(MATPLOTLIBTARGETS) \
	$(SCIPYTARGETS) \
	$(PETSCTARGETS) \
	$(VISITWRITERBUILDDATE) \
	$(ATSBUILDDATE) \
	$(FFTWBUILDDATE) \
	$(FFTWPPBUILDDATE) \
	$(ARRAYBUILDDATE) \
	$(ZLIBBUILDDATE) \
	$(PYBINDGENBUILDDATE) \
	$(QHULLDATE) \
	$(HDF5BUILDDATE) \
	$(SILOBUILDDATE) \
	$(EIGENBUILDDATE) \
	$(MPMATHBUILDDATE) \
	$(POLYTOPEBUILDDATE) \
	$(OPENSUBDIVTARGETS) \
	$(EXTRATHIRDPARTYTARGETS) \
#	$(HYPREBUILDDATE)
#	$(FLANNBUILDDATE) \
#	$(PCLBUILDDATE)

clean:
	rm -fr $(PYTHONDIR) $(PYTHONBUILDDATE) \
	       $(NUMPYDIR) $(NUMPYBUILDDATE) \
	       $(GNUPLOTDIR) $(GNUPLOTBUILDDATE) \
	       $(PSYCODIR) $(PSYCOBUILDDATE) \
	       $(BOOSTDIR) $(BOOSTDATE) \
               $(CMAKEDIR) \
               $(GCCXMLDIR) $(GCCXMLDATE) \
	       $(PSYCODIR) \
	       $(MPI4PYDIR) $(MPI4PYBUILDDATE) \
               $(PARMETISDIR) $(PARMETISBUILDDATE) \
	       $(SCIPYDIR) $(SCIPYBUILDDATE) \
	       $(GSLDIR) \
	       $(PYGSLDIR) \
               $(MATPLOTLIBDIR) $(MATPLOTLIBBUILDDATE) \
	       $(HASHMAPDATE) \
	       $(VISITWRITERDIR) $(VISITWRITERBUILDDATE) \
	       $(ATSWRITERDIR) $(ATSWRITERBUILDDATE) \
	       $(FFTWDIR) $(FFTWBUILDDATE) \
	       $(FFTWPPDIR) $(FFTWPPBUILDDATE) \
	       $(ARRAYDIR) $(ARRAYBUILDDATE) \
	       $(ZLIBDIR) $(ZLIBBUILDDATE) \
	       $(PYBINDGENDIR) $(PYBINDGENBUILDDATE) \
	       $(QHULLDIR) $(QHULLDATE) \
	       $(HDF5DIR) $(HDF5BUILDDATE) \
	       $(SILODIR) $(SILOBUILDDATE) \
	       $(EIGENDIR) $(EIGENBUILDDATE) \
	       $(MPMATHDIR) $(MPMATHBUILDDATE) \
	       $(POLYTOPEDIR) $(POLYTOPEBUILDDATE) \
	       $(HYPREDIR) $(HYPREBUILDDATE) \
	       $(FLANNDIR) $(FLANNBUILDDATE) \
	       $(PCLDIR) $(PCLBUILDDATE) \
	       $(OPENSUBDIVDIR) $(OPENSUBDIVBUILDDATE) \
	       $(EXTRATHIRDPARTYTARGETS) \
	       .*.date

.SUFFIXES:	.date .gz .tgz

#-------------------------------------------------------------------------------
# Python
#-------------------------------------------------------------------------------
$(PYTHONDIST):
	$(CURL) "$(PYTHONURL)" > $(PYTHONDIST)

$(PYTHONBUILDDATE):	$(PYTHON) $(PYTHONDIST)
	@echo "################################################################################"
	@echo "# Building $(PYTHONDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(PYTHONDIR)
	$(GZIP) -dc $(PYTHONDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(PYTHONDIR); \
	env PATH='$(SPHERALDIR)/bin:$$PATH' CONFIG_SHELL=$(CONFIG_SHELL) CC='$(PYTHONCC)' CXX='$(PYTHONCXX)' CFLAGS='$(PYTHONCFLAGS)' LDFLAGS='$(LDFLAGS)' HAS_HG=no ./configure --disable-ipv6 --with-gcc='$(PYTHONCC)' --with-cxx-main=$(PYTHONCXX) $(PYTHONCONFFLAGS) --prefix=$(SPHERALDIR)"
	$(MAKE) -C $(PYTHONDIR) CC='$(PYTHONCC)' CXX='$(PYTHONCXX)' CFLAGS='$(PYTHONCFLAGS)' LDFLAGS='$(LDFLAGS)'
	$(MAKE) -C $(PYTHONDIR) install
	@if ! test -e $(PYTHONEXE); then \
	echo "ERROR building local python distribution."; \
	rm -f $(PYTHONBUILDDATE); \
	exit 1; \
	fi
	$(SRCTOP)/helpers/runexecstack $(PYTHONEXE)
	rm -fr $(PYTHONDIR)
	touch $(PYTHONBUILDDATE)

#-------------------------------------------------------------------------------
# Numpy python extension
#-------------------------------------------------------------------------------
$(NUMPYDIST):
	$(CURL) "$(NUMPYURL)" > $(NUMPYDIST)

$(NUMPYBUILDDATE):	$(PYTHONBUILDDATE) $(NUMPYDIST)
	@echo "################################################################################"
	@echo "# Building $(NUMPYDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	$(MAKE) $(NUMPYDIST)
	rm -fr $(NUMPYDIR)
	$(GZIP) -dc $(NUMPYDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(NUMPYDIR); $(PYTHONEXE) setup.py install;'
	@if test -n "`$(PYTHONEXE) -c 'import numpy'`"; then \
	echo "ERROR building numpy python extension."; \
	rm -f $(NUMPYBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(NUMPYDIR)
	touch $(NUMPYBUILDDATE)

#-------------------------------------------------------------------------------
# Gnuplot python extension
#-------------------------------------------------------------------------------
$(GNUPLOTDIST):
	$(CURL) "$(GNUPLOTURL)" > $(GNUPLOTDIST)

$(GNUPLOTBUILDDATE):	$(PYTHONBUILDDATE) $(NUMPYBUILDDATE) $(GNUPLOTDIST)
	@echo "################################################################################"
	@echo "# Building $(GNUPLOTDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(GNUPLOTDIR)
	$(GZIP) -dc $(GNUPLOTDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(GNUPLOTDIR); $(PYTHONEXE) setup.py install;'
	@if test -n "`$(PYTHONEXE) -c 'import Gnuplot'`"; then \
	echo "ERROR building Gnuplot python extension."; \
	rm -f $(GNUPLOTBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(GNUPLOTDIR)
	touch $(GNUPLOTBUILDDATE)

#-------------------------------------------------------------------------------
# Psyco python accelerator
# Note that psyco only works on intel processors, so we try and detect that
# here.
#-------------------------------------------------------------------------------
$(PSYCODIST):
	$(CURL) "$(PSYCOURL)" > $(PYSCODIST)

$(PSYCOBUILDDATE):	$(PYTHONBUILDDATE) $(PSYCODIST)
	@echo "################################################################################"
	@echo "# Building $(PSYCODIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	$(SHELL) -ec 'cd $(PSYCOBUILDDIR); $(PYTHONEXE) setup.py install;'
	@if test -n "`$(PYTHONEXE) -c 'import psyco'`"; then \
	echo "ERROR building Psyco python extension."; \
	rm -f $(PSYCOBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PSYCODIR)
	touch $(PSYCOBUILDDATE)

# #-------------------------------------------------------------------------------
# # Elementtree python extension
# #-------------------------------------------------------------------------------
# $(ELEMENTTREEDIST):
# 	$(CURL) "$(ELEMENTTREEURL)" > $(ELEMENTTREEDIST)

# $(ELEMENTTREEBUILDDATE):	$(PYTHONBUILDDATE) $(ELEMENTTREEDIST)
# 	@echo "################################################################################"
# 	@echo "# Building $(ELEMENTTREEDIST)"
# 	@echo "################################################################################"
# 	rm -fr $(ELEMENTTREEDIR)
# 	$(UNZIP) $(ELEMENTTREEDIST)
# 	$(SHELL) -ec 'cd $(ELEMENTTREEDIR); $(PYTHONEXE) setup.py install;'
# 	@if test -n "`$(PYTHONEXE) -c 'import elementtree'`"; then \
# 	echo "ERROR building elementtree python extension."; \
# 	rm -f $(ELEMENTTREEBUILDDATE); \
# 	exit 1; \
# 	fi
# 	rm -fr $(ELEMENTTREEDIR)
# 	touch $(ELEMENTTREEBUILDDATE)

#-------------------------------------------------------------------------------
# Boost
#-------------------------------------------------------------------------------
$(BOOSTDIST):
	$(CURL) "$(BOOSTURL)" > $(BOOSTDIST)

$(BOOSTDATE):	$(BOOSTDIST) $(ZLIBBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(BOOSTDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(BOOSTDIR)
	$(BZIP2) -dc $(BOOSTDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(BOOSTDIR); \
	PATH=$(SPHERALDIR)/bin:$$PATH; export PATH; \
	./bootstrap.sh --without-libraries=atomic,coroutine,log,chrono,context,date_time,exception,filesystem,graph,graph_parallel,iostreams,locale,math,mpi,program_options,python,random,regex,serialization,signals,system,test,thread,timer,wave --prefix=$(SPHERALDIR); \
	./b2 install'
	@if ! test -e "$(SPHERALDIR)/include/boost" ; then \
	echo "ERROR building local boost distribution. (0, no include dir)"; \
	rm -f $(BOOSTDATE); \
	exit 1; \
	fi
	@if test "$(BOOSTPYTHONTARGET)" = "python" -a ! -e $(SPHERALDIR)/lib/$(BOOSTPYTHONLIB0) ; then \
	echo "ERROR building local boost distribution. (0, no $(BOOSTPYTHONLIB0))"; \
	rm -f $(BOOSTDATE); \
	exit 1; \
	fi
	@if test "$(BOOSTPYTHONTARGET)" = "python"; then \
	rm -f $(LIBDIR)/$(BOOSTPYTHONLIB) $(SPHERALDIR)/include/boost; \
	ln -s $(SPHERALDIR)/lib/$(BOOSTPYTHONLIB0) $(LIBDIR)/$(BOOSTPYTHONLIB); \
	fi
	@if test "$(BOOSTPYTHONTARGET)" = "python" -a ! -e $(LIBDIR)/$(BOOSTPYTHONLIB); then \
	echo "ERROR building local boost distribution. (2)"; \
	rm -f $(BOOSTDATE); \
	exit 1; \
	fi
	rm -fr $(BOOSTDIR)
	touch $(BOOSTDATE)

#-------------------------------------------------------------------------------
# gccxml
#-------------------------------------------------------------------------------
$(GCCXMLDATE):		$(CMAKEEXE) $(GCCXMLDIST)
	@echo "################################################################################"
	@echo "# Building $(GCCXMLDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(GCCXMLDIR)
	$(BZIP2) -dc $(GCCXMLDIST) | $(TAR) -xf -
	$(MKDIR) $(GCCXMLDIR)/build
	$(SHELL) -ec 'cd $(GCCXMLDIR)/build; \
	env CC='$(GCCXMLCC)' CXX='$(GCCXMLCXX)' $(CMAKEEXE) -DCMAKE_INSTALL_PREFIX:PATH=$(SPHERALDIR) ..; $(MAKE); $(MAKE) install'
	@if ! test -e $(SPHERALDIR)/bin/gccxml_cc1plus -a -e $(SPHERALDIR)/bin/gccxml; then \
	echo "ERROR building gccxml."; \
	rm -f $(GCCXMLDATE); \
	exit 1; \
	fi
	rm -fr $(GCCXMLDIR)
	touch $(GCCXMLDATE)

#-------------------------------------------------------------------------------
# Cmake
#-------------------------------------------------------------------------------
$(CMAKEDIST):
	$(CURL) "$(CMAKEURL)" > $(CMAKEDIST)

$(CMAKEEXE):	$(CMAKEDIST)
	@echo "################################################################################"
	@echo "# Building $(CMAKEDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(CMAKEDIR)
	$(GZIP) -dc $(CMAKEDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(CMAKEDIR); \
	env CC='$(CMAKECC)' CXX='$(CMAKECXX)' ./bootstrap --prefix=$(SPHERALDIR); \
	$(MAKE); \
	$(MAKE) install'
	@if ! test -e $(CMAKEEXE); then \
	echo "ERROR building cmake."; \
	exit 1; \
	fi
	rm -fr $(CMAKEDIR)

#-------------------------------------------------------------------------------
# ParMETIS
#-------------------------------------------------------------------------------
$(PARMETISDIST):
	$(CURL) "$(PARMETISURL)" > $(PARMETISDIST)

$(PARMETISBUILDDATE):  $(PARMETISDIST) $(CMAKEEXE)
	@echo "################################################################################"
	@echo "# Building $(PARMETISDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(PARMETISDIR)
	$(GZIP) -dc $(PARMETISDIST) | $(TAR) -xf -
	env PATH="$(SPHERALDIR)/bin:$$PATH" CC="$(MPICC)" CXX="$(MPICXX)" CFLAGS='$(MPICCFLAGS)' CXXFLAGS='$(MPICXXFLAGS)' $(MAKE) -C $(PARMETISDIR) config shared=1 prefix=$(SPHERALDIR) cc=$(MPICC) cxx=$(MPICXX)
	$(MAKE) MAKEFLAGS= -C $(PARMETISDIR) 
	$(MAKE) MAKEFLAGS= -C $(PARMETISDIR) install
	$(CP) -f $(PARMETISDIR)/metis/include/metis.h $(SPHERALDIR)/include
	@if ! test -e $(SPHERALDIR)/lib/libparmetis.$(DYLIBEXT) -a \
	           -e $(SPHERALDIR)/include/parmetis.h -a \
	           -e $(SPHERALDIR)/include/metis.h; then \
	echo "ERROR building ParMETIS."; \
	rm -f $(PARMETISBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PARMETISDIR)
	touch $(PARMETISBUILDDATE)

#-------------------------------------------------------------------------------
# Tau
#-------------------------------------------------------------------------------
$(TAUDIST):
	$(CURL) "$(TAUURL)" > $(TAUDIST)

$(TAUBUILDDATE):  $(TAUDIST)
	@echo "################################################################################"
	@echo "# Building $(TAUDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(TAUDIR) $(TAUDUMMYDATE)
	$(GZIP) -dc $(TAUDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(TAUDIR); ./configure $(TAUCONFIGUREFLAGS)"
	$(SHELL) -ec "$(MAKE) $(AIXSHELL) -C $(TAUDIR) install"
	$(SHELL) -ec "cd $(TAUDIR); ln -s $(TAUARCH)/lib/Makefile.tau* Makefile.tau"
	touch $(TAUDIR)/include/TAU.h
	touch $(TAUBUILDDATE)

#-------------------------------------------------------------------------------
# Tau -- dummy (non-compiled) version.
#-------------------------------------------------------------------------------
$(TAUDUMMYDATE):  $(TAUDIST)
	@echo "################################################################################"
	@echo "# Building $(TAUDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(TAUDIR) $(TAUBUILDDATE)
	$(GZIP) -dc $(TAUDIST) | $(TAR) -xf -
	touch $(TAUDIR)/include/TAU.h
	touch $(TAUDUMMYDATE)

#-------------------------------------------------------------------------------
# PyGTK
#-------------------------------------------------------------------------------
$(PYGTKDIST):
	$(CURL) "$(PYGTKURL)" > $(PYGTKDIST)

$(PYGTKBUILDDATE):  $(NUMPYBUILDDATE) $(PYGTKDIST)
	@echo "################################################################################"
	@echo "# Building $(PYGTKDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(PYGTKDIR) $(PYGTKDUMMYDATE)
	$(GZIP) -dc $(PYGTKDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(PYGTKDIR); \
	env CC='$(PYTHONCC)' CXX='$(PYTHONCXX)' PYTHON=$(PYTHONEXE) ./configure --enable-thread --prefix=$(SPHERALDIR); \
	$(MAKE); \
	$(MAKE) install"
	@if test -n "`$(PYTHONEXE) -c 'import pygtk'`"; then \
	echo "ERROR building PyGTK extension."; \
	rm -f $(PYGTKBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PYGTKDIR)
	touch $(PYGTKBUILDDATE)

#-------------------------------------------------------------------------------
# matplotlib
#-------------------------------------------------------------------------------
$(MATPLOTLIBDIST):
	$(CURL) "$(MATPLOTLIBURL)" > $(MATPLOTLIBDIST)

$(MATPLOTLIBBUILDDATE): $(MATPLOTLIBDIST)
	@echo "################################################################################"
	@echo "# Building $(MATPLOTLIBDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(MATPLOTLIBDIR) $(MATPLOTLIBDUMMYDATE)
	$(GZIP) -dc $(MATPLOTLIBDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(MATPLOTLIBDIR); \
	env PKG_CONFIG_PATH=$(SPHERALDIR)/lib/pkgconfig $(PYTHONEXE) setup.py build; \
	$(PYTHONEXE) setup.py install"
	@if test -n "`$(PYTHONEXE) -c 'import matplotlib'`"; then \
	echo "ERROR building matplotlib extension."; \
	rm -f $(MATPLOTLIBBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(MATPLOTLIBDIR)
	touch $(MATPLOTLIBBUILDDATE)

#-------------------------------------------------------------------------------
# SciPy
#-------------------------------------------------------------------------------
$(SCIPYDIST):
	$(CURL) "$(SCIPYURL)" > $(SCIPYDIST)

$(SCIPYBUILDDATE):  $(PYTHONBUILDDATE) $(NUMPYBUILDDATE) $(SCIPYDIST)
	@echo "################################################################################"
	@echo "# Building $(SCIPYDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(SCIPYDIR) $(SCIPYDUMMYDATE)
	$(GZIP) -dc $(SCIPYDIST) | $(TAR) -xf -
#	$(PATCH) $(SCIPYDIR)/setup.py $(SCIPYDIR)-setup.py-patch 
	$(SHELL) -ec "cd $(SCIPYDIR); \
	env PKG_CONFIG_PATH=$(SPHERALDIR)/lib/pkgconfig $(PYTHONEXE) setup.py build; \
	$(PYTHONEXE) setup.py install"
	@if test -n "`$(PYTHONEXE) -c 'import scipy'`"; then \
	echo "ERROR building scipy extension."; \
	rm -f $(SCIPYBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(SCIPYDIR)
	touch $(SCIPYBUILDDATE)

#-------------------------------------------------------------------------------
# F2Py -- Needed by SciPy.
#-------------------------------------------------------------------------------
# $(F2PYBUILDDATE):  $(F2PYDIST) $(PYTHONBUILDDATE) $(NUMPYBUILDDATE) 
# 	@echo "################################################################################"
# 	@echo "# Building $(F2PYDIST)"
# 	@echo "################################################################################"
# 	rm -fr $(F2PYDIR) $(F2PYDUMMYDATE)
# 	$(GZIP) -dc $(F2PYDIST) | $(TAR) -xf -
# 	$(SHELL) -ec "cd $(F2PYDIR); \
# 	env PKG_CONFIG_PATH=$(SPHERALDIR)/lib/pkgconfig $(PYTHONEXE) setup.py build; \
# 	$(PYTHONEXE) setup.py install"
# 	@if test -n "`$(PYTHONEXE) -c 'import f2py2e'`"; then \
# 	echo "ERROR building f2py extension."; \
# 	rm -f $(F2PYBUILDDATE); \
# 	exit 1; \
# 	fi
# 	rm -fr $(F2PYDIR)
# 	touch $(F2PYBUILDDATE)

#-------------------------------------------------------------------------------
# GNU Scientific Library -- Needed for special functions.
#-------------------------------------------------------------------------------
$(GSLDIST):
	$(CURL) "$(GSLURL)" > $(GSLDIST)

$(GSLBUILDDATE): $(GSLDIST)
	@echo "################################################################################"
	@echo "# Building $(GSLDIST) $(GSLBUILDDATE)" 
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(GSLDIR)
	$(GZIP) -dc $(GSLDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(GSLDIR); \
	env CC='$(CC)' CXX='$(CXX)' ./configure --prefix=$(SPHERALDIR) ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/$(GSLLIBSO) -o \
	           -e $(SPHERALDIR)/lib/$(GSLLIBA) ; then \
	echo "ERROR building GNU scientific library."; \
	rm -f $(GSLBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(GSLDIR)
	touch $(GSLBUILDDATE)

#-------------------------------------------------------------------------------
# Python interface to GNU Scientific Library.
#-------------------------------------------------------------------------------
$(PYGSLDIST):
	$(CURL) "$(PYGSLURL)" > $(PYGSLDIST)

$(PYGSLBUILDDATE): $(PYGSLDIST) $(GSLBUILDDATE) $(PYTHONBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(PYGSLDIST) $(PYGSLBUILDDATE)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(PYGSLDIR)
	$(GZIP) -dc $(PYGSLDIST) | $(TAR) -xf -
	env PATH=$(SPHERALDIR)/bin:$$PATH $(SHELL) -ec 'cd $(PYGSLDIR); $(PYTHONEXE) setup.py install;'
	@if test -n "`$(PYTHONEXE) -c 'import pygsl'`"; then \
	echo "ERROR building pygsl extension."; \
	rm -f $(PYGSLBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PYGSLDIR)
	touch $(PYGSLBUILDDATE)

#-------------------------------------------------------------------------------
# PETSc linear solver library -- Needed for SPH Gravity and MHD packages.
#-------------------------------------------------------------------------------
$(PETSCIDIST):
	$(CURL) "$(PETSCURL)" > $(PETSCIDIST)

$(PETSCBUILDDATE):	$(PETSCIDIST)
	@echo "################################################################################"
	@echo "# Fetching and building $(PETSCDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	@if ! test -e $(PETSC_DIR) ; then \
	$(MKDIR) $(PETSC_DIR); \
	fi
	@if ! test -e $(PETSC_DIR)/lib/$(PETSC_ARCH)/libpetsc.a ; then \
	$(SHELL) -ec "cd $(PETSC_DIR)/.. ; \
	$(GZIP) -dc $(SRCTOP)/thirdPartyLibs/$(PETSCDIST) | $(TAR) -xf - ; \
   cd $(PETSC_DIR); \
	$(PYTHONEXE) config/configure.py --prefix=$(SPHERALDIR) $(PETSCOPTS) ; \
	$(MAKE) ; \
	$(MAKE) test ; \
	$(MAKE) python ; \
	$(MAKE) install"; \
	fi
	touch $(PETSCBUILDDATE)

#-------------------------------------------------------------------------------
# petsc4py Python module -- used by SPHGravity and MHD.
#-------------------------------------------------------------------------------
$(PETSC4PYDIST):
	$(CURL) "$(PETSC4PYURL)" > $(PETSC4PYDIST)

$(PETSC4PYBUILDDATE): $(PETSC4PYDIST)
	@echo "################################################################################"
	@echo "# Fetching and building $(PETSC4PYDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	$(GZIP) -dc $(PETSC4PYDIST) | $(TAR) -xf -
#	$(PATCH) $(SCIPYDIR)/setup.py $(SCIPYDIR)-setup.py-patch 
	$(SHELL) -ec "cd $(PETSC4PYDIR); \
	env CC=$(SRCTOP)/helpers/mpicc $(PYTHONEXE) setup.py install"
	@if test -n "`$(PYTHONEXE) -c 'import petsc4py'`"; then \
	echo "ERROR building petsc4py extension."; \
	rm -f $(PETSC4PYBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PETSC4PYDIR)
	touch $(PETSC4PYBUILDDATE)

#-------------------------------------------------------------------------------
# Google's hash map implementation.
#-------------------------------------------------------------------------------
$(HASHMAPDIST):
	$(CURL) "$(HASHMAPURL)" > $(HASHMAPDIST)

$(HASHMAPBUILDDATE):	$(HASHMAPDIST)
	@echo "################################################################################"
	@echo "# Building $(HASHMAPDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	$(GZIP) -dc $(HASHMAPDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(HASHMAPDIR); \
	env CC='$(CC)' CXX='$(CXX)' ./configure --enable-namespace=google --prefix=$(SPHERALDIR); \
	$(MAKE); \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/include/google; then \
	echo "ERROR building google hash_map."; \
	exit 1; \
	fi
	rm -fr $(HASHMAPBUILDDATE) $(HASHMAPDIR)
	touch $(HASHMAPBUILDDATE)

#-------------------------------------------------------------------------------
# Hank's visit_writer 
#-------------------------------------------------------------------------------
$(VISITWRITERDIST):
	$(CURL) "$(VISITWRITERURL)" > $(VISITWRITERDIST)

$(VISITWRITERBUILDDATE):  $(VISITWRITERDIST) $(POLYTOPEBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(VISITWRITERDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	$(GZIP) -dc $(VISITWRITERDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(VISITWRITERDIR); $(MAKE)"
	@if ! test -e $(SPHERALDIR)/include/visit_writer.h -a \
	           -e $(SPHERALDIR)/lib/libvisit_writer.$(DYLIBEXT); then \
	echo "ERROR building visit_writer."; \
	exit 1; \
	fi
	rm -fr $(VISITWRITERBUILDDATE) $(VISITWRITERDIR)
	touch $(VISITWRITERBUILDDATE)

#-------------------------------------------------------------------------------
# ATS testing system.
#-------------------------------------------------------------------------------
$(ATSDIST):
	$(CURL) "$(ATSURL)" > $(ATSDIST)

$(ATSBUILDDATE):	$(ATSDIST) $(PYTHONBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(ATSDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(ATSDIR)
	$(GZIP) -dc $(ATSDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(ATSDIR); \
	$(PYTHONEXE) setup.py install; \
	cd LC; \
	$(PYTHONEXE) setup.py install'
	@if ! test -e $(SPHERALDIR)/bin/ats; then \
	echo "ERROR building ATS."; \
	rm -f $(ATSBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(ATSDIR)
	touch $(ATSBUILDDATE)

#-------------------------------------------------------------------------------
# FFTW fast-fourier transform
#-------------------------------------------------------------------------------
$(FFTWDIST):
	$(CURL) "$(FFTWURL)" > $(FFTWDIST)

$(FFTWBUILDDATE):	$(FFTWDIST)
	@echo "################################################################################"
	@echo "# Building $(FFTWDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(FFTWDIR)
	$(GZIP) -dc $(FFTWDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(FFTWDIR); \
	env CC='$(CC)' CXX='$(CXX)' CFLAGS='$(CFLAGS)' ./configure --prefix=$(SPHERALDIR) ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/include/fftw3.h -a \
	           -e $(SPHERALDIR)/lib/libfftw3.a; then \
	echo "ERROR building FFTW."; \
	exit 1; \
	fi
	rm -fr $(FFTWBUILDDATE) $(FFTWDIR)
	touch $(FFTWBUILDDATE)

#-------------------------------------------------------------------------------
# ARRAY class for use in fftw++
#-------------------------------------------------------------------------------
$(ARRAYDIST):
	$(CURL) "$(ARRAYURL)" > $(ARRAYDIST)

$(ARRAYBUILDDATE):	$(ARRAYDIST)
	@echo "################################################################################"
	@echo "# Building $(ARRAYDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	$(CP) -f $(ARRAYDIST) $(SPHERALDIR)/include
	@if ! test -e $(SPHERALDIR)/include/$(ARRAYDIST); then \
	echo "ERROR building ARRAY."; \
	exit 1; \
	fi
	touch $(ARRAYBUILDDATE)

#-------------------------------------------------------------------------------
# FFTW++ interface
#-------------------------------------------------------------------------------
$(FFTWPPDIST):
	$(CURL) "$(FFTWPPURL)" > $(FFTWPPDIST)

$(FFTWPPBUILDDATE):	$(FFTWPPDIST) $(ARRAYBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(FFTWPPDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(FFTWPPDIR)
	$(GZIP) -dc $(FFTWPPDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(FFTWPPDIR); \
	$(CXX) $(CXXFLAGS) -c -o fftw++.o fftw++.cc; \
	$(AR) rcs libfftw++.a fftw++.o'
	$(CP) -f $(FFTWPPDIR)/fftw++.h $(SPHERALDIR)/include
	$(CP) -f $(FFTWPPDIR)/libfftw++.a $(SPHERALDIR)/lib
	@if ! test -e $(SPHERALDIR)/include/fftw++.h -a \
	           -e $(SPHERALDIR)/lib/libfftw++.a; then \
	echo "ERROR building FFTW++."; \
	exit 1; \
	fi
	rm -fr $(FFTWPPBUILDDATE) $(FFTWPPDIR)
	touch $(FFTWPPBUILDDATE)

#-------------------------------------------------------------------------------
# ZLIB
#-------------------------------------------------------------------------------
$(ZLIBDIST):
	$(CURL) "$(ZLIBURL)" > $(ZLIBDIST)

$(ZLIBBUILDDATE):	$(ZLIBDIST)
	@echo "################################################################################"
	@echo "# Building $(ZLIBDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(ZLIBDIR)
	$(GZIP) -dc $(ZLIBDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(ZLIBDIR); \
	./configure --prefix=$(SPHERALDIR) ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/include/zlib.h -a \
	           -e $(SPHERALDIR)/include/zconf.h -a \
	           -e $(SPHERALDIR)/lib/libz.a; then \
	echo "ERROR building ZLIB."; \
	exit 1; \
	fi
	rm -fr $(ZLIBBUILDDATE) $(ZLIBDIR)
	touch $(ZLIBBUILDDATE)

#-------------------------------------------------------------------------------
# JPEG
#-------------------------------------------------------------------------------
$(JPEGDIST):
	$(CURL) "$(JPEGURL)" > $(JPEGDIST)

$(JPEGBUILDDATE):	$(JPEGDIST)
	@echo "################################################################################"
	@echo "# Building $(JPEGDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(JPEGDIR)
	$(GZIP) -dc $(JPEGDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(JPEGDIR); \
	env CC='$(CC)' CXX='$(CXX)' CFLAGS="-fPIC" ./configure --prefix=$(SPHERALDIR) ; \
	$(MAKE) ; \
	$(MAKE) install ; \
	$(CP) -f jpeglib.h jerror.h jconfig.h jmorecfg.h $(SPHERALDIR)/include ; \
	$(CP) -f .libs/libjpeg.a $(SPHERALDIR)/lib"
	@if ! test -e $(SPHERALDIR)/include/jpeglib.h -a \
	           -e $(SPHERALDIR)/include/jerror.h -a \
	           -e $(SPHERALDIR)/include/jconfig.h -a \
	           -e $(SPHERALDIR)/include/jmorecfg.h -a \
	           -e $(SPHERALDIR)/lib/libjpeg.a; then \
	echo "ERROR building JPEG."; \
	exit 1; \
	fi
	rm -fr $(JPEGBUILDDATE) $(JPEGDIR)
	touch $(JPEGBUILDDATE)

#-------------------------------------------------------------------------------
# PIL (Python Imaging Library)
#-------------------------------------------------------------------------------
$(PILDIST):
	$(CURL) "$(PILURL)" > $(PILDIST)

$(PILBUILDDATE):	$(PYTHONBUILDDATE) $(PILDIST) $(JPEGBUILDDATE) $(ZLIBBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(PILDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(PILDIR)
	$(GZIP) -dc $(PILDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(PILDIR); $(PYTHONEXE) setup.py install;'
	@if test -n "`$(PYTHONEXE) -c 'import PIL'`"; then \
	echo "ERROR building PIL python extension."; \
	rm -f $(PILBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PILDIR)
	touch $(PILBUILDDATE)

#-------------------------------------------------------------------------------
# pybindgen (Python bindings)
#-------------------------------------------------------------------------------
$(PYBINDGENDIST):
	$(CURL) "$(PYBINDGENURL)" > $(PYBINDGENDIST)

$(PYBINDGENBUILDDATE):	$(PYTHONBUILDDATE) $(PYBINDGENDIST)
	@echo "################################################################################"
	@echo "# Building $(PYBINDGENDIST)"
	@echo "################################################################################"
	rm -fr $(PYBINDGENDIR)
	$(GZIP) -dc $(PYBINDGENDIST) | $(TAR) -xf -
	$(CP) -r $(PYBINDGENDIR)/pybindgen $(PYTHONLIBDIR)
	@if test -n "`$(PYTHONEXE) -c 'import pybindgen'`"; then \
	echo "ERROR building PYBINDGEN python extension."; \
	rm -f $(PYBINDGENBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PYBINDGENDIR)
	touch $(PYBINDGENBUILDDATE)

#-------------------------------------------------------------------------------
# WildMagic
#-------------------------------------------------------------------------------
$(WILDMAGICDIST):
	$(CURL) "$(WILDMAGICURL)" > $(WILDMAGICDIST)

$(WILDMAGICBUILDDATE):	$(PYTHONBUILDDATE) $(WILDMAGICDIST)
	@echo "################################################################################"
	@echo "# Building $(WILDMAGICDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(WILDMAGICDIR)
	$(UNZIP) $(WILDMAGICDIST)
	rm -fr $(SPHERALDIR)/WildMagic5
	$(MKDIR) $(SPHERALDIR)/WildMagic5
	chmod +x $(BUILDWILDMAGIC) $(WILDMAGICDIR)/WildMagic5/*.sh $(WILDMAGICDIR)/WildMagic5/*/*.sh
	$(CP) $(BUILDWILDMAGIC) GeometricTools/WildMagic5
	$(CP) -f makefile.wm5 GeometricTools/WildMagic5
	$(SHELL) -ec 'env MAKE=$(MAKE) CXX=$(CXX) ./$(BUILDWILDMAGIC) $(WILDMAGICTARGET) build $(MACHTYPE)'
	$(CP) -r $(WILDMAGICDIR)/WildMagic5/SDK/* $(SPHERALDIR)/WildMagic5
	ln -f -s $(WILDMAGICLIBDIR)/libWm5Core.$(WMLIBEXT) $(LIBDIR)/libWm5Core.$(WMLIBEXT)
	ln -f -s $(WILDMAGICLIBDIR)/libWm5Mathematics.$(WMLIBEXT)  $(LIBDIR)/libWm5Mathematics.$(WMLIBEXT) 
	ln -f -s $(WILDMAGICLIBDIR)/libWm5Imagics.$(WMLIBEXT) $(LIBDIR)/libWm5Imagics.$(WMLIBEXT) 
	ln -f -s $(WILDMAGICLIBDIR)/libWm5Physics.$(WMLIBEXT) $(LIBDIR)/libWm5Physics.$(WMLIBEXT) 
	@if ! test -e $(WILDMAGICLIBDIR)/libWm5Core.$(WMLIBEXT) -a \
		   -e $(WILDMAGICLIBDIR)/libWm5Mathematics.$(WMLIBEXT) -a \
		   -e $(WILDMAGICLIBDIR)/libWm5Imagics.$(WMLIBEXT) -a \
		   -e $(WILDMAGICLIBDIR)/libWm5Physics.$(WMLIBEXT); then \
	echo "ERROR building WildMagic."; \
	rm -f $(WILDMAGICBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(WILDMAGICDIR)
	touch $(WILDMAGICBUILDDATE)

#-------------------------------------------------------------------------------
# Qhull
#-------------------------------------------------------------------------------
$(QHULLDIST):
	$(CURL) "$(QHULLURL)" > $(QHULLDIST)

$(QHULLDATE):	$(PYTHONBUILDDATE) $(QHULLDIST) $(CMAKEEXE)
	@echo "################################################################################"
	@echo "# Building $(QHULLDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(QHULLDIR) $(QHULLINCDIR)
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(QHULLINCDIR)
	$(GZIP) -dc $(QHULLDIST) | $(TAR) -xf -
	$(PATCH) $(QHULLDIR)/src/libqhull/qhull_a.h $(QHULLDIR)-qhull_a.h-patch
	$(SHELL) -ec "cd $(QHULLDIR)/build; \
	env CC='$(CC)' CXX='$(CXX)' CFLAGS='$(CFLAGS)' $(CMAKEEXE) -DCMAKE_INSTALL_PREFIX:PATH=$(SPHERALDIR) .. ; \
	$(MAKE); \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libqhullstatic.a; then \
	echo "ERROR building Qhull."; \
	rm -f $(QHULLDATE); \
	exit 1; \
	fi
	rm -fr $(QHULLDIR)
	touch $(QHULLDATE)

#-------------------------------------------------------------------------------
# mpi4py
#-------------------------------------------------------------------------------
$(MPI4PYDIST):
	$(CURL) "$(MPI4PYURL)" > $(MPI4PYDIST)

$(MPI4PYBUILDDATE):  $(PYTHONBUILDDATE) $(MPI4PYDIST)
	@echo "################################################################################"
	@echo "# Building $(MPI4PYDIST)"
	@echo "################################################################################"
	rm -fr $(MPI4PYDIR)
	$(GZIP) -dc $(MPI4PYDIST) | $(TAR) -xf -
	$(PATCH) $(MPI4PYDIR)/setup.cfg mpi4py-1.2.2_setup_dot_cfg.patch
	$(SHELL) -ec "cd $(MPI4PYDIR); $(PYTHONEXE) setup.py install"
	@if test -n "`$(PYTHONEXE) -c 'from mpi4py import MPI'`"; then \
	echo "ERROR building mpi4py python extension."; \
	rm -f $(MPI4PYBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(MPI4PYDIR)
	touch $(MPI4PYBUILDDATE)

#-------------------------------------------------------------------------------
# HDF5
#-------------------------------------------------------------------------------
$(HDF5DIST):
	$(CURL) "$(HDF5URL)" > $(HDF5DIST)

$(HDF5BUILDDATE):  $(HDF5DIST)
	@echo "################################################################################"
	@echo "# Building $(HDF5DIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(HDF5DIR)
	$(BZIP2) -dc $(HDF5DIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(HDF5DIR); \
	env CC='$(CC)' CXX='$(CXX)' CFLAGS='$(CFLAGS)' CXXFLAGS='$(CXXFLAGS)' ./configure --prefix=$(SPHERALDIR) $(HDF5FLAGS) ; \
	rm -f test.conf ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libhdf5.a -a \
	           -e $(SPHERALDIR)/include/hdf5.h; then \
	echo "ERROR building HDF5."; \
	rm -f $(HDF5BUILDDATE); \
	exit 1; \
	fi
	rm -fr $(HDF5DIR)
	touch $(HDF5BUILDDATE)

#-------------------------------------------------------------------------------
# SILO
#-------------------------------------------------------------------------------
$(SILODIST):
	$(CURL) -k "$(SILOURL)" > $(SILODIST)

$(SILOBUILDDATE):  $(SILODIST) $(HDF5BUILDDATE) $(PYTHONBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(SILODIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(SILODIR)
	$(GZIP) -dc $(SILODIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(SILODIR); \
	env PATH=$(SPHERALDIR)/bin:$$PATH CC='$(MPICC)' CXX='$(MPICXX)' CFLAGS='$(MPICCFLAGS)' CXXFLAGS='$(MPICXXFLAGS)' $(SILOFLAGS) ./configure --enable-fortran=no --enable-shared --prefix=$(SPHERALDIR) --with-hdf5='$(SPHERALDIR)/include,$(SPHERALDIR)/lib' --enable-pythonmodule --enable-silex=no --enable-browser=no ; \
	$(MAKE) ; \
	$(MAKE) install"
	ln -f -s $(SPHERALDIR)/lib/Silo.so $(LIBDIR)
	@if ! test -e $(SPHERALDIR)/lib/libsiloh5.a -a \
	           -e $(SPHERALDIR)/lib/libsiloh5.$(DYLIBEXT); then \
	echo "ERROR building SILO."; \
	rm -f $(SILOBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(SILODIR)
	touch $(SILOBUILDDATE)

#-------------------------------------------------------------------------------
# Eigen
#-------------------------------------------------------------------------------
$(EIGENDIST):
	$(CURL) -k "$(EIGENURL)" > $(EIGENDIST)

$(EIGENBUILDDATE):  $(EIGENDIST) $(HDF5BUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(EIGENDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(EIGENDIR)
	$(BZIP2) -dc $(EIGENDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(EIGENDIR); $(CP) -r Eigen $(SPHERALDIR)/include'
	@if ! test -e $(SPHERALDIR)/include/Eigen; then \
	echo "ERROR building EIGEN."; \
	rm -f $(EIGENBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(EIGENDIR)
	touch $(EIGENBUILDDATE)

#-------------------------------------------------------------------------------
# Voro++
#-------------------------------------------------------------------------------
$(VOROPPDIST):
	$(CURL) -k "$(VOROPPURL)" > $(VOROPPDIST)

$(VOROPPBUILDDATE):  $(VOROPPDIST)
	@echo "################################################################################"
	@echo "# Building $(VOROPPDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	@$(MKDIR) $(SPHERALDIR)/include/voro++
	rm -fr $(VOROPPDIR)
	$(GZIP) -dc $(VOROPPDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(VOROPPDIR)/src; $(CP) -r * $(SPHERALDIR)/include/voro++'
	@if ! test -e $(SPHERALDIR)/include/voro++/voro++.hh; then \
	echo "ERROR building voro++."; \
	rm -f $(VOROPPBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(VOROPPDIR)
	touch $(VOROPPBUILDDATE)

#-------------------------------------------------------------------------------
# Voro++ (2D)
#-------------------------------------------------------------------------------
$(VOROPP2DBUILDDATE):  $(VOROPP2DDIST)
	@echo "################################################################################"
	@echo "# Building $(VOROPP2DDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(VOROPP2DDIR) $(SPHERALDIR)/include/$(VOROPP2DINC)
	$(MKDIR) $(SPHERALDIR)/include/$(VOROPP2DINC)
	$(BZIP2) -dc $(VOROPP2DDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(VOROPP2DDIR)/src; $(MAKE) CC='$(CC)' CXX='$(CXX)' CFLAGS="-fPIC"; \
	$(CP) -r *.hh $(SPHERALDIR)/include/$(VOROPP2DINC); \
	$(CP) *.a  $(SPHERALDIR)/lib'
	@if ! test -e $(SPHERALDIR)/include/$(VOROPP2DINC)/voro++_2d.hh -a \
	           -e $(SPHERALDIR)/lib/libvoro++_2d.a; then \
	echo "ERROR building voro++ 2D."; \
	rm -f $(VOROPP2DBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(VOROPP2DDIR)
	touch $(VOROPP2DBUILDDATE)

#-------------------------------------------------------------------------------
# Gnu multiprecision library
#-------------------------------------------------------------------------------
$(GMPDIST):
	$(CURL) -k "$(GMPURL)" > $(GMPDIST)

$(GMPBUILDDATE):  $(GMPDIST)
	@echo "################################################################################"
	@echo "# Building $(GMPDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(GMPDIR)
	$(BZIP2) -dc $(GMPDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(GMPDIR); \
	env CC="$(CC)" CXX="$(CXX)" CFLAG="$(CFLAGS)" CXXFLAG="$(CXXFLAGS)" ./configure --prefix=$(SPHERALDIR) ;\
	$(MAKE); \
	$(MAKE) check; \
	$(MAKE) install'
	@if ! test -e $(SPHERALDIR)/include/gmp.h -a \
	           -e $(SPHERALDIR)/lib/libgmp.a; then \
	echo "ERROR building gmp."; \
	rm -f $(GMPBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(GMPDIR)
	touch $(GMPBUILDDATE)

#-------------------------------------------------------------------------------
# mpmath
#-------------------------------------------------------------------------------
$(MPMATHDIST):
	$(CURL) "$(MPMATHURL)" > $(MPMATHDIST)

$(MPMATHBUILDDATE):  $(PYTHONBUILDDATE) $(MPMATHDIST)
	@echo "################################################################################"
	@echo "# Building $(MPMATHDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/include $(SPHERALDIR)/bin
	rm -fr $(MPMATHDIR)
	$(GZIP) -dc $(MPMATHDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(MPMATHDIR); $(PYTHONEXE) setup.py install"
	@if test -n "`$(PYTHONEXE) -c 'import mpmath'`"; then \
	echo "ERROR building mpmath python extension."; \
	rm -f $(MPMATHBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(MPMATHDIR)
	touch $(MPMATHBUILDDATE)

#-------------------------------------------------------------------------------
# The HYPRE linear solver library
#-------------------------------------------------------------------------------
$(HYPREDIST):
	$(CURL) "$(HYPREURL)" > $(HYPREDIST)

$(HYPREBUILDDATE):	$(HYPREDIST)
	@echo "################################################################################"
	@echo "# Building $(HYPREDIST)"
	@echo "################################################################################"
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(SPHERALDIR)/include
	rm -fr $(HYPREDIR)
	$(GZIP) -dc $(HYPREDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(HYPREDIR)/src; \
	env CC='$(MPICC)' CXX='$(MPICXX)' CFLAGS='$(MPICCFLAGS) $(CFLAGS)' CXXFLAGS='$(MPICXXFLAGS) $(CXXFLAGS)' ./configure --prefix=$(SPHERALDIR) ;\
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libHYPRE.a -a \
	           -e $(SPHERALDIR)/include/HYPRE.h; then \
	echo "ERROR building Hypre."; \
	rm -f $(HYPREBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(HYPREDIR)
	touch $(HYPREBUILDDATE)

#-------------------------------------------------------------------------------
# FLANN
#-------------------------------------------------------------------------------
$(FLANNDIST):
	$(CURL) "$(FLANNURL)" > $(FLANNDIST)

$(FLANNBUILDDATE):	$(FLANNDIST) $(CMAKEEXE) $(FLANNDATE) $(PYTHONBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(FLANNDIST)"
	@echo "################################################################################"
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(SPHERALDIR)/include
	rm -fr $(FLANNDIR)
	$(UNZIP) $(FLANNDIST)
	$(MKDIR) $(FLANNDIR)/build
	$(SHELL) -ec "cd $(FLANNDIR)/build; \
	env CC='$(CC)' CXX='$(CXX)' CFLAGS='$(CFLAGS)' $(CMAKEEXE) -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=$(SPHERALDIR) .. ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libflann_s.a; then \
	echo "ERROR building Flann."; \
	rm -f $(FLANNBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(FLANNDIR)/build
	touch $(FLANNBUILDDATE)

#-------------------------------------------------------------------------------
# The Point Cloud Library (PCL)
#-------------------------------------------------------------------------------
$(PCLDIST):
	$(CURL) "$(PCLURL)" > $(PCLDIST)

$(PCLBUILDDATE):	$(PCLDIST) $(CMAKEEXE) $(BOOSTDATE) $(EIGENDATE) $(QHULLDATE) $(FLANNDATE)
	@echo "################################################################################"
	@echo "# Building $(PCLDIST)"
	@echo "################################################################################"
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(SPHERALDIR)/include
	rm -fr $(PCLDIR)
	$(BZIP2) -dc $(PCLDIST) | $(TAR) -xf -
	$(MKDIR) $(PCLDIR)/build
	$(SHELL) -ec "cd $(PCLDIR)/build; \
	env CC='$(CC)' CXX='$(CXX)' CFLAGS='$(CFLAGS)' $(CMAKEEXE) -DFLANN_INCLUDE_DIR=$(SPHERALDIR)/spheral/src/thirdPartyLibs/$(FLANNDIR) -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=$(SPHERALDIR) .. ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libpclstatic.a; then \
	echo "ERROR building Pcl."; \
	rm -f $(PCLBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PCLDIR)
	touch $(PCLBUILDDATE)

#-------------------------------------------------------------------------------
# Octtree
#-------------------------------------------------------------------------------
$(OCTTREEDIST):
	$(CURL) "$(OCTTREEURL)" > $(OCTTREEDIST)

$(OCTTREEBUILDDATE):	$(OCTTREEDIST)
	@echo "################################################################################"
	@echo "# Building $(OCTTREEDIST)"
	@echo "################################################################################"
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(OCTTREEINCDIR)
	rm -fr $(OCTTREEDIR) $(OCTTREEINCDIR)
	$(MKDIR) $(OCTTREEDIR)
	$(SHELL) -ec "cd $(OCTTREEDIR); \
	$(UNZIP) ../$(OCTTREEDIST)"
	$(CP) -f $(OCTTREEDIR)/*.h $(OCTTREEDIR)/*.tcc $(OCTTREEINCDIR)
	@if ! test -e $(OCTTREEINCDIR)/octree.h; then \
	echo "ERROR building Octtree."; \
	rm -f $(OCTTREEBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(OCTTREEDIR)
	touch $(OCTTREEBUILDDATE)

#-------------------------------------------------------------------------------
# polytope
#-------------------------------------------------------------------------------
$(TRIANGLEDIST):
	$(CURL) "$(TRIANGLEURL)" > $(TRIANGLEDIST)

$(TETGENDIST):
	$(CURL) "$(TETGENURL)" > $(TETGENDIST)

$(POLYTOPEDIST):
	$(CURL) "$(POLYTOPEURL)" > $(POLYTOPEDIST)

$(POLYTOPEBUILDDATE):	$(POLYTOPEDIST) $(CMAKEEXE) $(BOOSTDATE) $(SILOBUILDDATE) $(PYTHONBUILDDATE) $(PYBINDGENBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(POLYTOPEDIST)"
	@echo "################################################################################"
	rm -fr $(POLYTOPEDIR) $(POLYTOPEINCDIR) $(TRIANGLEDIR)
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(POLYTOPEINCDIR) $(TRIANGLEDIR)
	$(BZIP2) -dc $(POLYTOPEDIST) | $(TAR) -xf -
	@if test "$(USE_TRIANGLE)" = "1"; then \
		$(MAKE) $(TRIANGLEDIST); \
		$(SHELL) -ec "cd $(TRIANGLEDIR); $(UNZIP) ../$(TRIANGLEDIST)"; \
		$(CP) -f $(TRIANGLEDIR)/triangle.h $(POLYTOPEDIR)/src; \
		$(CP) -f $(TRIANGLEDIR)/triangle.c $(POLYTOPEDIR)/src; \
	fi
	@if test "$(USE_TETGEN)" = "1"; then \
		$(MAKE) $(TETGENDIST); \
		$(GZIP) -cd $(TETGENDIST) | $(TAR) -xf - ; \
		$(CP) -f $(TETGENDIR)/tetgen.h $(POLYTOPEDIR)/src; \
		$(CP) -f $(TETGENDIR)/tetgen.cxx $(POLYTOPEDIR)/src; \
	fi
	env PATH=$(SPHERALDIR)/bin:$(PATH) $(SHELL) -ec "cd $(POLYTOPEDIR); \
	$(MAKE) MAKEFLAGS= $(POLYTOPEFLAGS) config; \
	$(MAKE) MAKEFLAGS= $(POLYTOPEFLAGS) install"
	@if ! test -e $(SPHERALDIR)/lib/libpolytope.a -a \
	           -e $(SPHERALDIR)/lib/libvoro_2d.a -a \
	           -e $(SPHERALDIR)/lib/libvoro_3d.a -a \
                   -e $(SPHERALDIR)/lib/PolytopeModules.so -a \
	           -e $(SPHERALDIR)/include/polytope/polytope.hh; then \
	echo "ERROR building Polytope."; \
	rm -f $(POLYTOPEBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(POLYTOPEDIR) $(TRIANGLEDIR) $(TETGENDIR)
	ln -f -s $(SPHERALDIR)/lib/PolytopeModules.$(DYLIBEXT) $(LIBDIR)/PolytopeModules.$(DYLIBEXT)
	touch $(POLYTOPEBUILDDATE)

#-------------------------------------------------------------------------------
# OpenSubDiv
#-------------------------------------------------------------------------------
$(OPENSUBDIVDIST):
	$(CURL) "$(OPENSUBDIVURL)" > $(OPENSUBDIVDIST)

$(OPENSUBDIVBUILDDATE):	$(OPENSUBDIVDIST) $(CMAKEEXE)
	@echo "################################################################################"
	@echo "# Building $(OPENSUBDIVDIST)"
	@echo "################################################################################"
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(SPHERALDIR)/include
	rm -fr $(OPENSUBDIVDIR)
	$(UNZIP) $(OPENSUBDIVDIST)
	$(PATCH) $(OPENSUBDIVDIR)/CMakeLists.txt OpenSubdiv-2_5_0-CMakeLists_txt.patch
	$(PATCH) $(OPENSUBDIVDIR)/opensubdiv/osd/CMakeLists.txt OpenSubdiv-2_5_0-opensubdiv-osd-CMakeLists_txt.patch
	$(MKDIR) $(OPENSUBDIVDIR)/build
	env PATH=$(SPHERALDIR)/bin:$(PATH) $(SHELL) -ec "cd $(OPENSUBDIVDIR)/build; \
	$(CMAKEEXE) -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(SPHERALDIR) -DCMAKE_C_COMPILER='$(CC)' -DCMAKE_CXX_COMPILER='$(CXX)' -DCMAKE_C_FLAGS='$(CFLAGS)' -DCMAKE_CXX_FLAGS='$(CXXFLAGS)' -DNO_PYTHON=1 -DNO_EXAMPLES=1 -DNO_REGRESSION=1 -DNO_OMP=1 -DNO_CUDA=1 -DNO_OPENGL=1 -DNO_CLEW=1 .. ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libosdCPU.a -a \
	           -e $(SPHERALDIR)/lib/libosdutil.a ; then \
	echo "ERROR building Opensubdiv."; \
	rm -f $(OPENSUBDIVBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(OPENSUBDIVDIR)
	touch $(OPENSUBDIVBUILDDATE)

