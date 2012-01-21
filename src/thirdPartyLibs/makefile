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
CURL = curl -L

#-------------------------------------------------------------------------------
# Targets.
#-------------------------------------------------------------------------------
# Python
PYTHONDIST = Python-2.7.2.tar.bz2
PYTHONURL = http://www.python.org/ftp/python/2.7.2/Python-2.7.2.tar.bz2
PYTHONDIR = $(PYTHONDIST:.tar.bz2=)
PYTHONBUILDDATE = $(patsubst %.tar.bz2, .%.date, $(PYTHONDIST))

# Numpy
NUMPYDIST = numpy-1.5.1.tar.gz
NUMPYURL = http://sourceforge.net/projects/numpy/files/NumPy/1.5.1/$(NUMPYDIST)/download
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

# Element tree
ELEMENTTREEDIST = elementtree-1.2.7-20070827-preview.zip
ELEMENTTREEURL = http://effbot.org/media/downloads/$(ELEMENTTREEDIST)
ELEMENTTREEDIR = $(ELEMENTTREEDIST:.zip=)
ELEMENTTREEBUILDDATE = $(patsubst %.zip, .%.date, $(ELEMENTTREEDIST))

# Boost
BOOSTDIST = boost_$(BOOSTVER1).tar.bz2
BOOSTURL = http://sourceforge.net/projects/boost/files/boost/$(BOOSTVER2)/$(BOOSTDIST)/download
BOOSTDIR = $(BOOSTDIST:.tar.bz2=)
BOOSTBUILDDATE = $(patsubst %.tar.bz2, .%.date, $(BOOSTDIST))
BOOSTPYTHONLIB0 = libboost_python$(BOOSTEXT).$(DYLIBEXT)
BOOSTPYTHONLIB = libboost_python.$(DYLIBEXT)

# Cmake
CMAKEEXE = $(SPHERALDIR)/bin/cmake
CMAKEDIST = cmake-2.8.7.tar.gz
CMAKEURL = http://www.cmake.org/files/v2.8/$(CMAKEDIST)
CMAKEDIR = $(CMAKEDIST:.tar.gz=)

# Gccxml
GCCXMLDIR = $(GCCXMLDIST:.tar.bz2=)
GCCXMLDATE = $(patsubst %.tar.bz2, .%.date, $(GCCXMLDIST))

# Pympi
PYMPIDIST = pyMPI-2.4b4.tar.bz2
PYMPIURL = http://downloads.sourceforge.net/pympi/pyMPI-2.4b4.tar.gz?modtime=1138016139&big_mirror=0
PYMPIDIR = $(patsubst %.tar.bz2, %, $(PYMPIDIST))
PYMPIBUILDDATE = $(patsubst %.tar.bz2, .%.date, $(PYMPIDIST))

# Parmetis
PARMETISDIST = parmetis-4.0.tar.gz
PARMETISURL = http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/$(PARMETISDIST)
PARMETISDIR = $(PARMETISDIST:.tar.gz=)
PARMETISBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PARMETISDIST))
#PARMETISMAKEFILE = Makefile.ParMETIS

# Tau
TAUURL=http://tau.uoregon.edu/$(TAUDIST)
TAUDIR = $(TAUDIST:.tgz=-2.21.1)
TAUBUILDDATE = $(patsubst %.tgz, .%.date, $(TAUDIST))
TAUDUMMYDATE = $(patsubst %.tgz, .%.dummydate, $(TAUDIST))

# Pygtk
PYGTKDIST = pygtk-2.12.1.tar.gz
PYGTKURL = http://saimei.acc.umu.se/pub/GNOME/sources/pygtk/2.12/$(PYGTKDIST)
PYGTKDIR = $(PYGTKDIST:.tar.gz=)
PYGTKBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PYGTKDIST))

# Matplotlib
MATPLOTLIBDIST = matplotlib-1.0.0.tar.gz
MATPLOTLIBURL = https://sourceforge.net/projects/matplotlib/files/matplotlib-1.0/$(MATPLOTLIBDIST)/download
MATPLOTLIBDIR = $(MATPLOTLIBDIST:.tar.gz=)
MATPLOTLIBBUILDDATE = $(patsubst %.tar.gz, .%.date, $(MATPLOTLIBDIST))

# Scipy
SCIPYDIST = scipy-0.6.0.tar.gz
SCIPYURL = http://easynews.dl.sourceforge.net/sourceforge/scipy/$(SCIPYDIST)
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
VISITWRITERDIR = $(VISITWRITERDIST:.tar.gz=)
VISITWRITERBUILDDATE = $(patsubst %.tar.gz, .%.date, $(VISITWRITERDIST))

# ATS testing system.
ATSDIST = ats-5.2.tar.gz
ATSURL = http://ats.googlecode.com/files/$(ATSDIST)
ATSDIR = $(ATSDIST:.tar.gz=)
ATSBUILDDATE = $(patsubst %.tar.gz, .%.date, $(ATSDIST))

# FFTW
FFTWDIST = fftw-3.2.2.tar.gz
FFTWURL = http://www.fftw.org/fftw-3.2.2.tar.gz
FFTWDIR = $(FFTWDIST:.tar.gz=)
FFTWBUILDDATE = $(patsubst %.tar.gz, .%.date, $(FFTWDIST))

# FFTW++
FFTWPPDIST = fftw++-1.08.tar.gz
FFTWPPURL = http://sourceforge.net/projects/fftwpp/files/$(FFTWPPDIST)/download
FFTWPPDIR = $(FFTWPPDIST:.tar.gz=)
FFTWPPBUILDDATE = $(patsubst %.tar.gz, .%.date, $(FFTWPPDIST))

# Array
ARRAYDIST = Array.h
ARRAYURL = http://www.math.ualberta.ca/~bowman/$(ARRAYDIST)
ARRAYBUILDDATE = $(patsubst %.h, .%.date, $(ARRAYDIST))

# ZLIB
ZLIBDIST = zlib-1.2.5.tar.bz2
ZLIBURL = http://www.zlib.net/$(ZLIBDIST)
ZLIBDIR = $(ZLIBDIST:.tar.bz2=)
ZLIBBUILDDATE = $(patsubst %.tar.bz2, .%.date, $(ZLIBDIST))

# JPEG
JPEGDIST = jpegsrc.v7.tar.gz
JPEGURL = http://www.ijg.org/files/$(JPEGDIST)
JPEGDIR = jpeg-7
JPEGBUILDDATE = $(patsubst %.tar.gz, .%.date, $(JPEGDIST))

# PIL
PILDIST = Imaging-1.1.6.tar.gz
PILURL = http://effbot.org/downloads/$(PILDIST)
PILDIR = $(PILDIST:.tar.gz=)
PILBUILDDATE = $(patsubst %.tar.gz, .%.date, $(PILDIST))

# pybindgen
PYBINDGENDIST = pybindgen-0.15.1.zip
PYBINDGENURL = http://pybindgen.googlecode.com/files/$(PYBINDGENDIST)
PYBINDGENDIR = $(PYBINDGENDIST:.zip=)
PYBINDGENBUILDDATE = $(patsubst %.zip, .%.date, $(PYBINDGENDIST))

# WildMagic
WILDMAGICDIST = WildMagic5p7.zip
WILDMAGICURL = http://www.geometrictools.com/Downloads/$(WILDMAGICDIST)
WILDMAGICDIR = GeometricTools
WILDMAGICBUILDDATE = $(patsubst %.zip, .%.date, $(WILDMAGICDIST))

# Qhull
QHULLDIST = qhull-2011.2-src.tgz
QHULLURL = http://www.qhull.org/download/$(QHULLDIST)
QHULLDIR = $(QHULLDIST:-src.tgz=)
QHULLINCDIR = $(SPHERALDIR)/include/qhull
QHULLBUILDDATE = $(patsubst %-src.tgz, .%.date, $(QHULLDIST))

# mpi4py
MPI4PYDIST = mpi4py-1.2.2.tar.gz
MPI4PYURL = http://mpi4py.googlecode.com/files/$(MPI4PYDIST)
MPI4PYDIR = $(MPI4PYDIST:.tar.gz=)
MPI4PYBUILDDATE = $(patsubst %.tar.gz, .%.date, $(MPI4PYDIST))

# HDF5
HDF5DIST = hdf5-1.8.8.tar.bz2
HDF5URL = http://www.hdfgroup.org/ftp/HDF5/current/src/$(HDF5DIST)
HDF5DIR = $(HDF5DIST:.tar.bz2=)
HDF5BUILDDATE = $(patsubst %.tar.bz2, .%.date, $(HDF5DIST))

# SILO
SILODIST = silo-4.8-bsd.tar.gz
SILOURL = https://wci.llnl.gov/codes/silo/silo-4.8/$(SILODIST)
SILODIR = $(SILODIST:.tar.gz=)
SILOBUILDDATE = $(patsubst %.tar.gz, .%.date, $(SILODIST))

# Eigen
EIGENDIST = 3.0.2.tar.bz2
EIGENURL = http://bitbucket.org/eigen/eigen/get/$(EIGENDIST)
EIGENDIR = eigen-eigen-3.0.2
EIGENBUILDDATE = .eigen.date

# Voro++
VOROPPDIST = voro++_0.3.1.tar.gz
VOROPPURL = http://math.lbl.gov/voro++/download/dir/$(VOROPPDIST)
VOROPPDIR = voro++
VOROPPBUILDDATE = $(patsubst %.tar.gz, .%.date, $(VOROPPDIST))

# Voro++ 2d specialization
VOROPP2DDIST = voro++_2d.tbz
VOROPP2DDIR = voro++_2d
VOROPP2DINC = voro++_2d
VOROPP2DBUILDDATE = $(patsubst %.tbz, .%.date, $(VOROPP2DDIST))

# GNU multiprecision library.
GMPDIST = gmp-5.0.2.tar.bz2
GMPURL = ftp://ftp.gmplib.org/pub/gmp-5.0.2/$(GMPDIST)
GMPDIR = $(GMPDIST:.tar.bz2=)
GMPBUILDDATE = $(patsubst %.tar.bz2, .%.date, $(GMPDIST))

#-------------------------------------------------------------------------------
all:	$(PYTHONBUILDDATE) \
	$(ELEMENTTREEBUILDDATE) \
	$(BOOSTBUILDDATE) \
	$(GCCXMLTARGETS) \
	$(NUMPYTHIRDPARTYTARGETS) \
	$(PSYCOTHIRDPARTYTARGET) \
	$(MPITHIRDPARTYTARGETS) \
	$(MATPLOTLIBTARGETS) \
	$(SCIPYTARGETS) \
	$(GSLBUILDDATE) \
	$(PYGSLBUILDDATE) \
	$(MHDTHIRDPARTYTARGETS) \
	$(PETSCTARGETS) \
	$(TAUTARGET) \
	$(VISITWRITERBUILDDATE) \
	$(ATSBUILDDATE) \
	$(FFTWBUILDDATE) \
	$(FFTWPPBUILDDATE) \
	$(ARRAYBUILDDATE) \
	$(ZLIBBUILDDATE) \
	$(JPEGBUILDDATE) \
	$(PILBUILDDATE) \
	$(PYBINDGENBUILDDATE) \
	$(WILDMAGICBUILDDATE) \
	$(QHULLBUILDDATE) \
	$(HDF5BUILDDATE) \
	$(SILOBUILDDATE) \
	$(EIGENBUILDDATE) \
	$(VOROPPBUILDDATE) \
	$(VOROPP2DBUILDDATE)

clean:
	rm -fr $(PYTHONDIR) $(PYTHONBUILDDATE) \
	       $(NUMPYDIR) $(NUMPYBUILDDATE) \
	       $(GNUPLOTDIR) $(GNUPLOTBUILDDATE) \
	       $(ELEMENTTREEDIR) $(ELEMENTTREEBUILDDATE) \
	       $(PSYCODIR) $(PSYCOBUILDDATE) \
	       $(ELEMENTTREEDIR) $(ELEMENTTREEBUILDDATE) \
	       $(BOOSTDIR) $(BOOSTBUILDDATE) \
               $(CMAKEDIR) \
               $(GCCXMLDIR) $(GCCXMLDATE) \
	       $(PSYCODIR) $(PSYCOTHIRDPARTYTARGET) \
               $(PYMPIDIR) $(PYMPIBUILDDATE) \
	       $(MPI4PYDIR) $(MPI4PYBUILDDATE) \
               $(PARMETISDIR) $(PARMETISBUILDDATE) \
	       $(SCIPYDIR) $(SCIPYBUILDDATE) \
	       $(GSLDIR) $(GSLBUILDDATE) \
	       $(PYGSLDIR) $(PYGSLBUILDDATE) \
               $(MATPLOTLIBDIR) $(MATPLOTLIBBUILDDATE) \
	       $(TAUDIR) $(TAUBUILDDATE) $(TAUDUMMYDATE) \
	       $(HASHMAPDATE) \
	       $(VISITWRITERDIR) $(VISITWRITERBUILDDATE) \
	       $(ATSWRITERDIR) $(ATSWRITERBUILDDATE) \
	       $(FFTWDIR) $(FFTWBUILDDATE) \
	       $(FFTWPPDIR) $(FFTWPPBUILDDATE) \
	       $(ARRAYDIR) $(ARRAYBUILDDATE) \
	       $(ZLIBDIR) $(ZLIBBUILDDATE) \
	       $(JPEGDIR) $(JPEGBUILDDATE) \
	       $(PILDIR) $(PILBUILDDATE) \
	       $(PYBINDGENDIR) $(PYBINDGENBUILDDATE) \
	       $(WILDMAGICDIR) $(WILDMAGICBUILDDATE) \
	       $(QHULLDIR) $(QHULLBUILDDATE) \
	       $(HDF5DIR) $(HDF5BUILDDATE) \
	       $(SILODIR) $(SILOBUILDDATE) \
	       $(EIGENDIR) $(EIGENBUILDDATE) \
	       $(VOROPPDIR) $(VOROPPBUILDDATE) \
	       $(VOROPP2DDIR) $(VOROPP2DBUILDDATE) \
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
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/bin
	rm -fr $(PYTHONDIR)
	$(BZIP2) -dc $(PYTHONDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(PYTHONDIR); \
	env CONFIG_SHELL=$(CONFIG_SHELL) CC=$(PYTHONCC) CXX=$(PYTHONCXX) CFLAGS='$(PYTHONCFLAGS)' LDFLAGS='$(LDFLAGS)' HAS_HG=no ./configure --disable-ipv6 --with-gcc=$(PYTHONCC) --with-cxx-main=$(PYTHONCXX) --prefix=$(SPHERALDIR) $(PYTHONCONFFLAGS)"
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
	$(MAKE) $(NUMPYDIST)
	rm -fr $(NUMPYDIR)
	$(GZIP) -dc $(NUMPYDIST) | $(TAR) -xf -
	$(PATCH) $(NUMPYDIR)/numpy/distutils/system_info.py $(NUMPYDIR)-system_info.py.patch 
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
	$(SHELL) -ec 'cd $(PSYCOBUILDDIR); $(PYTHONEXE) setup.py install;'
	@if test -n "`$(PYTHONEXE) -c 'import psyco'`"; then \
	echo "ERROR building Psyco python extension."; \
	rm -f $(PSYCOBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PSYCODIR)
	touch $(PSYCOBUILDDATE)

#-------------------------------------------------------------------------------
# Elementtree python extension
#-------------------------------------------------------------------------------
$(ELEMENTTREEDIST):
	$(CURL) "$(ELEMENTTREEURL)" > $(ELEMENTTREEDIST)

$(ELEMENTTREEBUILDDATE):	$(PYTHONBUILDDATE) $(ELEMENTTREEDIST)
	@echo "################################################################################"
	@echo "# Building $(ELEMENTTREEDIST)"
	@echo "################################################################################"
	rm -fr $(ELEMENTTREEDIR)
	$(UNZIP) $(ELEMENTTREEDIST)
	$(SHELL) -ec 'cd $(ELEMENTTREEDIR); $(PYTHONEXE) setup.py install;'
	@if test -n "`$(PYTHONEXE) -c 'import elementtree'`"; then \
	echo "ERROR building elementtree python extension."; \
	rm -f $(ELEMENTTREEBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(ELEMENTTREEDIR)
	touch $(ELEMENTTREEBUILDDATE)

#-------------------------------------------------------------------------------
# Boost
#-------------------------------------------------------------------------------
$(BOOSTDIST):
	$(CURL) "$(BOOSTURL)" > $(BOOSTDIST)

$(BOOSTBUILDDATE):	$(BOOSTDIST)
	@echo "################################################################################"
	@echo "# Building $(BOOSTDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR) $(SPHERALDIR)/bin
	rm -fr $(BOOSTDIR)
	$(BZIP2) -dc $(BOOSTDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(BOOSTDIR); \
	PATH=$(SPHERALDIR)/bin:$$PATH; export PATH; \
	./bootstrap.sh --without-libraries=date_time,filesystem,graph,graph_parallel,iostreams,math,mpi,program_options,python,regex,serialization,signals,system,test,thread,wave --prefix=$(SPHERALDIR); \
	./bjam install;'
	@if test "$(BOOSTPYTHONTARGET)" = "python" -a ! -e $(SPHERALDIR)/lib/$(BOOSTPYTHONLIB0) ; then \
	echo "ERROR building local boost distribution. (1, no $(BOOSTPYTHONLIB0))"; \
	rm -f $(BOOSTBUILDDATE); \
	exit 1; \
	fi
	@if test "$(BOOSTPYTHONTARGET)" = "python"; then \
	rm -f $(LIBDIR)/$(BOOSTPYTHONLIB) $(SPHERALDIR)/include/boost; \
	ln -s $(SPHERALDIR)/lib/$(BOOSTPYTHONLIB0) $(LIBDIR)/$(BOOSTPYTHONLIB); \
	fi
	@if test "$(BOOSTPYTHONTARGET)" = "python" -a ! -e $(LIBDIR)/$(BOOSTPYTHONLIB); then \
	echo "ERROR building local boost distribution. (2)"; \
	rm -f $(BOOSTBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(BOOSTDIR)
	touch $(BOOSTBUILDDATE)

#-------------------------------------------------------------------------------
# gccxml
#-------------------------------------------------------------------------------
$(GCCXMLDATE):		$(CMAKEEXE) $(GCCXMLDIST)
	@echo "################################################################################"
	@echo "# Building $(GCCXMLDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(SPHERALDIR)/bin
	rm -fr $(GCCXMLDIR)
	$(BZIP2) -dc $(GCCXMLDIST) | $(TAR) -xf -
	$(MKDIR) $(GCCXMLDIR)/build
	$(SHELL) -ec 'cd $(GCCXMLDIR)/build; \
	env CC=$(GCCXMLCC) CXX=$(GCCXMLCXX) $(CMAKEEXE) -DCMAKE_INSTALL_PREFIX:PATH=$(SPHERALDIR) ..; $(MAKE); $(MAKE) install'
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
	rm -fr $(CMAKEDIR)
	$(GZIP) -dc $(CMAKEDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(CMAKEDIR); \
	env CC=$(CMAKECC) CXX=$(CMAKECXX) ./bootstrap --prefix=$(SPHERALDIR); \
	$(MAKE); \
	$(MAKE) install'
	@if ! test -e $(CMAKEEXE); then \
	echo "ERROR building cmake."; \
	exit 1; \
	fi
	rm -fr $(CMAKEDIR)

#-------------------------------------------------------------------------------
# pyMPI
#-------------------------------------------------------------------------------
#$(PYMPIDIST):
#	$(CURL) "$(PYMPIURL)" > $(PYMPIDIST)

$(PYMPIBUILDDATE):  $(PYTHONBUILDDATE) $(WILDMAGICBUILDDATE) $(PYMPIDIST)
	@echo "################################################################################"
	@echo "# Building $(PYMPIDIST)"
	@echo "################################################################################"
	rm -fr $(PYMPIDIR)
	$(BZIP2) -dc $(PYMPIDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(PYMPIDIR); \
	PATH=$(SPHERALDIR)/bin:$$PATH; export PATH; \
	env CONFIG_SHELL=$(CONFIG_SHELL) CC='$(MPICXX)' CXX='$(MPICXX)' LDFLAGS='$(LDFLAGS) $(LIBS) $(OPTFLAGS)' ./configure --prefix=$(SPHERALDIR); \
	$(MAKE) MAKEFLAGS= ; \
	$(MAKE) install"
	@if ! test -e $(PYMPIEXE); then \
	echo "ERROR building pyMPI distribution."; \
	rm -f $(PYMPIBUILDDATE); \
	exit 1; \
	fi
	$(SRCTOP)/helpers/runexecstack $(PYMPIEXE)
	rm -fr $(PYMPIDIR)
	touch $(PYMPIBUILDDATE)

#-------------------------------------------------------------------------------
# ParMETIS
#-------------------------------------------------------------------------------
$(PARMETISDIST):
	$(CURL) "$(PARMETISURL)" > $(PARMETISDIST)

$(PARMETISBUILDDATE):  $(PARMETISDIST) $(CMAKEEXE)
	@echo "################################################################################"
	@echo "# Building $(PARMETISDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(SPHERALDIR)/lib $(SPHERALDIR)/include
	rm -fr $(PARMETISDIR)
	$(GZIP) -dc $(PARMETISDIST) | $(TAR) -xf -
	env PATH="$(SPHERALDIR)/bin:$$PATH" CC="$(MPICC)" CXX="$(MPICXX)" $(MAKE) -C $(PARMETISDIR) config shared=1 prefix=$(SPHERALDIR)
	$(MAKE) MAKEFLAGS= -C $(PARMETISDIR) 
	$(MAKE) MAKEFLAGS= -C $(PARMETISDIR) install
	@if ! test -e $(SPHERALDIR)/lib/libparmetis.so -a \
	           -e $(SPHERALDIR)/include/parmetis.h; then \
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
	rm -fr $(PYGTKDIR) $(PYGTKDUMMYDATE)
	$(GZIP) -dc $(PYGTKDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(PYGTKDIR); \
	env CC=$(PYTHONCC) CXX=$(PYTHONCXX) PYTHON=$(PYTHONEXE) ./configure --enable-thread --prefix=$(SPHERALDIR); \
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
	rm -fr $(GSLDIR)
	$(GZIP) -dc $(GSLDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(GSLDIR); \
	env CC=$(CC) CXX=$(CXX) ./configure --prefix=$(SPHERALDIR) ; \
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
	rm -fr $(PYGSLDIR)
	$(GZIP) -dc $(PYGSLDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(PYGSLDIR); $(PYTHONEXE) setup.py install;'
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
	@$(MKDIR) $(SPHERALDIR)/include
	$(GZIP) -dc $(HASHMAPDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(HASHMAPDIR); \
	env cc='$(CC)' CXX='$(CXX)' ./configure --enable-namespace=google --prefix=$(SPHERALDIR); \
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
$(VISITWRITERBUILDDATE):  $(VISITWRITERDIST) $(WILDMAGICBUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(VISITWRITERDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(SPHERALDIR)/include
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
# FFTW++ interface
#-------------------------------------------------------------------------------
$(FFTWPPDIST):
	$(CURL) "$(FFTWPPURL)" > $(FFTWPPDIST)

$(FFTWPPBUILDDATE):	$(FFTWPPDIST)
	@echo "################################################################################"
	@echo "# Building $(FFTWPPDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(SPHERALDIR)/include
	rm -fr $(FFTWPPDIR)
	$(GZIP) -dc $(FFTWPPDIST) | $(TAR) -xf -
	$(CP) -f $(FFTWPPDIR)/fftw++.h $(SPHERALDIR)/include
	$(CP) -f $(FFTWPPDIR)/fftw++.cc $(SPHERALTOP)/FractalGravity
	@if ! test -e $(SPHERALDIR)/include/fftw++.h -a \
	           -e $(SPHERALTOP)/FractalGravity/fftw++.cc; then \
	echo "ERROR building FFTW++."; \
	exit 1; \
	fi
	rm -fr $(FFTWPPBUILDDATE) $(FFTWPPDIR)
	touch $(FFTWPPBUILDDATE)

#-------------------------------------------------------------------------------
# ARRAY class for use in fftw++
#-------------------------------------------------------------------------------
$(ARRAYDIST):
	$(CURL) "$(ARRAYURL)" > $(ARRAYDIST)

$(ARRAYBUILDDATE):	$(ARRAYDIST)
	@echo "################################################################################"
	@echo "# Building $(ARRAYDIST)"
	@echo "################################################################################"
	@$(MKDIR) $(SPHERALDIR)/include
	$(CP) -f $(ARRAYDIST) $(SPHERALDIR)/include
	@if ! test -e $(SPHERALDIR)/include/$(ARRAYDIST); then \
	echo "ERROR building ARRAY."; \
	exit 1; \
	fi
	touch $(ARRAYBUILDDATE)

#-------------------------------------------------------------------------------
# ZLIB
#-------------------------------------------------------------------------------
$(ZLIBDIST):
	$(CURL) "$(ZLIBURL)" > $(ZLIBDIST)

$(ZLIBBUILDDATE):	$(ZLIBDIST)
	@echo "################################################################################"
	@echo "# Building $(ZLIBDIST)"
	@echo "################################################################################"
	rm -fr $(ZLIBDIR)
	$(BZIP2) -dc $(ZLIBDIST) | $(TAR) -xf -
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
	@$(MKDIR) $(SPHERALDIR)/include
	rm -fr $(JPEGDIR)
	$(GZIP) -dc $(JPEGDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(JPEGDIR); \
	env CC=$(CC) CXX=$(CXX) CFLAGS="-fPIC" ./configure --prefix=$(SPHERALDIR) ; \
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
#$(PYBINDGENDIST):
#	$(CURL) "$(PYBINDGENURL)" > $(PYBINDGENDIST)

$(PYBINDGENBUILDDATE):	$(PYTHONBUILDDATE) $(PYBINDGENDIST)
	@echo "################################################################################"
	@echo "# Building $(PYBINDGENDIST)"
	@echo "################################################################################"
	rm -fr $(PYBINDGENDIR)
	$(UNZIP) $(PYBINDGENDIST)
	$(CP) -r $(PYBINDGENDIR)/pybindgen $(PYTHONLIBDIR)
	@if test -n "`$(PYTHONEXE) -c 'import pybindgen'`"; then \
	echo "ERROR building PYBINDGEN python extension."; \
	rm -f $(PYBINDGENBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(PYBINDGENDIR)
	touch $(PYBINDGENBUILDDATE)

# $(PYBINDGENBUILDDATE):	$(PYTHONBUILDDATE) $(PYBINDGENDIST)
# 	@echo "################################################################################"
# 	@echo "# Building $(PYBINDGENDIST)"
# 	@echo "################################################################################"
# 	rm -fr $(PYBINDGENDIR)
# 	$(BZIP2) -dc $(PYBINDGENDIST) | $(TAR) -xf -
# 	$(SHELL) -ec 'cd $(PYBINDGENDIR); ./waf configure --prefix=$(SPHERALDIR); ./waf install'
# 	@if test -n "`$(PYTHONEXE) -c 'import pybindgen'`"; then \
# 	echo "ERROR building PYBINDGEN python extension."; \
# 	rm -f $(PYBINDGENBUILDDATE); \
# 	exit 1; \
# 	fi
# 	rm -fr $(PYBINDGENDIR)
# 	touch $(PYBINDGENBUILDDATE)

#-------------------------------------------------------------------------------
# WildMagic
#-------------------------------------------------------------------------------
$(WILDMAGICDIST):
	$(CURL) "$(WILDMAGICURL)" > $(WILDMAGICDIST)

$(WILDMAGICBUILDDATE):	$(PYTHONBUILDDATE) $(WILDMAGICDIST)
	@echo "################################################################################"
	@echo "# Building $(WILDMAGICDIST)"
	@echo "################################################################################"
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

$(QHULLBUILDDATE):	$(PYTHONBUILDDATE) $(QHULLDIST) $(CMAKEEXE)
	@echo "################################################################################"
	@echo "# Building $(QHULLDIST)"
	@echo "################################################################################"
	rm -fr $(QHULLDIR) $(QHULLINCDIR)
	@$(MKDIR) -p $(SPHERALDIR)/bin $(SPHERALDIR)/man $(SPHERALDIR)/lib $(QHULLINCDIR)
	$(GZIP) -dc $(QHULLDIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(QHULLDIR)/build; \
	env CC=$(CC) CXX=$(CXX) CFLAGS='$(CFLAGS)' $(CMAKEEXE) -DCMAKE_INSTALL_PREFIX:PATH=$(SPHERALDIR) .. ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libqhullstatic.a; then \
	echo "ERROR building Qhull."; \
	rm -f $(QHULLBUILDDATE); \
	exit 1; \
	fi
	rm -fr $(QHULLDIR)
	touch $(QHULLBUILDDATE)

#-------------------------------------------------------------------------------
# mpi4py
#-------------------------------------------------------------------------------
$(MPI4PYDIST):
	$(CURL) "$(MPI4PYURL)" > $(MPI4PYDIST)

$(MPI4PYBUILDDATE):  $(PYTHONBUILDDATE) $(WILDMAGICBUILDDATE) $(MPI4PYDIST)
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
	@$(MKDIR) $(LIBDIR)
	rm -fr $(HDF5DIR)
	$(BZIP2) -dc $(HDF5DIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(HDF5DIR); \
	env CC=$(CC) CXX=$(CXX) ./configure --prefix=$(SPHERALDIR) ; \
	$(MAKE) ; \
	$(MAKE) install"
	@if ! test -e $(SPHERALDIR)/lib/libhdf5.a -a \
	           -e $(SPHERALDIR)/include/hdf5.h -a \
	           -e $(SPHERALDIR)/lib/libhdf5.$(DYLIBEXT); then \
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

$(SILOBUILDDATE):  $(SILODIST) $(HDF5BUILDDATE)
	@echo "################################################################################"
	@echo "# Building $(SILODIST)"
	@echo "################################################################################"
	@$(MKDIR) $(LIBDIR)
	rm -fr $(SILODIR)
	$(GZIP) -dc $(SILODIST) | $(TAR) -xf -
	$(SHELL) -ec "cd $(SILODIR); \
	env CC=$(CC) CXX=$(CXX) $(SILOFLAGS) ./configure --enable-fortran=no --enable-silex --enable-shared --prefix=$(SPHERALDIR) --with-hdf5='$(SPHERALDIR)/include,$(SPHERALDIR)/lib' ; \
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
	@$(MKDIR) $(SPHERALDIR)/include
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
	rm -fr $(VOROPP2DDIR) $(SPHERALDIR)/include/$(VOROPP2DINC)
	$(MKDIR) $(SPHERALDIR)/include/$(VOROPP2DINC)
	$(BZIP2) -dc $(VOROPP2DDIST) | $(TAR) -xf -
	$(SHELL) -ec 'cd $(VOROPP2DDIR)/src; $(MAKE) CC=$(CC) CXX=$(CXX) CFLAGS="-fPIC"; \
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
	@$(MKDIR) $(SPHERALDIR)/include $(SPHERALDIR)/lib
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

