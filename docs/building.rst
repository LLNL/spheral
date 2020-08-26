.. include:: <isoamsa.txt>

###############################################
Obtaining, building, and installing Spheral
###############################################

Spheral minimally requires a C++11 compliant compiler.  If you use git to clone the Spheral source note Spheral includes a git submodule: `BLT <https://github.com/LLNL/blt>`_.  In order to ensure such submodules are properly downloaded when cloning Spheral be sure to use the ``--recurse-submodules`` or ``--recursive`` git options:

::

  git clone --recursive https://github.com/jmikeowen/Spheral

If you forget to use the ``--recursive`` argument or if you checkout from a different branch you must run:

::

  git submodule update --init --recursive

Basic Build
-----------

The default build will try to build Spheral with Python and MPI support. This will be enough to run Spheral from the locally built python and pass the test suite.

I like to keep my build & install files separate from the git cloned source, so I usually build in a completely separate directory as follows:

::

  git clone --recursive https://github.com/jmikeowen/Spheral
  mkdir -p Spheral_release/BUILD && cd Spheral_release/BUILD
  cmake -DCMAKE_INSTALL_PREFIX=`chdir ..; pwd` ../../spheral
  make -j<N> install
  ../python/bin/python2.7 -c "import Spheral"

In this example we performed our build in the directory ``Spheral_release/BUILD``, and installed all binaries and libraries in ``Spheral_release``.  Note this includes the TPL libraries, which are downloaded under ``Spheral_release/BUILD`` and installed to ``Spheral_release``.  The final line is simply a test that the installed Python client can load the Spheral Python modules.

The somewhat obtuse command ``-DCMAKE_INSTALL_PREFIX=`chdir ..; pwd``` just specifies the install directory as the full path to ``Spheral_release``.  Alternatively you can specify this path explicitly, such as ``-DCMAKE_INSTALL_PREFIX=/usr/local/Spheral_release``, if that were the correct path.

Note many users of CMake like to place the build directory as a subdirectory of the cloned code, so many examples you'll see online use "``cmake ..``".  All that matters really is that the final path on the CMake command line point to the top of the source tree.

C++ Only Build
--------------

If you do not need to build the Python interface you can build the compiled C++ libraries alone with ``-DENABLE_CXXONLY=On``.  This will skip the Python wrapping stage of the build. 

By default Spheral builds the libraries as shared objects.  If instead you would like to build the C++ libraries as static libs use ``-DENABLE_STATIC_CXXONLY=On``.

Third party libraries and Spheral
---------------------------------

Upon first install third party libraries (TPL) tar files and source will be installed through an external network. TPLs are cached within the Spheral build directory tree for future builds off network. To completely turn off installation of TPL's use ``-DINSTALL_TPLS=Off``.

For just the C++ compiled Spheral a number of TPLs are required (and automatically installed by default):

- Boost
- Python
- Eigen
- Pybind11
- Polytope
- HDF5
- Silo
- Qhull
- ANEOS
- Opensubdiv

There are also a number of Python packages that are automatically installed during the third party build, including:

- Scipy
- Sobol
- Cython
- Twine
- sphinx_rtd_theme
- h5py
- decorator
- Matplotlib
- mpi4py
- Numpy-stl

By default the Spheral TPLs are installed to the path specified on the CMake command line with ``-DCMAKE_INSTALL_PREFIX=<full path>``.

To pass a custom path to a TPL build use the ``-D<TPL-Name-Here>_DIR=<full path>`` option. CMake will assume the directory passed has the same structure as if it were built the library itself.  For example, to specify a preinstalled version of Boost:

::

  cmake -DBOOST_DIR=$HOME/my_boost_build_dir ...

To force a download and build from a specific TAR for a TPL use ``-D<TPL-Name-Here>_URL=<url path>``. This can be either a local file path or URL:

::

   cmake -DBOOST_URL=https://my.tarfiles.com/boost/boost-1.2.3.4.tar ...

OpenMP/MPI
----------

OpenMP and MPI support is handled through BLT.  Use the option flags ``-DENABLE_OPENMP`` and ``-DENABLE_MPI`` respectively, choosing ``ON`` or ``OFF`` as appropriate.  

CMake variables
--------------------

In this section we list the CMake variables that can be tweaked for a Spheral build.  Where appropriate the options are listed, with the default value in *italics*.

``CMAKE_BUILD_TYPE``   (Debug, *Release*, RelWithDebInfo, MinSizeRel)
  Choose the type of build -- for more information see the `Cmake documentation <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`_.

``CMAKE_INSTALL_PREFIX``
  The top-level path for installing Spheral include files, libraries, and any Python modules or documentation.  This is synonymous with and replaces the older ``SPHERAL_INSTALL_DIR``.

``ENABLE_CXXONLY`` (On, *Off*)
  Do not build python wrappers for Spheral.

``ENABLE_STATIC_CXXONLY`` (On, *Off*)
  Do not build python wrappers for Spheral. Build static library files for Spheral.

``BUILD_TPLS`` (*On*, Off)
  Option to install TPLs or not during configuration stage.

``<TPL-Name-Here>_DIR``
  Directory of previously built TPL.

``<TPL-Name-Here>_URL`` 
  URL or local path to zip/tar file of TPL to download and install.
  Defaults to CMake defined URL/cache,  see cmake/InstallLibraries.cmake

``<TPL-Name-Here>_BUILD`` (*On*, Off)
  Tell the build system to build or not to build the given TPL. A ``<TPL>_DIR`` must be provided, otherwise it will search a default directory.

``ENABLE_OPENMP`` (*On*, Off)
  Support for OpenMP.

``ENABLE_MPI`` (*On*, Off)
  Support for MPI.

``ENABLE_1D`` (*On*, Off)
  Build Spheral with 1D support.

``ENABLE_2D`` (*On*, Off)
  Build Spheral with 2D support.

``ENABLE_3D`` (*On*, Off)
  Build Spheral with 3D support.

``ENABLE_ANEOS`` (*On*, Off)
  Install the ANEOS (Analytics Equation of State) package along with the Spheral interface to it.  This is a legacy equation of state frequently used for geophysical materials.  See descriptions in the `iSALE <https://github.com/isale-code/M-ANEOS>`_ documentation.

``ENABLE_HELMHOLTZ`` (*On*, Off)
  Compile the included Helmholtz equation of state, typically used in astrophysical calculations. See a discussion `here <http://cococubed.asu.edu/code_pages/eos.shtml>`_.

``ENABLE_TIMER`` (*On*, Off)
  Enable timer information from Spheral.

``DBC_MODE`` (None, All, Pre)
  Set the compile time design by contract (DBC) mode for Spheral.  Design by contract statements are very useful developer tools, whereby the developer can insert tests in the code as they write it.  These statements are both useful for tracking down bugs with fine-grained testing throughout the code, as well as useful documentation in the code about what sort of conditions are expected to hold.

  +------+---------------------------------------------------------------------------------+
  | None | Design by contract not enforced                                                 |
  +------+---------------------------------------------------------------------------------+
  | All  | All design by contract (``REQUIRE``, ``ENSURE``, ``CHECK``) statements active   |
  +------+---------------------------------------------------------------------------------+
  | Pre  | Only prerequisites (``REQUIRE``) statements active                              |
  +------+---------------------------------------------------------------------------------+

  Note the default depends on the ``CMAKE_BUILD_TYPE``:

  ``CMAKE_BUILD_TYPE=Debug`` |xrArr| default ``DBC_MODE`` is ``All``

  In all other cases the default is ``None``.

  It is worth noting ``DBC_MODE=All`` is quite expensive at run time (of order 4x more), so this is not intended to be active for a release/production compilation of Spheral.

``ENABLE_WARNINGS`` (On, *Off*)
  Enable compiler warnings.

``ENABLE_BOUNDCHECKING`` (On, *Off*)
  If building with the Gnu compilers enable STL bound checking by passing -D_GLIBCXX_DEBUG=1 to the compiler. 
  Note, this is a very expensive option at runtime!

``ENABLE_DOCS`` (On, *Off*)
  Choose whether or not to build this documentation.

``SPHINX_EXECUTABLE``
  Specify where the Sphinx executable is that should be used to build documentation.  If not given, assumes the Spheral built Sphinx will be used.

``SPHINX_THEME`` (*sphinx_rtd_theme*)
  Give the Sphinx theme to use when generating documentation.  Default based on read the docs theme.

``SPHINX_THEME_DIR``
  Where to look for Sphinx themes.
