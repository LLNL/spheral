CMake Configurations
====================

This section outlines the common CMake options Spheral provides for configuring how to build Spheral.

C++ Only Build
--------------

If you do not need to build the Python interface you can build the compiled C++ libraries alone with ``-DENABLE_CXXONLY=On``.  This will skip the Python wrapping stage of the build. 

By default Spheral builds the libraries as shared objects.  If instead you would like to build the C++ libraries as static libs use ``-DENABLE_STATIC_CXXONLY=On``.

Third party libraries and Spheral
---------------------------------

For just the C++ compiled Spheral a number of TPLs are required:

- Zlib
- Boost
- Python
- Eigen
- Polytope
- HDF5
- Silo
- Qhull
- M-ANEOS
- Opensubdiv
- Polyclipper
- Conduit
- Axom
- Adiak
- Caliper

There are also a number of libraries / python packages that are required for compiling the python bindings and executing Spheral at runtime:

- Python
- pip
- setuptools
- pybind11
- pyb11generator
- sphinx
- sphinx_rtd_theme
- Scipy
- Sobol
- Cython
- Twine
- h5py
- decorator
- Matplotlib
- mpi4py
- Numpy-stl

Custom TPL Installation
.......................

You can build the Spheral TPLs manually or even with your own spack installation to bypass the use of ``tpl-manager.py``. Custom built TPL installations can be passed to Spheral's CMake with ``-D<tpl-name>_DIR=<tpl-install-prefix>``.

::

  cmake -DBOOST_DIR=$HOME/my_boost_build_dir ...


OpenMP/MPI
----------

OpenMP and MPI support is handled through BLT.  Use the option flags ``-DENABLE_OPENMP`` and ``-DENABLE_MPI`` respectively, choosing ``ON`` or ``OFF`` as appropriate.  

CMake variables
--------------------

In this section we list the CMake variables that can be tweaked for a Spheral build.  Where appropriate the options are listed, with the default value in **bold**.

``CMAKE_BUILD_TYPE``   (Debug, **Release**, RelWithDebInfo, MinSizeRel)
  Choose the type of build -- for more information see the `CMake documentation <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`_.

``CMAKE_INSTALL_PREFIX``
  The top-level path for installing Spheral include files, libraries, and any Python modules or documentation.

``ENABLE_CXXONLY`` (On, **Off**)
  Do not build python wrappers for Spheral.

``ENABLE_STATIC_CXXONLY`` (On, **Off**)
  Do not build python wrappers for Spheral. Build and export static library files for Spheral.

``ENABLE_SHARED`` (**On**, Off)
  Build Spheral C++ libraries as shared libraries.

``ENABLE_DEV_BUILD`` (On, **Off**)
  Builds separate internal C++ libraries for faster code development.

``SPHERAL_ENABLE_VVI`` (On, *Off*)
  Builds Spheral with the experimental Value-View Interface for GPU porting capabilities.

``<TPL-Name-Here>_DIR``
  Directory of previously built TPL.

``ENABLE_STATIC_TPL`` (On, **Off**)
  Link static libraries instead of shared for HDF5 and Conduit.

``ENABLE_OPENMP`` (**On**, Off)
  Support for OpenMP.

``ENABLE_MPI`` (**On**, Off)
  Support for MPI.

``ENABLE_1D`` (**On**, Off)
  Build Spheral with 1D support.

``ENABLE_2D`` (**On**, Off)
  Build Spheral with 2D support.

``ENABLE_3D`` (**On**, Off)
  Build Spheral with 3D support.

``ENABLE_ANEOS`` (**On**, Off)
  Install the ANEOS (Analytics Equation of State) package along with the Spheral interface to it.  This is a legacy equation of state frequently used for geophysical materials.  See descriptions in the `iSALE <https://github.com/isale-code/M-ANEOS>`_ documentation.

``ENABLE_HELMHOLTZ`` (**On**, Off)
  Compile the included Helmholtz equation of state, typically used in astrophysical calculations. See a discussion `here <http://cococubed.asu.edu/code_pages/eos.shtml>`_.

``ENABLE_OPENSUBDIV`` (**On**, Off)
  Install the Opensubdiv library along with the Spheral interface to it.  Opensubdiv is a `Pixar provided library <https://github.com/PixarAnimationStudios/OpenSubdiv>`_, which Spheral uses to implement refinement of polyhedra for some specialized problem generation capabilities.

``ENABLE_TIMER`` (On, **Off**)
  Enable Caliper timer information for Spheral.

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

  ``CMAKE_BUILD_TYPE=Debug`` default ``DBC_MODE`` is ``All``

  In all other cases the default is ``None``.

  It is worth noting ``DBC_MODE=All`` is quite expensive at run time (of order 4x more), so this is not intended to be active for a release/production compilation of Spheral.

``ENABLE_WARNINGS`` (On, **Off**)
  Enable compiler warnings.

``ENABLE_BOUNDCHECKING`` (On, **Off**)
  If building with the Gnu compilers enable STL bound checking by passing -D_GLIBCXX_DEBUG=1 to the compiler. 
  Note, this is a very expensive option at runtime!

``ENABLE_NAN_EXCEPTIONS`` (On, **Off**)
  Raise exceptions in the C++ code when floating-point exceptions occur.  Gnu compilers only.

``ENABLE_DOCS`` (On, **Off**)
  Choose whether or not to build this documentation.

``SPHERAL_NETWORK_CONNECTED`` (**On**, Off)
  Spheral assumes there exists a network connection. Disable this to force pip to build python environments using only ``SPHERAL_PIP_CACHE_DIR``.

``SPHERAL_PIP_CACHE_DIR`` (**~/.cache/spheral_pip**)
  Default location Spheral will search for cached pip packages.

``SPHINX_EXECUTABLE``
  Specify where the Sphinx executable is that should be used to build documentation.  If not given, assumes the Spheral built Sphinx will be used.

``SPHINX_THEME`` (**sphinx_rtd_theme**)
  Give the Sphinx theme to use when generating documentation.  Default based on read the docs theme.

``SPHINX_THEME_DIR``
  Where to look for Sphinx themes.


