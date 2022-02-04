.. include:: <isoamsa.txt>

###############################################
Obtaining, building, and installing Spheral
###############################################

Cloning Spheral
===============

Spheral minimally requires a C++11 compliant compiler.  If you use git to clone the Spheral source note Spheral includes git submodules: `BLT <https://github.com/LLNL/blt>`_ and `Uberenv <https://github.com/LLNL/uberenv>`_.  In order to ensure such submodules are properly downloaded when cloning Spheral be sure to use the ``--recurse-submodules`` or ``--recursive`` git options:

::

  git clone --recursive https://github.com/jmikeowen/Spheral

If you forget to use the ``--recursive`` argument or if you checkout from a different branch you must run:

::

  git submodule update --init --recursive

A Two Stage Build System
========================

The Spheral build system works in two stages. 
 - First: Building and setting up Third Party Libraries (TPL)s.
 - Second: Building and installing Spheral.

Setting up Spheral TPLs
=======================

Spheral uses `Spack <https://github.com/llnl/spack>`_ under the hood to handle Third Party Library dependencies. Spack will track dependencies between TPLs and version constraints of those TPLs as the project develops. Spack is also particularly good at handling TPL setup across various configurations of compilers and build time options.

tpl-manager.py
--------------

To handle setting up a Spack instance for Spheral with the appropriate TPLs we supply a tool called ``tpl-manager`` located at ``scripts/devtools/tpl-manager.py`` . TPL manager can be used in two ways; it can be used to build TPL's for a single compiler configuration (ideal for users) or to setup the full range of TPLs for all of our supported compiler configurations (this is ideal for developers).

``tpl-manager`` will generate `CMake host-config <https://llnl-blt.readthedocs.io/en/develop/tutorial/host_configs.html>`_ files that can be used to configure your builds to use the TPLs installed by ``tpl-manager``.

.. note::
   You do not need to use ``tpl-manager`` to setup TPLs for Spheral. TPLs can be built individually and passed to the Spheral CMake system or built through your own spack installation. See `Custom TPL Installation`_ and `Custom Spack Installation`_ for more details.

Single Spec Builds
------------------

Installing TPLs for a single spec you can use ``--spec`` in your ``tpl-manager`` command:

::

  ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1

This will install the Spheral Spack instance into the adjacent directory to your Spheral root dir. You can use ``--spheral-spack-dir`` if you would like to setup the spack instance somewhere else. 

Above we are telling ``tpl-manager`` to build our TPLs with gcc at version 8.3.1. By default this will build with ``mpi`` support, however we can disable ``mpi`` support by appending ``~mpi`` to our spec.

.. note::
   For more information on ``spec`` syntax please see the spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.

Developer tpl-manager.py features
---------------------------------

.. note::
   Much of this discussion centers around using ``tpl-manager.py`` on LC (Livermore Computing) machines, as this is where much Spheral development occurs.  ``spack`` and ``tpl-manager.py`` are not tied to this environment however.  It is simple enough to use ``tpl-manager.py`` in your own environment as needed, however you must be aware of which compilers and libraries (MPI, OpenMP, etc.) you wish to utilize.

Setup Full TPL List
...................

The simplest way to have ``tpl-manager`` build all TPLs for all supported Spheral configurations is by running this command from your Spheral root directory:

::

  ./scripts/devtools/tpl-manager.py --spec-list scripts/devtools/spec-list.json

spec-list.json
..............

``spec-list.json`` contains a list of specs supported for common development system types. You can add or edit the specs in this file to build a variety of TPLs for various compiler combinations. The specs are grouped by the ``$SYS_TYPE`` environment variable for LC (Livemore Computing) machines and will default to x86_64 for eveything else.

If you are not building on an LC system you may want to create your own ``.json`` file with defaults appropriate for your system.

Mirrors
.......

A mirror is not necessary, however Spack mirrors are useful as they provide a way to cache downloaded tar files of TPLs and even store zipped binaries of previously built libraries. This is extremely useful on machines where regular development happens such as LLNL TOSS and BlueOS systems. By default ``tpl-manager`` tries to load a Spack mirror from ``/usr/gapps/Spheral/spheral-spack-tpls/mirror``, as this is the predefined mirror location for LC machines. If this is not available, or you would like to use another mirror, you can set ``--mirror-dir``. You can also disable mirrors entirely by setting ``--no-mirror``.

Help
----

``tpl-manager`` supports ``-h`` or ``--help`` if you need to reference the available options.


Basic Spheral Build (host-config-build.py)
==========================================

After running ``tpl-manager`` you will see files in your Spheral root dir following the format ``<sys_type>-<spec>.cmake``. For example if you run ``tpl-manager`` with ``spec-list.json`` on a toss_3_x86_64_ib system you will see:

::

  toss_3_x86_64_ib-clang@9.0.0^mvapich2.cmake
  toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake
  toss_3_x86_64_ib-gcc@8.3.1^mvapich2.cmake

However if you ran ``tpl-manager`` with only a single ``--spec`` e.g. ``gcc@8.1.0^mvapcich2``, you will only see:

::

  toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake

The default build will try to build Spheral with Python, MPI, and OpenMP support. This will be enough to run Spheral from the locally built python and pass the test suite.

.. note::
  A basic build & install from this point would look as follows:

  ::
    
    ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake
    cd build_toss_3_x86_64_ib-gcc@8.1.0^mvapich2/build
    make -j <N> install
    cd ../install/
    ./spheral -c "import Spheral"

  The following sections detail these commands further.

host-config-build.py
--------------------

A new devtool that takes the place of the previous build scripts in ``scripts/lc-builds/...`` is the ``host-config-build`` tool. Located at ``scripts/devtools/host-config-build.py``. ``host-config-build`` takes one of these generated host-config files from ``tpl-manager`` to setup your CMake build with the appropriate TPLs. ``host-config-build`` by default also sets up a basic build/install directory structure. 

::

  ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake"

``--host-config`` is a **required** argument of the tool, by default this will create two directories in the format of:

::

  build_<host-config-name>/build
  build_<host-config-name>/install

If you wish your build directory to live somewhere else, run ``host-config-build`` from that directory and use ``--source-dir`` to point at the root Spheral dir. You can set up a custom install location by passing ``--install-dir`` if you do not like creating build/install trees inside your source dir.

.. note::
   ``host-config-build.py`` is simply a wrapper around CMake, and you can pass ordinary CMake commands through this script as well. It is also possible to directly run CMake rather than ``host-config-build.py``; ``host-config-build.py`` reports the CMake line it is using, so this can be a good starting point if you need to run CMake manually yourself.  See :ref:`cmake_customization` for more details.

Build & Install
---------------

::

  cd build_toss_3_x86_64_ib-gcc@8.1.0^mvapich2/build
  make -j <N> install


After running ``host-config-build`` you can enter the ``build`` directory and ``make -j N install`` to build and install Spheral. You can build and develop as you would normally from this directory. Alternatively the ``host-config-build.py`` tools also provides arguments to automate the build/install process of Spheral if you desire. 

--build
.......

If you would like the script to handle running a build and install for you ``--build`` exists. This will configure your CMake as usual and then launch a build and install. 

--lc-modules
............

If you use build you may need some system modules in your environment during the build and install step. you can pass these to ``host-config-build`` with ``--lc-modules`` as so:

::

  ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake --build --lc-modules "gcc/8.1.0"

If ``--build`` is not passed ``--lc-modules`` will not do anything, you will need to ensure the correct modules are in your path before building manually.

.. _cmake_customization:

Customize CMake Options
-----------------------

With ``host-config-build`` we are still able to pass and override CMake arguments (See: `Spheral / CMake Configurations`_). To do this add your CMake ``-D<XXXXX>`` options to your ``host-config-build`` arguments. This is particularly useful if you want to change the ``CMAKE_BUILD_TYPE`` or use a TPL that was not installed by ``tpl-manager``.

The example below show how you would take our ``gcc@8.1.0^mvapich2`` host-config used above, and configure with ``Release`` and a custom ``PYB11Generator`` install.

::

  ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake" -DCMAKE_BUILD_TYPE=Release -Dpyb11generator_DIR=<PYB11generator_install_prefix>/lib/python2.7/site-packages/

Help
----

``host-config-build`` supports ``-h`` or ``--help`` if you need to reference the available options.


Basic Spheral Build (Manual)
============================

``host-config-build.py`` is a tool for convenience if you are comfortable with using CMake and wish to setup your own build/install directory structure that is still very easy to do.

::

  mkdir -p Spheral_release/BUILD && cd Spheral_release/BUILD
  cmake -C ../../Spheral/toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake \
        -DCMAKE_INSTALL_PREFIX=`cd ..; pwd` ../../Spheral
  make -j<N> install
  ../spheral -c "import Spheral"

In this example we performed our build in the directory ``Spheral_release/BUILD``, and installed all binaries and libraries in ``Spheral_release``. The final line is simply a test that the installed Python client can load the Spheral Python modules.

The CMake command ``-C ../../Spheral/toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake`` is how we tell the Spheral CMake system to use the TPLs we installed with ``tpl-manager.py`` for ``gcc v8.1.0``.

The somewhat obscure command ``-DCMAKE_INSTALL_PREFIX=`chdir ..; pwd``` just specifies the install directory as the full path to ``Spheral_release``. Alternatively you can specify this path explicitly, such as ``-DCMAKE_INSTALL_PREFIX=/usr/local/Spheral_release``, if that were the correct path.

.. note::
   Although Spheral is simply a set of Python modules, it installs in a Python virtual environment, so the script ``spheral`` installed at the top level of the install tree is designed to load the virtual environment on invocation, and then unload it on completion.

.. note::
   Many users of CMake like to place the build directory as a subdirectory of the cloned code, so many examples you'll see online use "``cmake ..``".  All that matters really is that the final path on the CMake command line point to the top of the source tree.

Spheral / CMake Configurations
==============================

C++ Only Build
--------------

If you do not need to build the Python interface you can build the compiled C++ libraries alone with ``-DENABLE_CXXONLY=On``.  This will skip the Python wrapping stage of the build. 

By default Spheral builds the libraries as shared objects.  If instead you would like to build the C++ libraries as static libs use ``-DENABLE_STATIC_CXXONLY=On``.

Third party libraries and Spheral
---------------------------------

Upon first install third party libraries (TPL) tar files and source will be installed through an external network. TPLs are cached within the Spheral build directory tree for future builds off network. To completely turn off installation of TPL's use ``-DBUILD_TPLS=Off``.

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

In this section we list the CMake variables that can be tweaked for a Spheral build.  Where appropriate the options are listed, with the default value in *italics*.

``CMAKE_BUILD_TYPE``   (Debug, *Release*, RelWithDebInfo, MinSizeRel)
  Choose the type of build -- for more information see the `CMake documentation <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`_.

``CMAKE_INSTALL_PREFIX``
  The top-level path for installing Spheral include files, libraries, and any Python modules or documentation.

``ENABLE_CXXONLY`` (On, *Off*)
  Do not build python wrappers for Spheral.

``ENABLE_STATIC_CXXONLY`` (On, *Off*)
  Do not build python wrappers for Spheral. Build static library files for Spheral.

``<TPL-Name-Here>_DIR``
  Directory of previously built TPL.

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

``ENABLE_OPENSUBDIV`` (*On*, Off)
  Install the Opensubdiv library along with the Spheral interface to it.  Opensubdiv is a `Pixar provided library <https://github.com/PixarAnimationStudios/OpenSubdiv>`_, which Spheral uses to implement refinement of polyhedra for some specialized problem generation capabilities.

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

Appendecies
===========

Linux Ubuntu Notes
------------------

When building on any system a few basic utilities are assumed to be installed.  It's impossible to cover all the possible build environments, but a common case is a Linux Ubuntu install.  In our experience we need at least the following packages beyond the base system default, which can be easily accomplished using ``apt install``)::

    sudo apt install git autotools-dev autoconf pkg-config uuid gettext libgdbm-dev libexpat-dev cmake g++ gfortran zlib1g-dev libssl-dev libbz2-dev libreadline-dev build-essential libncurses5-dev libgdbm-dev libnss3-dev libffi-dev wget tk tk-dev libsqlite3-dev texlive-latex-recommended texlive-latex-extra dvipng

Most of these requirements are for building a full-featured Python installation.  If you also want to build the MPI parallel enabled version of Spheral you need an MPI implementation such as OpenMPI or MPICH -- OpenMPI for instance can be installed by adding the Ubuntu package ``mpich`` or ``openmpi-bin`` to the above list.

Checking/updating CMake version
...............................

Unfortunately most recent versions of Ubuntu Linux (and derivatives such as Mint) come with an older version of CMake by default (typically something like CMake v3.10).  This is too out of date for a Spheral build, and therefore needs to be updated before configuring and building Spheral.  First, just to make sure you have this issue you should check the version of cmake that comes with your distribution::

    cmake --version

If the result is something less than version 3.18, it's worth updating before starting to configure Spheral.  How to accomplish this varies by platform, but for the common case of Ubuntu (and similar ``apt`` based distributions) something like the following should suffice.

1. First, remove any existing cmake installation using apt::

     sudo apt remove --purge cmake

2. Follow the directions on the `Kitware site at this link <https://apt.kitware.com/>`_ to add their repository for installing packages.

3. Install a current version of cmake with::

     sudo apt install cmake

Check the final version again to make sure you have what you expect::

     cmake --version

WSL2 Notes
-----------

The Windows Subsystem for Linux (WSL) is a useful method of development on Windows 10 based systems.  If going this route we recommend having at least WSL2 for best results -- the original version of WSL (WSL1) also functioned, but is `significantly` slower for jobs such as compilation.

For the most part using an Ubuntu based WSL environment works just using the Ubuntu notes above.  However, one aspect of WSL2 needs to be adjusted.  The build process requires a fair amount of memory (in particular for a few of the Python binding modules), so we recommend having at least 32GB of swap space available.  On WSL2 this is accomplished by creating a ``.wslconfig`` file in your Windows home directory containing at least the following lines::

    [wsl2]
    swap=32GB

Custom Spack Installation
-------------------------

Building Spheral TPLs with your own Spack installation will require deeper knowledge of `how Spack works <https://spack.readthedocs.io/en/latest/>`_. All of the steps to set up Spheral with your own spack installation are not detailed here, however you will want to at least:

 - Point your spack instances `repo <https://spack.readthedocs.io/en/latest/repositories.html?highlight=repo#spack-repo>`_ at the ``scripts/spack/packages/`` dir. This contains all of our changes to spack packages that have not yet made it to upstream Spack.
 - You will want to model your ``compiler.yaml`` and ``packages.yaml`` files off of those found in ``scripts/spack/configs/`` (`Spack Configuration Files <https://spack.readthedocs.io/en/latest/configuration.html#configuration>`_).
   
Further notes on setting up Spack and how it is used with the Spheral dev-tools scripts can be found in `Development Documentation: Spheral Spack / Uberenv <Development_Documentation.html#spheral-spack-uberenv>`_.

Other Distros
-------------

To build Spheral TPLs on an x86_64 Linux OS that isn`t Ubuntu 20.04 you will need to edit the file ``scripts/spack/configs/x86_64/compilers.yaml``. For example an Arch based Manjaro distro`s ``compilers.yaml`` looks like this.

::

  compilers:
  - compiler:
      spec: gcc@11.1.0
      paths:
        cc: /usr/bin/gcc
        cxx: /usr/bin/g++
        f77: /usr/bin/gfortran
        fc: /usr/bin/gfortran
      flags: {}
      operating_system: manjaro21
      target: x86_64
      modules: []
      environment: {}
      extra_rpaths: []

Most notably you will need to edit the ``spec`` of your compiler, ``paths``, and ``operating_system``.

Please also ensure that the packages listed in `Linux Ubuntu Notes`_ are installed on your system through your respective package manager.
