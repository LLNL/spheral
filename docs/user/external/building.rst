.. include:: <isoamsa.txt>

###################
External User Guide
###################

Quickstart for External Users
#############################

This page details commands for setting up and running Spheral on an Ubuntu 20.04 system. 

Update and install necessary package dependencies.
::

  sudo apt update
  sudo apt upgrade
  sudo apt install build-essential git gfortran mpich autotools-dev autoconf sqlite pkg-config uuid gettext cmake libncurses4-dev libgdbm-dev libffi-dev libssl-dev libexpat-dev libreadline-dev

.. warning::
   For alternative Linux distros your mileage may vary, ensure you are installing compatible packages to the ones listed above.

Create a directory structure and clone spheral.
::

  mkdir -p SPHERAL && cd SPHERAL
  git clone --recursive https://github.com/llnl/spheral
  cd spheral

Build our TPL dependencies from source with the Spheral tpl-management tool (``tpl-manager.py``).
::

  python3 scripts/devtools/tpl-manager.py --spec gcc

.. note::
   This command sequence assumes ``gcc`` is installed and will use the version in your system path. If you wish to use a different compiler (such as ``clang``) ensure it is in your path and replace ``gcc`` with the compiler name of your choice (e.g. ``clang``).

.. note::
  This command will generate a ``.cmake`` file with the naming convention ``<system-type>-<compiler-spec>``. The following commands will refer to this format as ``<host-config>`` for generalization across operating systems and architectures. You will need to substitute the correct format in the following commands. 

Set up a build/install directory structure and configure cmake.
::

  python3 scripts/devtools/host-config-build.py --host-config <host-config>.cmake

Build and install Spheral.
::

  cd build_<host-config>/build
  make -j install

Run a basic smoke test for Spheral
::

  cd ../install
  ./spheral -c "import Spheral"

Run our full test suite.
::

  ./.venv/bin/ats -e spheral test/integration.ats

These commands are explained in further detail below.

Building, and Installing Spheral for External Users
###################################################

Cloning Spheral
===============

If you use git to clone the Spheral source be aware Spheral includes git submodules: `BLT <https://github.com/LLNL/blt>`_ and `Uberenv <https://github.com/LLNL/uberenv>`_.  In order to ensure such submodules are properly downloaded when cloning Spheral be sure to use the ``--recursive`` git option:

::

  git clone --recursive https://github.com/jmikeowen/Spheral

If you forget to use the ``--recursive`` argument or if you checkout from a different branch you must run:

::

  git submodule update --init --recursive

A Two Stage Build System
========================

The Spheral build system works in two stages. 
 - Stage 1: Building and setting up Third Party Libraries (TPL)s.
 - Stage 2: Building and installing Spheral.

 .. note::
   Stage 1 is technically optional however highly recommended, especially for first time users. If you do not wish to skip Stage 1 you can build the TPLs manually and pass them to CMake when configuring Spheral. (See ...)

Setting up Spheral Third Party Libraries (TPLs)
===============================================

Spheral uses `Spack <https://github.com/llnl/spack>`_ under the hood to handle Third Party Library dependencies. Spack will track dependencies between TPLs and version constraints of TPLs as the project develops. Spack is also particularly good at handling TPL setup across various configurations, compilers and build time options.

tpl-manager.py
--------------

Spheral provides an internal tool ``tpl-manager.py`` to attempt to deobfuscate the spack process from the user. ``tpl-manager`` s primary responsibilities are:

 - Setup a local ``Spack`` instance for Spheral.
 - Generate a dependency tree of third party libraries (relative to the provided configuration).
 - Build and install all dependent libraries in local ``Spack`` instance.
 - Generate a `CMake host-config <https://llnl-blt.readthedocs.io/en/develop/tutorial/host_configs.html>`_  file for configuring Spheral builds.

``tpl-manager`` is located at ``scripts/devtools/tpl-manager.py``. ``tpl-manager`` can be used in two ways:

 - Build TPL's for a single compiler configuration (ideal for users).
 - Build and setup the full range of TPLs for all of our supported compiler configurations (ideal for developers).

 .. note::
    You do not need to use ``tpl-manager`` to setup TPLs for Spheral. TPLs can be built individually and passed to the Spheral CMake system or built through your own spack installation. See `Customize CMake Options`_ for more details.

Running tpl-manager
-------------------

``tpl-manager`` takes a ``--spec`` argument to determine what compiler to use and what configuration we want to build Spheral in.

::

  ./scripts/devtools/tpl-manager.py --spec gcc

This will install the local Spheral Spack instance into the adjacent directory of your Spheral root dir. You can use ``--spheral-spack-dir`` if you would like to setup the spack instance somewhere else. 

Above we are telling ``tpl-manager`` to build our TPLs with the ``gcc`` that is in our path. By default this will build with ``+mpi`` support, however we can disable ``mpi`` support for the TPLs and Spheral by appending ``~mpi`` to our spec.
::

  ./scripts/devtools/tpl-manager.py --spec gcc~mpi

.. note::
   By default we have ``python`` bindings enabled (``+python``) and docs disabled (``~docs``). Therefore the spec ``gcc+mpi+python~docs`` will build the same TPL set as ``gcc``.For more information on ``spec`` syntax please see the spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.

.. note::
   Spheral minimally requires a C++11 compliant compiler.

Help
----

``tpl-manager`` supports ``-h`` or ``--help`` if you need to reference the available options.


Configuring a Spheral Build (host-config-build.py)
==================================================

After running ``tpl-manager`` you will see a file in your Spheral root directory following the format ``<sys_type>-<spec>.cmake``. 

For example if you run ``tpl-manager`` with only a single ``--spec`` e.g. ``gcc~mpi``, you will only see:

::

  <sys_type>-gcc~mpi.cmake

.. note::
  A basic build & install from this point would look as follows:

  ::
    
    ./scripts/devtools/host-config-build.py --host-config <sys_type>-gcc.cmake
    cd build_<sys_type>-gcc/build
    make -j <N> install
    cd ../install/
    ./spheral -c "import Spheral"

  The following sections detail these commands further.

host-config-build.py
--------------------

The ``host-config-build`` tool is located at ``scripts/devtools/host-config-build.py``. ``host-config-build`` takes a host-config file and sets up Spheral's CMake with the appropriate TPLs. ``host-config-build`` by default also sets up a basic build/install directory structure. 

::

  ./scripts/devtools/host-config-build.py --host-config <sys_type>-gcc.cmake"

``--host-config`` is a **required** argument of the tool, by default this will create two directories in the format of:

::

  build_<host-config>/build
  build_<host-config>/install

.. note::
   ``host-config-build.py`` is simply a wrapper around CMake. It is possible to directly run CMake rather than ``host-config-build.py``; ``host-config-build.py`` reports the CMake line it is using, so this can be a good starting point if you need to run CMake manually yourself.  See :ref:`manual_build` for more details.


Changing the build and install locations with ``host-config-build``
-------------------------------------------------------------------

Custom Build Directory
......................

If you wish your build directory to live somewhere else, run ``host-config-build`` from that directory and use ``--source-dir`` to point at the root Spheral dir.

Custom Install Directory
........................

You can setup a custom install location by using ``--install-dir``.

--build
.......

If you would like the script to handle running a build and install for you ``--build`` exists. This will configure your CMake as usual and then launch a build and install. 

Customize CMake Options
-----------------------

With ``host-config-build`` we are still able to pass and override CMake arguments (See: `Spheral / CMake Configurations <cmake_config.html>`_). To do this add your CMake ``-D<XXXXX>`` options to your ``host-config-build`` arguments. This is particularly useful if you want to change the ``CMAKE_BUILD_TYPE`` or use a TPL that was not installed by ``tpl-manager``.

The example below show how you would take our ``gcc`` host-config used above, and configure with ``Release`` and a custom ``PYB11Generator`` install.

::

  ./scripts/devtools/host-config-build.py --host-config <sys_type>-gcc.cmake" -DCMAKE_BUILD_TYPE=Release -Dpyb11generator_DIR=<PYB11generator_install_prefix>/lib/python2.7/site-packages/

Help
----

``host-config-build`` supports ``-h`` or ``--help`` if you need to reference the available options.

Build & Install
===============

::

  cd build_<sys_type>-gcc/build
  make -j <N> install


After running ``host-config-build`` you can enter the ``build`` directory and ``make -j N install`` to build and install Spheral. You can build and develop as you would normally from this directory. Alternatively the ``host-config-build.py`` tools also provides arguments to automate the build/install process of Spheral if you desire. 


.. _manual_build:

Manual Spheral Build
--------------------

``host-config-build.py`` is a tool for convenience if you are comfortable with using CMake and wish to setup your own build/install directory structure that is still very easy to do.

::

  mkdir -p Spheral_release/BUILD && cd Spheral_release/BUILD
  cmake -C ../../Spheral/<sys_type>-gcc.cmake \
        -DCMAKE_INSTALL_PREFIX=`cd ..; pwd` ../../Spheral
  make -j<N> install
  ../spheral -c "import Spheral"

In this example we performed our build in the directory ``Spheral_release/BUILD``, and installed all binaries and libraries in ``Spheral_release``. The final line is simply a test that the installed Python client can load the Spheral Python modules.

The CMake command ``-C ../../Spheral/<sys_type>-gcc.cmake`` is how we tell the Spheral CMake system to use the TPLs we installed with ``tpl-manager.py`` for ``gcc``.

The somewhat obscure command ``-DCMAKE_INSTALL_PREFIX=`chdir ..; pwd``` just specifies the install directory as the full path to ``Spheral_release``. Alternatively you can specify this path explicitly, such as ``-DCMAKE_INSTALL_PREFIX=/usr/local/Spheral_release``.

.. note::
   Although Spheral is simply a set of Python modules, it installs in a Python virtual environment, so the script ``spheral`` installed at the top level of the install tree is designed to load the virtual environment on invocation, and then unload it on completion.

.. note::
   Many users of CMake like to place the build directory as a subdirectory of the cloned code, so many examples you'll see online use "``cmake ..``".  All that matters really is that the final path on the CMake command line point to the top of the source tree.

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

Running Tests
#############

Basic Smoke Test
================

After a build and install it's recommended to perform a quick smoke test with Spheral to see if the Spheral environment was installed and all of the libraries were built and linked together correctly.

From your install directory run:
::

  ./spheral -c "import Spheral"


ATS Testing Suite
=================

Spheral uses ATS to execute a suite of parallel tests. To run this on an external system we need to use Spheral's virtual-env installation of ATS, as external users will not have access to some LC available scripts.

From the install directory run:
::

  ./.venv/bin/ats -e spheral tests/integration.ats 

ATS Filters and Options
-----------------------

We will need to pass some filters to ``ATS`` depending on how we built Spheral.

Non MPI Filter
..............

If Spheral was built without ``MPI`` support we will need to pass the filter to our ``ats`` command. This stops the ATS suite from attempting to run any tests that rely on more than one process/rank.
:: 

  --filter='"np<2"'
  
Debug Build Filter
..................

If Spheral was built in Debug mode it is recommended to pass the below filter if you value your time.
::

  --filter='"level<100"'

These filters stack when calling them. So if you are running``ats`` on a non-mpi debug buid your tests command would be as such:

::

  ./.venv/bin/ats -e spheral tests/integration.ats --filter='"np<2"' --filter='"level<100"'
