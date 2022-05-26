.. include:: <isoamsa.txt>

#############
LC User Guide
#############

Quickstart for LC Users
#######################

This page details commands for setting up and running Spheral on LC systems.

Create a directory structure and clone spheral.
::

  mkdir -p SPHERAL && cd SPHERAL
  git clone --recursive https://github.com/llnl/spheral
  cd spheral

Build our TPL dependencies from source with the Spheral tpl-management tool (``tpl-manager.py``).
::

  ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1

.. warning::
  This command needs to be run under an allocation as any TPLs that need to be built will be built in parallel.

.. note::
  This command will generate a ``.cmake`` file with the naming convention ``<system-type>-<compiler-spec>``. The following commands will refer to this format as ``<host-config>`` for generalization across operating systems and architectures. You will need to substitute the correct format in the following commands. 

Set up a build/install directory structure and configure cmake.
::

  ./scripts/devtools/host-config-build.py --host-config <host-config>.cmake

Build and install Spheral.
::

  cd build_<host-config>/build
  make -j N install

.. warning::
 ``make`` should be run under an allocation as it will take considerable resources on LC to build and install Spheral.

Run a basic smoke test for Spheral
::

  cd ../install
  ./spheral -c "import Spheral"

Run our full test suite.
::

  ./spheral-atstest test/integration.ats

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

.. warning::
   If running on LC be sure to launch ``tpl-manager`` in a resource allocation, as ``tpl-manager`` will take advantage of parallel builds when compiling libraries.

Running tpl-manager
-------------------

``tpl-manager`` takes a ``--spec`` argument to determine what compiler to use and what configuration we want to build Spheral in.

::

  ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1

This will install the local Spheral Spack instance into the adjacent directory of your Spheral root dir. You can use ``--spheral-spack-dir`` if you would like to setup the spack instance somewhere else. 

Above we are telling ``tpl-manager`` to build our TPLs with ``gcc`` at version ``8.3.1``. By default this will build with ``+mpi`` support, however we can disable ``mpi`` support for the TPLs and Spheral by appending ``~mpi`` to our spec.
::

  ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1~mpi

.. note::
   By default we have ``python`` bindings enabled (``+python``) and docs disabled (``~docs``). Therefore the spec ``gcc@8.3.1+mpi+python~docs`` will build the same TPL set as ``gcc@8.3.1``.For more information on ``spec`` syntax please see the spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.

.. note::
   ``gcc@8.3.1`` is the recommended compiler for building Spheral on LC (Livermore Computing). If you are not building Spheral on an LC machine you can pass just ``--spec gcc`` and it will try to pick up the ``gcc`` in your path. Spheral minimally requires a C++11 compliant compiler.


Developer tpl-manager.py features
---------------------------------

Much of this discussion centers around using ``tpl-manager.py`` on LC (Livermore Computing) machines, as this is where most Spheral development occurs.  ``spack`` and ``tpl-manager.py`` are not tied to this environment however.  It is simple enough to use ``tpl-manager.py`` in your own environment as needed, however you must be aware of which compilers and libraries (MPI, OpenMP, etc.) you wish to utilize.

Setup Full TPL List
...................

The simplest way to have ``tpl-manager`` build all TPLs for all supported Spheral configurations is by running this command from your Spheral root directory:

::

  ./scripts/devtools/tpl-manager.py --spec-list scripts/devtools/spec-list.json

spec-list.json
..............

``spec-list.json`` contains a list of specs supported for common development system types. You can add or edit the specs in this file to build a variety of TPLs for various compiler combinations. The specs are grouped by the ``$SYS_TYPE`` environment variable for LC (Livemore Computing) machines and will default to x86_64 for eveything else.

If you are not building on an LC system you may want to create your own ``.json`` file with defaults appropriate for your system.

Help
----

``tpl-manager`` supports ``-h`` or ``--help`` if you need to reference the available options.


Configuring a Spheral Build (host-config-build.py)
==================================================

After running ``tpl-manager`` you will see a file in your Spheral root directory following the format ``<sys_type>-<spec>.cmake``. 

For example if you run ``tpl-manager`` with ``spec-list.json`` on a ``toss_3_x86_64_ib`` system you will see:

::

  toss_3_x86_64_ib-clang@9.0.0.cmake
  toss_3_x86_64_ib-gcc@8.1.0.cmake
  toss_3_x86_64_ib-gcc@8.3.1.cmake

However if you ran ``tpl-manager`` with only a single ``--spec`` e.g. ``gcc@8.1.0~mpi``, you will only see:

::

  toss_3_x86_64_ib-gcc@8.1.0~mpi.cmake

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

The ``host-config-build`` tool is located at ``scripts/devtools/host-config-build.py``. ``host-config-build`` takes a host-config file and sets up Spheral's CMake with the appropriate TPLs. ``host-config-build`` by default also sets up a basic build/install directory structure. 

::

  ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake"

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

.. warning::
   If running on LC and using ``--build`` be sure to launch in a resource allocation, as ``--build`` will take advantage of parallel compilation.

--lc-modules
............

If you use build you may need some system modules in your environment during the build and install step. you can pass these to ``host-config-build`` with ``--lc-modules`` as so:

::

  ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake --build --lc-modules "gcc/8.1.0"

If ``--build`` is not passed ``--lc-modules`` will not do anything, you will need to ensure the correct modules are in your path before building manually.


Customize CMake Options
-----------------------

With ``host-config-build`` we are still able to pass and override CMake arguments (See: `Spheral / CMake Configurations <cmake_config.html>`_). To do this add your CMake ``-D<XXXXX>`` options to your ``host-config-build`` arguments. This is particularly useful if you want to change the ``CMAKE_BUILD_TYPE`` or use a TPL that was not installed by ``tpl-manager``.

The example below show how you would take our ``gcc@8.1.0^mvapich2`` host-config used above, and configure with ``Release`` and a custom ``PYB11Generator`` install.

::

  ./scripts/devtools/host-config-build.py --host-config toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake" -DCMAKE_BUILD_TYPE=Release -Dpyb11generator_DIR=<PYB11generator_install_prefix>/lib/python2.7/site-packages/

Help
----

``host-config-build`` supports ``-h`` or ``--help`` if you need to reference the available options.

Build & Install
===============

::

  cd build_toss_3_x86_64_ib-gcc@8.1.0^mvapich2/build
  make -j <N> install


After running ``host-config-build`` you can enter the ``build`` directory and ``make -j N install`` to build and install Spheral. You can build and develop as you would normally from this directory. Alternatively the ``host-config-build.py`` tools also provides arguments to automate the build/install process of Spheral if you desire. 


.. _manual_build:

Manual Spheral Build
--------------------

``host-config-build.py`` is a tool for convenience if you are comfortable with using CMake and wish to setup your own build/install directory structure that is still very easy to do.

::

  mkdir -p Spheral_release/BUILD && cd Spheral_release/BUILD
  cmake -C ../../Spheral/toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake \
        -DCMAKE_INSTALL_PREFIX=`cd ..; pwd` ../../Spheral
  make -j<N> install
  ../spheral -c "import Spheral"

In this example we performed our build in the directory ``Spheral_release/BUILD``, and installed all binaries and libraries in ``Spheral_release``. The final line is simply a test that the installed Python client can load the Spheral Python modules.

The CMake command ``-C ../../Spheral/toss_3_x86_64_ib-gcc@8.1.0^mvapich2.cmake`` is how we tell the Spheral CMake system to use the TPLs we installed with ``tpl-manager.py`` for ``gcc v8.1.0``.

The somewhat obscure command ``-DCMAKE_INSTALL_PREFIX=`chdir ..; pwd``` just specifies the install directory as the full path to ``Spheral_release``. Alternatively you can specify this path explicitly, such as ``-DCMAKE_INSTALL_PREFIX=/usr/local/Spheral_release``.

.. note::
   Although Spheral is simply a set of Python modules, it installs in a Python virtual environment, so the script ``spheral`` installed at the top level of the install tree is designed to load the virtual environment on invocation, and then unload it on completion.

.. note::
   Many users of CMake like to place the build directory as a subdirectory of the cloned code, so many examples you'll see online use "``cmake ..``".  All that matters really is that the final path on the CMake command line point to the top of the source tree.

Appendecies
===========

Custom Spack Installation
-------------------------

Building Spheral TPLs with your own Spack installation will require deeper knowledge of `how Spack works <https://spack.readthedocs.io/en/latest/>`_. All of the steps to set up Spheral with your own spack installation are not detailed here, however you will want to at least:

 - Point your spack instances `repo <https://spack.readthedocs.io/en/latest/repositories.html?highlight=repo#spack-repo>`_ at the ``scripts/spack/packages/`` dir. This contains all of our changes to spack packages that have not yet made it to upstream Spack.
 - You will want to model your ``compiler.yaml`` and ``packages.yaml`` files off of those found in ``scripts/spack/configs/`` (`Spack Configuration Files <https://spack.readthedocs.io/en/latest/configuration.html#configuration>`_).
   
Further notes on setting up Spack and how it is used with the Spheral dev-tools scripts can be found in `Development Documentation: Spheral Spack / Uberenv <Development_Documentation.html#spheral-spack-uberenv>`_.
