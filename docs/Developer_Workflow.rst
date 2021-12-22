###################
Developer Workflow
###################

Cloning Spheral
===============

Spheral uses submodules in it's source directories in order to clone the project you must use ``--recursive``:

::

  git clone --recursive https://github.com/llnl/spheral

Setting up Spheral TPLs
=======================

Spheral uses Spack to handle Third Party Library dependencies. Spack will track dependencies between TPLs and version constraints of those TPL as the project develops. Spack is also particularly good at handling TPL setup across various configurations of compilers and build time options.

tpl-manager.py
--------------

To handle setting up a Spack instance for our Spheral builds with all of the appropriate TPLs we supply a tool called ``tpl-manager`` located at ``scripts/devtools/tpl-manager.py`` . TPL manager can be used to set up the full suite of our TPLs for all supported Spheral build configurations or it can be used to build TPL's for a single / new configuration.

``tpl-manager`` will also generate ``host-config`` files that will be needed to configure your builds to use the TPLs installed by ``tpl-manager``. See : **Using the Generated Host-Configs**

Setup Full TPL List
-------------------

The simplest way to have ``tpl-manager`` build all TPLs for all supported Spheral configurations is by running the command below from your Spheral root directory:

::

  ./scripts/devtools/tpl-manager.py --spec-list scripts/devtools/sec-list.json

This will install the Spheral Spack instance into the adjacent directory to your Spheral root dir. You can use ``--spheral-spack-dir``, if you would like to setup the spack instance somewhere else. 

**Note :** Do not set this directory within you Spheral root directory. This upsets CMake.

Single Spec Builds
------------------

If you wish to only set up TPLs for a single spec you can use ``--spec`` in your ``tpl-manager`` command:

.scripts/devtools/tpl-manager.py --spec clang@13.0.0~mpi

::

  .scripts/devtools/tpl-manager.py --spec clang@13.0.0~mpi

Above we are telling ``tpl-manager`` to build our TPLs with the clang compiler at version 13.0.0 and we are disabling ``mpi`` support with ``~mpi``.

spec-list.json
--------------


``spec-list.json`` contains a list of specs supported for common development system types. You can add or edit the specs in this file to build a variety of TPLs for various compiler combinations. The specs are grouped by the ``$SYS_TYPE`` environment variable.

Mirrors
-------

A mirror is not necessary however Spack mirrors are useful as they provide a way to cache downloaded tar files of TPLs and even store zipped binaries of previously built libraries. This is extremely useful on machines where regular development happens such as LLNL TOSS and BlueOS systems. By default ``tpl-manager`` tries to load a Spack mirror from ``/usr/WS2/davis291/SPHERAL/spheral-tpl/mirror``. As this is the predefined mirror location for LC machines. If this is not available, or you would like to use another mirror you can set ``--mirror-dir``. You can also disable mirrors entirely by setting ``--use-mirror False``.

Help
----

``tpl-manager`` supports ``-h`` or ``--help`` if you need to reference the available options.



Using the Generated Host-Configs
================================

After running ``tpl-manager`` you will see files in your Spheral root dir following the format ``<Hostname>-<sys_type>-<spec>.cmake``. For example if you run ``tpl-manager`` with ``spec-list.json`` on RZGenie you will see:

::

  rzgenie-toss_3_x86_64_ib-clang@9.0.0^mvapich2.cmake
  rzgenie-toss_3_x86_64_ib-gcc@8.1.0 build_type=Debug ^mvapich2.cmake
  rzgenie-toss_3_x86_64_ib-gcc@8.3.1^mvapich2.cmake

host-config-build.py
--------------------

A new devtool that takes the place of the previous build scripts in ``scripts/lc-builds/...`` is the ``host-config-build`` tool. Located at ``scripts/devtools/host-config-build.py``. ``host-config-build`` takes one of these generated host-config files from ``tpl-manager`` to set up you build environment with the appropriate TPLs.

--host-config
.............

You can use ``host-config-build`` to setup up your directory structure and cmake config, similarly to how the old ``scripts/lc-builds/...`` scripts worked by passing ``--host-config`` to the tool, this is a **required** argument :

::

  ./scripts/devtools/host-config-build.py --host-config "$PWD/rzgenie-toss_3_x86_64_ib-gcc@8.1.0 build_type=Debug ^mvapich2.cmake"

By default this will create two directories in the format of:

::

  build_<host-config-name>/build
  build_<host-config-name>/install

If you wish your build directory to live somewhere else, run ``host-config-build`` from that directory and use ``--source-dir`` to point at the root Spheral dir. You can set up a custom install location by passing ``--install-dir``.

Build & Install
---------------

After running ``host-config-build`` you can enter the ``build`` directory and ``make -j N install`` to build and install Spheral. From here you can build and develop manually as you usually would. 

--build
.......

If you would like the script to handle running a build and install for you for whatever reason ``--build`` exists. This will configure your CMake as usual and then launch a build and install. 

--lc-modules
............

If you use build you may need some system modules in your environment during the build and install step. you can pass these to ``host-config-build`` with ``--lc-modules`` as so:

::

  ./scripts/devtools/host-config-build.py --host-config "$PWD/rzgenie-toss_3_x86_64_ib-gcc@8.1.0 build_type=Debug ^mvapich2.cmake" --build --lc-modules "gcc/8.1.0"

If ``--build`` is not passed ``--lc-modules`` will not do anything, you will need to ensure the correct modules are in your path before building manually.

Customize CMake Options
-----------------------

With ``host-config-build`` we are still able to pass and override CMake arguments. To do this add your CMake ``-D<XXXXX>`` options to your ``host-config-build`` arguments. This is particularly useful if you want to change the ``CMAKE_BUILD_TYPE`` or use a TPL that was not installed by ``tpl-manager``.

The example below show how you would take our ``gcc@8.1.0 build_type=Debug ^mvapich2`` host-config used above, and configure with ``Release`` and a custom ``PYB11Generator`` install.

::

  ./scripts/devtools/host-config-build.py --host-config "$PWD/rzgenie-toss_3_x86_64_ib-gcc@8.1.0 build_type=Debug ^mvapich2.cmake" -DCMAKE_BUILD_TYPE=Release -Dpyb11generator_DIR=<PYB11generator_install_prefix>/lib/python2.7/site-packages/

Help
----

``host-config-build`` supports ``-h`` or ``--help`` if you need to reference the available options.
