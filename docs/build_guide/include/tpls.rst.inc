[intro-start]
The Spheral build system works in two stages. 

 - Stage 1: Building and setting up Third Party Libraries (TPL)s.
 - Stage 2: Building and installing Spheral.

Spheral uses `Spack <https://github.com/llnl/spack>`_ under the hood to handle Third Party Library dependencies. Spack will track dependencies between TPLs and version constraints of TPLs as the project develops. Spack is also particularly good at handling TPL setup across various configurations, compilers and build time options.
[intro-end]


[tpl_manager-start]
tpl-manager.py
==============

Spheral provides a tool (``tpl-manager.py``) in an attempt to simplify the spack process for the user. ``tpl-manager`` 's primary responsibilities are to:

 - Set up a local ``Spack`` instance for Spheral.
 - Generate a dependency tree of third party libraries (relative to the provided configuration).
 - Build and install all dependent libraries in the local ``Spack`` instance.
 - Generate a `CMake host-config <https://llnl-blt.readthedocs.io/en/develop/tutorial/host_configs.html>`_  file for configuring Spheral builds.

``tpl-manager`` is located at ``scripts/devtools/tpl-manager.py``. ``tpl-manager`` can be used in two ways:

 - Build TPL's for a single compiler configuration (ideal for users).
 - Build and setup the full range of TPLs for all of our supported compiler configurations (ideal for developers).

[tpl_manager-end]

[ex_tpl_manager_notes-start]

 .. note:
    The External User Guide does not cover building the full range of tpl configurations. Please see the LC User Guide for more informtaion.

[ex_tpl_manager_notes-end]
[lc_tpl_manager_notes-start]

 .. note:
    The External User Guide does not cover building the full range of tpl configurations. Please see the LC User Guide for more informtaion.

[lc_tpl_manager_notes-end]



[ex_running_tpl_manager-start]
Running tpl-manager
===================

``tpl-manager`` requires ``python3`` and should be run from the root Spheral directory.

``tpl-manager`` takes a ``--spec`` argument to determine what compiler to use and what configuration we want to build Spheral in.

::

  python3 scripts/devtools/tpl-manager.py --spec gcc

This will install the local Spheral Spack instance into the adjacent directory of your Spheral root dir. You can use ``--spheral-spack-dir`` if you would like to setup the spack instance somewhere else. 

Above we are telling ``tpl-manager`` to build our TPLs with the ``gcc`` that is in our path. By default this will build with ``+mpi`` support, however we can disable ``mpi`` support for the TPLs and Spheral by appending ``~mpi`` to our spec.
::

  python3 scripts/devtools/tpl-manager.py --spec gcc~mpi

.. note::
   By default we have ``python`` bindings enabled (``+python``) and docs disabled (``~docs``). Therefore the spec ``gcc+mpi+python~docs`` will build the same TPL set as just ``gcc``.For more information on ``spec`` syntax please see the spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.

.. note::
   Spheral minimally requires a C++14 compliant compiler.

ERROR: invalid spack config dir
===============================

If you are trying to run ``tpl-manager.py`` on an operating system other than Ubuntu20.04, you will
see an error to the effect of:
``[ERROR: invalid spack config dir: /<path>/scripts/spack/configs/<OperatingSystem><Version> ]``

We define configuration files for Ubuntu20.04, as well as our common LLNL operating systems.
You will need to create a set of files for your own system.

The configuration files tell spack where the packages installed in :ref:`Required System Packages`
are located and what version they are. We have provided a ``generic`` set of config files to 
help in setting this up for you.

#. Copy the directory ``scripts/spack/configs/generic`` to ``scripts/spack/config/<OperatingSystem><Version>`` (you want to match the name of the directory to the one ``tpl-manager.py`` expects to find).


#. For each package within the ``packages.py`` file of your new folder edit the version number to be the same as
   what is installed on your system. There are a number of ways to retrieve versions for a given package:

   * Most version numbers should be searchable through your package manager, there are however, some system libraries that may not be managed by your package manager.

   * If the package has an executable, often you can run with some form of ``-V`` or ``--version``. e.g. for mpich:``mpiexec --version`` will report the MPI version.

   * For packages that only provide libraries, often the system library will be symlinked to one with the version as the extension. e.g. ``ls -lha /usr/lib/x86_64_gnu-linux/libreadlines.so`` will show it is symlinked to ``libreadlines.so.8.1``.


#. Each package requires the ``prefix:`` of the installation be provided. In most cases ``/usr`` is sufficient. Typically packages installed with a package manager will place files in: ``/usr/bin``, ``/usr/share``, ``/usr/lib``, ``/usr/lib64``, sometimes ``/usr/lib/x86_64.../``. Here the common prefix is ``/usr``. 

   * If you are building Spheral on a system where you don't have permissions to run package manager and install to ``/usr``, then you might have installed the system packages somewhere else. In that case, replace the ``prefix:`` path for those given packages as necessary.


[ex_running_tpl_manager-end]

[lc_running_tpl_manager-start]
Running tpl-manager
===================

``tpl-manager`` requires ``python3`` and should be run from the root Spheral directory.

``tpl-manager`` takes a ``--spec`` argument to determine what compiler to use and what configuration we want to build Spheral in.

::

  ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1

This will install the local Spheral Spack instance into the adjacent directory of your Spheral root dir. You can use ``--spheral-spack-dir`` if you would like to setup the spack instance somewhere else. 

.. warning::
   If running on LC be sure to launch ``tpl-manager`` in a resource allocation, as ``tpl-manager`` will take advantage of parallel builds when compiling libraries.  An example appropriate for use on rzgenie would be::

    srun --exclusive ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1

Above we are telling ``tpl-manager`` to build our TPLs with ``gcc`` at version ``8.3.1``. By default this will build with ``+mpi`` support, however we can disable ``mpi`` support for the TPLs and Spheral by appending ``~mpi`` to our spec.
::

  ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1~mpi

.. note::
   By default we have ``python`` bindings enabled (``+python``) and docs disabled (``~docs``). Therefore the spec ``gcc@8.3.1+mpi+python~docs`` will build the same TPL set as just ``gcc@8.3.1``.For more information on ``spec`` syntax please see the spack documentation on `specs-dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_.

.. note::
   ``gcc@8.3.1`` is the recommended compiler for building Spheral on LC (Livermore Computing). Spheral minimally requires a C++14 compliant compiler.

CUDA Support
------------

To build Spheral TPLs and configure Spheral to build with CUDA on LC systems we need to add the spack variant ``+cuda`` to our spec. We also need to define the cuda hardware architecture we are targetting Spheral to run on. This can be done with ``cuda_arch=XX``. On most LC NVIDIA systems we will be targetting ``sm_70`` architecture. Therefore a gcc spec with cuda support will look as such.::

  ./scripts/devtools/tpl-manager.py --spec "gcc@8.3.1+cuda cuda_arch=70"

.. note::
   We need to use "" here as spack requires spacing to define ``cuda_arch``. This ensures tpl-manager doesn't register the spec as separate arguments.

Setup Full TPL List
-------------------

The simplest way to have ``tpl-manager`` build all TPLs for all supported Spheral configurations is by running this command from your Spheral root directory:

::

  ./scripts/devtools/tpl-manager.py --spec-list scripts/devtools/spec-list.json

spec-list.json
--------------

``spec-list.json`` contains a list of specs supported for common development system types. You can add or edit the specs in this file to build a variety of TPLs for various compiler combinations. The specs are grouped by the ``$SYS_TYPE`` environment variable for LC (Livemore Computing) machines and will default to x86_64 for eveything else.

[lc_running_tpl_manager-end]

[help-start]
Help
----

``tpl-manager`` supports ``-h`` or ``--help`` if you need to reference the available options.
[help-end]

[manual_tpl-start]
Manually installing TPLs
========================

Although it is HIGHLY recommended, you do not need to use ``tpl-manager`` to setup TPLs for Spheral. TPLs can be built manually or built through your own spack installation and passed to the Spheral CMake system. See :ref:`manual_config` for more details.
[manual_tpl-end]
