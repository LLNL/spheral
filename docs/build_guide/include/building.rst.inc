..
   BUILD AND INSTALL
   ----------------------------------------

[build-start]
Build & Install
###############

::

  cd build_<sys_type>-gcc/build
  make -j <N> install

Where ``<N>`` should be replaced with the number of parallel make processes you want to launch (typically the number of cores on your CPU).

After configuring Spheral, enter the ``build`` directory and run ``make -j <N> install``. Alternatively the ``host-config-build`` tool provides options to automate the build and install process. 
[build-end]

[lc_alloc_warning-start]
.. warning::
   ``make -j <N>`` will launch a number of compilation jobs in parallel, so on shared machines such as the Livermore Computing (LC) environment it is highly recommended to launch this make process on a specially allocated node.

[lc_alloc_warning-end]

[install-start]
Although Spheral is simply a set of Python modules, during the install stage Spheral will also set up a Python virtual environment. The install will also download and install all required runtime python libraries needed to execute the full Spheral test suite.

In the install directory users will find a ``spheral`` script. This can be run from the install directory as ``./spheral`` and will drop users into the avtive virtual environment with acces to all of the Spheral libraries and dependencies.
[install-end]
