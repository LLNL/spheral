#########################
Development Documentation
#########################


Python Module Installation Strategy
===================================

The full list of Spherals python module dependencies:

::

  python==2.7.18
  pip==9.0.1
  setuptools
  wheel
  virtualenv==20.2.2
  pyb11generator==1.0.12
  sphinx==1.8.5
  sphinx_rtd_theme==0.5.0
  docutils==0.17.1
  numpy==1.16.6
  numpy-stl==2.11.2
  matplotlib==2.2.5
  decorator==4.4.2
  h5py==2.10.0
  twine==1.15.0
  cython==0.29.21
  sobol==0.9
  scipy==1.2.3
  pipreqs==0.4.10
  importlib_metadata==2.1.1
  mpi4py==3.3.0
  gnuplot==1.8
  ats==5.2

However many of these are not required to actually build Spheral, they are purely dependencies of the Spheral python interface.

----

Build Time Packages
-------------------

These packages are **always** required when building the Spheral python interface:

::

  python==2.7.18
  pip==9.0.1
  setuptools==?
  wheel==?
  virtualenv==20.2.2
  pyb11generator==1.0.12

Optionally if we are building Spheral with ``ENABLE_DOCS=On`` we will need:

::

  sphinx==1.8.5
  sphinx_rtd_theme==0.5.0
  docutils==0.17.1

Implementation
--------------

- These packages are used during the build as they are necessary utilities to compile Spheral source code and other targets.
- Each package is designed to be install / used from a directory outside of the python ``site-packages`` directory. Therefore during the build Spheral creates a ``PYTHONENV`` list of each of the packages respective directories and calls the python executable with the ``PYTHONENV`` variable as necessary.

CMake Interface
---------------

- With the new python package management system we expose these python packages the same way we expose other TPL 	dependencies in Spheral. They can be enabled and disabled through CMake options ``<package_name>_DIR`` and ``<package_name>_BUILD`` this is useful as we move towards a Spack based TPL system where python packages are treated as separate dependencies from a given python installation.

e.g. Using a locally built / installed version of a ``pyb11generator`` package

::

  cmake ..... -Dpyb11generator_BUILD=Off -Dpyb11generator_DIR=/<local_pyb11gen_site_package_dir>

----

Runtime Packages
----------------

::

  numpy==1.16.6
  numpy-stl==2.11.2
  matplotlib==2.2.5
  decorator==4.4.2
  h5py==2.10.0
  twine==1.15.0
  cython==0.29.21
  sobol==0.9
  scipy==1.2.3
  pipreqs==0.4.10
  importlib_metadata==2.1.1

The above packages are always required by Spheral at runtime and are added to a ``requirements.txt`` file that is generated during the CMake configuration stage. If ``ENABLE_MPI=On`` then ``mpi4py`` is also added to this list of packages.

The generated ``requirements.txt`` file is used in the ``spheral-setup-env.sh`` script to install these directly to the Spheral virtual environment we create.

- Note :  To continue to support off-line installations, all package ``.tar`` files are still saved to the ``${CACHE_DIR}`` and all install commands check this directory for previously downloaded ``.tar`` files.

-----

Custom Runtime Packages
-----------------------

::

  gnuplot==1.8
  ats==5.2

These packages do not have Pypi versions of themselves and are therefore installed alongside the Runtime Packages in the Spheral virtual environment setup script.

----

Spack Integration
-----------------

As we move to a Spack based TPL system it is useful to mimic the separation of python and python-package installs. This system allows us to expose our build time python dependencies for Spack to handle. However once we are past the build stage of Spheral we do not want Spack to manage our runtime dependencies as that requires users to use Spack as an interface to Spheral: launching the spack environment, loading the spack/spheral module, dirtying the users environment etc. This approach lets us control the Spheral runtime environment by continuing to providing a ``./spheral`` command at the install directory.


====


Spheral Spack / Uberenv
=======================



Creating a build cache / mirror :
---------------------------------

Set up Local Uberenv and Spack Environment
..........................................

::

  " Go to the spheral directory we are working from and create some files adjacently
  cd <spheral_src_dir>
  mkdir -p ../uberenv-tpl/spack-env
  mkdir -p ../uberenv-tpl/mirror   " <--- only if setting up a mirror. 

  " We use uberenv to setup our local spack instance with our configs and (spheral) packages
  python3 scripts/uberenv/uberenv.py --setup-only --prefix=../uberenv-tpl

  " Initialize the spack instance, then create and enter a spack environment
  . ../uberenv-tpl/spack/share/spack/setup-env.sh
  cd ../uberenv-tpl/spack-env
  spack env create -d .
  spacktivate .

To install just the dependencies for building a spheral tpl mirror:

::

  spack install --only dependencies spheral@develop%gcc@8.1.0
  spack install --only dependencies spheral@develop%gcc@8.3.1
  ...
  spack install --only dependencies spheral@<other_specs>

We should now have all of the TPLs we need for the specs we want on this machines architecture installed in this environment `uberenv/spack-env`.

Creating a gpg Key
..................

**!!-WARNING-!!** If a Spack gpg key has previously been made read the full section before doing anything...

To create a **NEW** gpg key :

::

  spack gpg create davis291 davis291@llnl.gov
  mkdir $HOME/private_gpg_backup
  cp $SPACK_ROOT/opt/spack/gpg/*.gpg $HOME/private_gpg_backup
  cp $SPACK_ROOT/opt/spack/gpg/pubring.* <mirror_dir>

  " This will probably be needed for a /usr/gapps install, otherwise I haven't used it.
  chgrp $GROUP $MIRROR_DIR/pubring.*

Creating Our Source Mirror
..........................

::

  spack mirror create -d <mirror_dir> --all
  chmod -R g+rws <mirror_dir>

This will create a source cache mirror which stores the source tar of all of our packages installed in our current spack environment.

Creating Our Binary Mirror
..........................

First we need to add our created mirror to our current spack instance and "trust" the keys :

::

  spack mirror add spheral-tpl <mirror_dir>
  spack buildcache keys --install --trust
  mkdir -p <mirror_dir>/build_cache

Updating binary mirrors from CI / separate Spack instances
----------------------------------------------------------

- **Note**: To install to the mirrors `build_cache` from a separate Spack instance (than the one it was created with) we need to copy the original gpg keys into the current Spack instance. 

- If we are using the same Spack instance initially used to create the mirror then this is not necessary.

- This will mostly be used for updating build caches from CI

:: 
 
  mkdir -p $SPACK_ROOT/opt/spack/gpg-backup
  mv <spack_root>/opt/ spack/gpg/* ../spack/opt/spack/gpg-backup/
  cp ~/private_gpg_backup/* ../spack/opt/spack/gpg/


To add all installed packages except `Spheral` to a given mirrors `build_cache` :

::

  " This loops through all of the non external installed packages (excluding Spheral) and adds 
  " them to build cache
  for ii in $(spack find --format "yyy {name} {version} /{hash}" |
              grep -v -E "^(develop^master)" |
              grep -v spheral |
              grep "yyy" |
              cut -f4 -d" ")
  do
    spack buildcache create -k davis291 --allow-root --force -d <mirror_dir> --only=package $ii
  done

  chmod -R g+rws <mirror_dir>/build_cache



Using a Binary Mirror for Spheral dev-build:
--------------------------------------------

Follow the above instructions to **Set up local uberenv and spack environment** until :

::

  spacktivate .

Now we need to add the mirror and keys to our current Spack instance / environment :

::

  spack mirror add spheral-tpl <mirror_dir>
  spack gpg trust `find <mirror_dir> -name "*.pub"`
  " or use "spack gpg trust <mirror_dir>/build_cache/_pgp/XXXXXXXXXX.pub"

  " The below command should do what "spack gpg trust" does, but does not work as expected...
  spack buildcache keys --install --trust --force
  " I've had luck with:
  spack buildcache keys --install --trust
  " However that hasn't been reproducable 

Now we should be able to install our spheral spec using binary TPLs from the mirror :

::

  spack dev-build -d <spheral_src_dir> spheral@develop%gcc@8.1.0


Reference
---------

https://github.com/LLNL/raja-suite-dev-env/blob/8f133efcc56f81638577762863402b3a7cedb910/README.md

https://spack-tutorial.readthedocs.io/en/latest/tutorial_binary_cache.html
