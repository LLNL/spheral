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
