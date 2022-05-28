.. include:: <isoamsa.txt>

Quickstart
##########

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
