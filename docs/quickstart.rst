
QuickStart Instructions (Ubuntu 20.04)
######################################

The following commands can be used to configure and build Spheral on Ubuntu 20.04. This has been tested on a fresh install.

::

  sudo apt update
  sudp apt upgrade
  sudo apt install build-essential git gfortran mpich autotools-dev autoconf sqlite pkg-config uuid gettext cmake libncurses4-dev libgdbm-dev libffi-dev libssl-dev libexpat-dev libreadline-dev
  mkdir -p SPHERAL && cd SPHERAL
  git clone --recursive -b feature/gitlab-ci https://github.com/llnl/spheral
  cd spheral
  python3 scripts/devtools/tpl-manager.py --spec gcc@9.3.0
  python3 scripts/devtools/host-config-build.py --host-config linux-ubuntu20.04-gcc@9.3.0.cmake
  cd build_linux-ubuntu20.04-gcc@9.3.0/build
  make -j install

These commands are explained in further detail in `Obtaining, building, and installing Spheral <building.html>`_.

.. note::
  To build on a distro that is not Ubuntu 20.04 please also see `Other Distros <building.html#other-distros>`_.
