Updating Ubuntu
###############

This guide assumes the use of an Ubuntu 20.04 system and using ``apt`` as the package manager. For other distrobutions please install the corresponding packages. 

.. note::
  Future steps (especially those detailed in :ref:`ex_tpl`) are assuming packages are installed under ``/usr/bin``, ``/usr/lib`` etc.

::

  sudo apt update
  sudo apt upgrade
  sudo apt install build-essential git gfortran mpich autotools-dev autoconf sqlite pkg-config uuid gettext cmake libncurses-dev libgdbm-dev libffi-dev libssl-dev libexpat-dev libreadline-dev liblapack-dev libbz2-dev locales python python3 unzip libtool wget curl tk-dev

