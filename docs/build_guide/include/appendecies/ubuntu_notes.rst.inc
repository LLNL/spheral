Linux Ubuntu Notes
==================

When building on any system a few basic utilities are assumed to be installed.  It's impossible to cover all the possible build environments, but a common case is a Linux Ubuntu install.  In our experience we need at least the following packages beyond the base system default, which can be easily accomplished using ``apt install``)::

    sudo apt install git autotools-dev autoconf pkg-config uuid gettext libgdbm-dev libexpat-dev cmake g++ gfortran zlib1g-dev libssl-dev libbz2-dev libreadline-dev build-essential libncurses5-dev libgdbm-dev libnss3-dev libffi-dev wget tk tk-dev libsqlite3-dev texlive-latex-recommended texlive-latex-extra dvipng

Most of these requirements are for building a full-featured Python installation.  If you also want to build the MPI parallel enabled version of Spheral you need an MPI implementation such as OpenMPI or MPICH -- OpenMPI for instance can be installed by adding the Ubuntu package ``mpich`` or ``openmpi-bin`` to the above list.

