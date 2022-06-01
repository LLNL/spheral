Checking/updating CMake Version
===============================

Unfortunately most recent versions of Ubuntu Linux (and derivatives such as Mint) come with an older version of CMake by default (typically something like CMake v3.10).  This is too out of date for a Spheral build, and therefore needs to be updated before configuring and building Spheral.  First, just to make sure you have this issue you should check the version of cmake that comes with your distribution::

    cmake --version

If the result is something less than version 3.18, it's worth updating before starting to configure Spheral.  How to accomplish this varies by platform, but for the common case of Ubuntu (and similar ``apt`` based distributions) something like the following should suffice.

1. First, remove any existing cmake installation using apt::

     sudo apt remove --purge cmake

2. Follow the directions on the `Kitware site at this link <https://apt.kitware.com/>`_ to add their repository for installing packages.

3. Install a current version of cmake with::

     sudo apt install cmake

Check the final version again to make sure you have what you expect::

     cmake --version

