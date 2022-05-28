.. include:: <isoamsa.txt>

Quickstart
##########

Create a directory structure and clone spheral.
::

  mkdir -p SPHERAL && cd SPHERAL
  git clone --recursive https://github.com/llnl/spheral
  cd spheral

Build our TPL dependencies from source with the Spheral tpl-management tool (``tpl-manager.py``).
::

  ./scripts/devtools/tpl-manager.py --spec gcc@8.3.1

.. warning::
  This command needs to be run under an allocation as any TPLs that need to be built will be built in parallel.

.. note::
  This command will generate a ``.cmake`` file with the naming convention ``<system-type>-<compiler-spec>``. The following commands will refer to this format as ``<host-config>`` for generalization across operating systems and architectures. You will need to substitute the correct format in the following commands. 

Set up a build/install directory structure and configure cmake.
::

  ./scripts/devtools/host-config-build.py --host-config <host-config>.cmake

Build and install Spheral.
::

  cd build_<host-config>/build
  make -j N install

.. warning::
 ``make`` should be run under an allocation as it will take considerable resources on LC to build and install Spheral.

Run a basic smoke test for Spheral
::

  cd ../install
  ./spheral -c "import Spheral"

Run our full test suite.
::

  ./spheral-atstest test/integration.ats

These commands are explained in further sections.
