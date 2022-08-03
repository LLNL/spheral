

..
   BASIC SMOKE TEST
   ----------------------------------------

[smoke_test-start]
Basic Smoke Test
================

After a build and install it's recommended to perform a quick smoke test with Spheral to see if the Spheral environment was installed and all of the libraries were built and linked together correctly.

From your install directory run:
::

  ./spheral -c "import Spheral"

This will run the Spheral Python process and try importing the Spheral modules.  If successful you should see a banner printed with information about the version of Spheral, after which Python exits.

[smoke_test-end]



..
   RUNNING ATS
   ----------------------------------------

[ex_ats_test-start]
ATS Testing Suite
=================

Spheral uses ATS to execute a suite of parallel tests. To run this on an external system we need to use Spheral's virtual-env installation of ATS, as external users will not have access to some LC available scripts.

From the install directory run:
::

  ./.venv/bin/ats -e spheral tests/integration.ats 

[ex_ats_test-end]
[lc_ats_test-start]
ATS Testing Suite
=================

Spheral uses ATS (the Python Automated Testing System) to execute a suite of parallel tests. During install a script is generated ``spheral-atstest`` which handles calling the ats script in ``/usr/gapps/Spheral/bin`` and passing the generated ``spheral`` executor script to it.

From the install directory run:
::

  ./spheral-atstest tests/integration.ats 

[lc_ats_test-end]



..
   ATS FILTER OPTIONS
   ----------------------------------------

[filter_options-start]
ATS Filters and Options
-----------------------

Appropriate filters are injected into the spheral-atstest script based on build configurations. However if you are using ats directly or wish to override some of these filters you can pass them as defined below.

Non MPI Filter
..............

If Spheral was built without ``MPI`` support we will need to pass a filter to our ``ats`` command to prevent parallel tests from being performed:
:: 

  --filter='"np<2"'
  
This stops the ATS suite from attempting to run any tests that rely on more than one process/rank.

Debug Build Filter
..................

If Spheral was built in Debug mode it is recommended to pass the below filter if you value your time.
::

  --filter='"level<100"'

These filters stack when invoked. So if you are running the test suite on a non-mpi debug build the command would be:

[filter_options-end]

[ex_filter_options_cmd-start]
::

  ./.venv/bin/ats -e spheral tests/integration.ats --filter='"np<2"' --filter='"level<100"'

[ex_filter_options_cmd-end]
[lc_filter_options_cmd-start]
::

  ./spheral-atstest tests/integration.ats --filter='"np<2"' --filter='"level<100"'

[lc_filter_options_cmd-end]

[lc_cuda_filter-start]
CUDA Filter
...........

When building Spheral with CUDA support and without MPI support on Blue-OS systems we need to pass --smpi_off to ATS. When using spheral-atstest this flag is injected for you. However, for veboseness the flag is:
::

  --addOp "--smpi_off"

.. note::
   We need to pass ``--addOp`` as a way to forward the ``-smpi_off`` option to ATS. This is because spheral-atstest uses an intermediate layer to control ATS execution.

If you are using ATS directly you only need to pass:
::

  --smpi_off

[lc_cuda_filter-end]
