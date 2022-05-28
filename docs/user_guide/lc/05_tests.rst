Running Tests
#############

Basic Smoke Test
================

After a build and install it's recommended to perform a quick smoke test with Spheral to see if the Spheral environment was installed and all of the libraries were built and linked together correctly.

From your install directory run:
::

  ./spheral -c "import Spheral"


ATS Testing Suite
=================

Spheral uses ATS to execute a suite of parallel tests. During install a script is generated ``spheral-atstest`` which handles calling the ats script in ``/usr/gapps/Spheral/bin`` and passing the generated ``spheral`` executor script to it.

From the install directory run:
::

  ./spheral-atstest tests/integration.ats 

ATS Filters and Options
-----------------------

We will need to pass some filters to ``ATS`` depending on how we built Spheral.

Non MPI Filter
..............

If Spheral was built without ``MPI`` support we will need to pass the filter to our ``ats`` command. This stops the ATS suite from attempting to run any tests that rely on more than one process/rank.
:: 

  --filter='"np<2"'
  
Debug Build Filter
..................

If Spheral was built in Debug mode it is recommended to pass the below filter if you value your time.
::

  --filter='"level<100"'

These filters stack when invoked. So if you are running the test suite on a non-mpi debug build the command would be:

::

  ./spheral-atstest tests/integration.ats --filter='"np<2"' --filter='"level<100"'


