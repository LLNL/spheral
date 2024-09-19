Diagnostics
###########

Spheral uses Caliper to preform diagnostics, such as timing. To enable this functionality in the code, Spheral needs to be configured with ``ENABLE_TIMER=ON``. Otherwise the diagnostic regions are no-ops for improved preformance.
::

  ./scripts/devtools/host-config-build.py <sys_type>-<spec>.cmake -DENABLE_TIMER=ON


Querying using Caliper
======================

Caliper is configured and started through the :kbd:`cali::ConfigManager`.
The :kbd:`cali::ConfigManager` is wrapped in a :kbd:`TimerMgr` singleton class, which has a python interface.
:kbd:`TimerMgr` is initialized and started in the :kbd:`InitTimers` routine which is called in :kbd:`commandLine()` in ``src/SimulationControl/SpheralOptionParser.py``.
By default, the Caliper configuration is set to ``spot`` and outputs Caliper files (``.cali``).
For the default configuration, the Caliper files are named based on what file is being run, for example:
::

   python Noh-cylindrical-2d.py

will produce timing files called
::

   Noh-cylindrical-2d_YEAR_MONTH_DATE_TIME.cali

where the file name includes the current date and time.
The Caliper file name can be specified using the command line
::

   python Noh-cylindrical-2d.py --caliperFilename 'new_test_name.cali'

Non-default Caliper configurations can be set at the command line using ``--caliperConfig`` like so
::

   python Noh-cylindrical-2d.py --caliperConfig 'runtime-report(output=time.txt),calc.inclusive,region.count'

.. note::
   The above configuration produces timing results similar to the previous :kbd:`Spheral::Timer` method. This results in a file named ``time.txt`` with cumulative times for the nested regions as well as a count of how many times each region ran.

Similarly, a non-default Caliper configuration can be read in from a JSON file using ``--caliperConfigJSON`` and providing the file name.
Lastly, Caliper timers can be turned off using ``--caliperConfig none``.

There are many different Caliper configurations to view various information. Here are some extra links for those who want to read or experiment with other features in Caliper that can be incorporated into Spheral:

  * `Configuration basics <https://software.llnl.gov/Caliper/CaliperBasics.html#more-on-configurations>`_
  * `Builtin Configuration <https://software.llnl.gov/Caliper/BuiltinConfigurations.html>`_
  * `Manual Configuration <https://software.llnl.gov/Caliper/configuration.html>`_
  * `Output Format <https://software.llnl.gov/Caliper/OutputFormats.html>`_


Adding Region Timers in C++
===========================

So far there are two different types of regions in Spheral, using the following macros:
::

  TIME_FUNCTION
  TIME_BEGIN("timer_name")
  TIME_END("timer_name")

- :kbd:`TIME_FUNCTION` can be added to the very beginning of a function and creates a region for the entire function using the function's name. :kbd:`TIME_FUNCTION` uses just the function name and no class or parameter information, so be careful when using this method with functions that could share names.

- :kbd:`TIME_BEGIN("timer_name")` and :kbd:`TIME_END("timer_name")` create a region between the two different calls and use the string (in this case timer_name) as the name.


Adding Region Timers in Python
==============================

Region timers can be added inside the python code using the following function calls:
::

   TimerMgr.timer_start("some_function")
   some_function_call()
   TimerMgr.timer_end("some_function")

.. note::
   IMPORTANT: All timers must have both a start and end call. Otherwise, memory issues will occur.
