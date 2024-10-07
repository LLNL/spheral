Code Performance Diagnostics
############################

Spheral uses Caliper to preform code diagnostics, such as timing. To enable this functionality in the code, Spheral needs to be configured with ``ENABLE_TIMER=ON``. Otherwise, the timing regions are no-ops for improved preformance.
::

  ./scripts/devtools/host-config-build.py <sys_type>-<spec>.cmake -DENABLE_TIMER=ON


Using Caliper
=============

Caliper is configured and started through the ``cali::ConfigManager``.
The ``cali::ConfigManager`` is wrapped in a ``TimerMgr`` singleton class, which has a python interface.

.. note::
   ``TimerMgr`` is initialized in ``src/SimulationControl/SpheralTimingParser.py`` which is called during ``commandLine()`` in ``src/SimulationControl/SpheralOptionParser.py``. This is because ``commandLine()`` is almost always invoked directly near the start of a problem. However, if ``commandLine()`` is not called, the timer manager would need to be configured and started directly using the ``TimerMgr`` class. See :ref:`below <manual_caliper>` for more details.

By default, the Caliper configuration is set to ``spot`` and outputs Caliper files (``.cali``).

There are many different Caliper configurations to view various information. Here are some extra links for those who want to read or experiment with other features in Caliper that can be incorporated into Spheral:

  * `Configuration basics <https://software.llnl.gov/Caliper/CaliperBasics.html#more-on-configurations>`_
  * `Builtin Configuration <https://software.llnl.gov/Caliper/BuiltinConfigurations.html>`_
  * `Manual Configuration <https://software.llnl.gov/Caliper/configuration.html>`_
  * `Output Format <https://software.llnl.gov/Caliper/OutputFormats.html>`_

Caliper and Adiak Options
-------------------------

.. option:: --caliperFilename

            Name of Caliper timing file. Should include file extensions. Optional, default: ``name_of_file_YEAR_MONTH_DATE_TIME.cali``.

.. option:: --caliperConfig CONFIG_STR

            Specify a built-in Caliper configuration or turn off timers with ``none``. Optional, default: ``spot``.

            **Example**:
            ::

               ./spheral ex_prog.py --caliperConfig 'runtime-report(output=time.txt),calc.inclusive,region.count'

.. note::
   The configuration in the example above produces timing results similar to the previous ``Spheral::Timer`` method. This results in a file named ``time.txt`` with cumulative times for the nested regions as well as a count of how many times each region ran.

.. option:: --caliperConfigJSON JSON_FILE

            Specify a JSON file containing a non-default Caliper configuration. Optional.

.. option:: --adiakData ADIAK_DATA_STR

            Specify any Adiak data directly in the command line. Must be a string in key:value format, separated by commas. Optional.

            **Example**:
            ::

               ./spheral ex_prog.py --adiakData "test_name: the_cheat, test_num:10"

.. note::
   By default, all ``commandLine()`` inputs are added as Adiak metadata. ``--adiakData`` are for metadata that does not come through Spheral command line arguments. Adiak metadata can also be added through the python interface. See :ref:`below <python_adiak>` for more details.


Adding Region Timers in C++
---------------------------

The following macros are used to create timing regions in the Spheral C++ interface:

- ``TIME_FUNCTION`` can be added to the very beginning of a function and creates a region for the entire function using the function's name. ``TIME_FUNCTION`` uses just the function name and no class or parameter information, so be careful when using this method with functions that could share names.

- ``TIME_BEGIN("timer_name")`` and ``TIME_END("timer_name")`` create a region between the two different calls and use the string (in this case ``timer_name``) as the name.


Adding Region Timers in Python
------------------------------

Region timers can be added inside the python code using the following function calls:
::

   from SpheralUtilities import TimerMgr
   TimerMgr.timer_start("timer_name")
   some_function_call()
   TimerMgr.timer_end("timer_name")

.. note::
   All timers must have both a start and end call. Otherwise, memory issues will occur.

.. _python_adiak:

Adding Adiak Metadata in Python
-------------------------------

Adiak metadata can be added inside python code using the following function calls:

.. code-block:: python

                adiak_values("value_name", value)

.. _manual_caliper:

Starting Caliper Manually
========================

As mentioned above, the Caliper timing manager is normally configured and started in the ``commandLine()`` routine. However, Caliper can be directly configured and started through the python interface, if desired. This can be done by putting the following into the python file:
::

   from SpheralUtilities import TimerMgr
   caliper_config = "some_configuration(output=some_filename.txt)"
   TimerMgr.add(caliper_config)
   TimerMgr.start()
