Diagnostics
###########

Spheral uses Caliper to preform diagnostics, such as timing. To enable this functionality in the code, Spheral needs to be configured with ENABLE_TIMER=On. Otherwise the diagnostic regions are no-ops for improved preformance.
::

  ./scripts/devtools/host-config-build.py <sys_type>-<spec>.cmake -DENABLE_TIMER=On


Querying using Caliper
======================

By defualt, even when configured with ENABLE_TIMER=On, there is no information being recorded. Caliper uses command line options to report data, the simplest of which is ``CALI_CONFIG=runtime-report`` which reports the timing for the regions and prints them out in the terminal. For example:
::

  CALI_CONFIG=runtime-report python Noh-cylindrical-2d.py

.. note::
  To obtain a similar result to that of the removed Spheral::Timer use :kbd:`CALI_CONFIG=runtime\-report(output=time.txt),calc.inclusive,region.count` this will result in a file named time.txt with cumulative times for the nested regions as well as a count of how many times each region ran.

There are many different options that can be used with ``CALI_CONFIG`` to view various information. Here are some extra-links for those who want to read or experiment with other features in Caliper that can be more closely incorperated with Spheral in the future:
::

  https://software.llnl.gov/Caliper/CaliperBasics.html#more-on-configurations
  https://software.llnl.gov/Caliper/BuiltinConfigurations.html
  https://software.llnl.gov/Caliper/configuration.html
  https://software.llnl.gov/Caliper/OutputFormats.html?highlight=cali+query


Adding Regions
==============

So far there are two different types of regions in Spheral, using the following macros:
::

  TIME_FUNCTION
  TIME_BEGIN("timer_name")
  TIME_END("timer_name")

- ``TIME_FUNCTION`` can be added to the very beginning of a function and creates a region for the entire function using the function's name. ``TIME_FUNCTION`` uses just the function name and no class or parameter information, so be careful when using this method with functions that could share names.

- ``TIME_BEGIN("timer_name")`` and ``TIME_END("timer_name")`` create a region between the two different calls and use the string (in this case timer_name) as the name. This is the most similar method to the removed Spheral::Timer.
