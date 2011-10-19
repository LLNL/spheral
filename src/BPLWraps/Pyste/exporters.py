# Copyright Bruno da Silva de Oliveira 2003. Use, modification and 
# distribution is subject to the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at 
# http:#www.boost.org/LICENSE_1_0.txt)


# a list of Exporter instances
exporters = []

# JMO begin
# a list of C++ names that we claim are exported somewhere else
external_exported_names = []
# JMO end

current_interface = None # the current interface file being processed
importing = False    # whetever we are now importing a pyste file.
                     # exporters created here shouldn't export themselves
