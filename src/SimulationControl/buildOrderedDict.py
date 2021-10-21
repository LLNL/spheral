#-------------------------------------------------------------------------------
# Convenient constructor for building an ordered dictionary (remembers the order
# elements were inserted in).
# Once we shift to Python 3 this method will no longer be necessary, as regular
# Python dictionaries now obey this principle.
#-------------------------------------------------------------------------------
import collections

def buildOrderedDict(*args):
    result = collections.OrderedDict()
    for (key, val) in args:
        result[key] = val
    return result
