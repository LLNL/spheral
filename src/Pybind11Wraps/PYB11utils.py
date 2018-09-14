import inspect

#-------------------------------------------------------------------------------
# PYB11parseArgs
#
# Return (argType, argName, <default_value>)
#-------------------------------------------------------------------------------
def PYB11parseArgs(meth):
    stuff = inspect.getargspec(meth)
    result = []
    if stuff.defaults:
        nargs = len(stuff.defaults)
        for argName, val in zip(stuff.args[-nargs:], stuff.defaults):
            if isinstance(val, tuple):
                assert len(val) == 2
                argType, default = val
            else:
                argType, default = val, None
            result.append((argType, argName, default))
    return result

#-------------------------------------------------------------------------------
# PYB11attrs
#
# Read the possible PYB11 generation attributes from the obj
#-------------------------------------------------------------------------------
def PYB11attrs(obj):
    d = {"pyname"       : obj.__name__,
         "cppname"      : obj.__name__,
         "ignore"       : False,
         "namespace"    : "",
         "singleton"    : False,
         "virtual"      : False,
         "pure_virtual" : False,
         "const"        : False,
         "static"       : False,
         "readwrite"    : False,          # Attribute
         "readonly"     : False}          # Attribute
    for key in d:
        if hasattr(obj, "PYB11" + key):
            d[key] = eval("obj.PYB11%s" % key)
    return d
