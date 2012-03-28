import sys
from pybindgen import *

## #-------------------------------------------------------------------------------
## # Find the object associated with the given name.
## #-------------------------------------------------------------------------------
## def findObject(name):
##     for mod in sys.modules.values():
##         try:
##             for xname, val in mod.__dict__.items():
##                 if xname == name:
##                     return val
##         except:
##             pass
##     raise ValueError, "Could not find %s" % name

## #-------------------------------------------------------------------------------
## # Find a set of objects.
## #-------------------------------------------------------------------------------
## def findObjects(*args):
##     return [findObject(name) for name in args]

## #-------------------------------------------------------------------------------
## # Publish an object to the world.
## #-------------------------------------------------------------------------------
## def publish(thing, name):
##     globals()[name] = thing
##     if not publishedRegistry in globals():
##         globals()[publishedRegistry] = []
##     globals()[publishedRegistry].append(name)

## def publishObjects(stuff):
##     for thing, name in stuff:
##         publish(thing, name)

## #-------------------------------------------------------------------------------
## # Load all the currently published types from Spheral.
## #-------------------------------------------------------------------------------
## def findAllPublishedObjects():
##     names = globals()[publishedRegistry]
##     result = [findObject(name) for name in names]
##     return result

#-------------------------------------------------------------------------------
# Add an object by name to a module/namespace, and publish it to the world.
#-------------------------------------------------------------------------------
def addObject(mod, name, *args, **kwargs):
    x = mod.add_class(name, *args, **kwargs)
    if not 'wrapObjs' in mod.__dict__:
        mod.wrapObjs = {}
    if "custom_template_class_name" in kwargs:
        pubname = kwargs["custom_template_class_name"]
    elif "python_name" in kwargs:
        pubname = kwargs["python_name"]
    else:
        pubname = name
    mod.wrapObjs[pubname] = x
    return x

#-------------------------------------------------------------------------------
# Add a reference symbol to a type.
#-------------------------------------------------------------------------------
def ref(name):
    return "%s&" % name

#-------------------------------------------------------------------------------
# Add a pointer symbol to a type.
#-------------------------------------------------------------------------------
def ptr(name):
    return "%s*" % name

def const_ptr(name):
    return "const %s*" % name

#-------------------------------------------------------------------------------
# Return the normal way we want to handle a pointer parameter.
#-------------------------------------------------------------------------------
def inputptrparam(cppobj, argname):
    return Parameter.new(ptr(cppobj), argname, transfer_ownership=False) # , direction=Parameter.DIRECTION_IN)

#-------------------------------------------------------------------------------
# Return the normal way we want to handle a reference parameter.
#-------------------------------------------------------------------------------
def refparam(cppobj, argname, default_value=None):
    return Parameter.new(ref(cppobj), argname, direction=Parameter.DIRECTION_INOUT, default_value=default_value)

def constrefparam(cppobj, argname, default_value=None):
    return Parameter.new(ref("const " + cppobj), argname, direction=Parameter.DIRECTION_INOUT, default_value=default_value)
