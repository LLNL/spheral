import sys
from pybindgen import *

#-------------------------------------------------------------------------------
# Add an object by name to a module/namespace, and publish it to the world.
#-------------------------------------------------------------------------------
def addObject(mod, name, *args, **kwargs):
    x = mod.add_class(name, *args, **kwargs)
    return x

#-------------------------------------------------------------------------------
# Find the wrapped object in a module.
#-------------------------------------------------------------------------------
def findObject(scope, name):
    for stuff in scope.items():
        assert len(stuff) == 2
        if stuff[0] == name:
            return stuff[1]
    raise RuntimeError, "Unable to find %s in the specified scope." % name

#-------------------------------------------------------------------------------
# For the given pybindgen module, return the (headers, spaces, cppclass')
#-------------------------------------------------------------------------------
def extractPybindgenObjs(mod):
    includes, spaces, objs = [], [], {}
    
    def extractSubmoduleStuff(submod):
        includes += submod.includes
        if space.cpp_namespace:
            spaces.append(cpp_namespaces)
    

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
