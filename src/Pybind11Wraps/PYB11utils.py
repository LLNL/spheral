import inspect, StringIO

#-------------------------------------------------------------------------------
# Key function to sort lists by source code order.
#-------------------------------------------------------------------------------
def PYB11sort_by_line(stuff):
    name, obj = stuff
    source, lineno = inspect.findsource(obj)
    return lineno

#-------------------------------------------------------------------------------
# PYB11classes
#
# Get the classes to bind from a module
#-------------------------------------------------------------------------------
def PYB11classes(modobj):
    result = [(name, cls) for (name, cls) in inspect.getmembers(modobj, predicate=inspect.isclass)
              if name[:5] != "PYB11"]
    result.sort(key = PYB11sort_by_line)
    return result

#-------------------------------------------------------------------------------
# PYB11classTemplateInsts
#
# Get the template class instantiations to bind from a module
#-------------------------------------------------------------------------------
def PYB11classTemplateInsts(modobj):
    from PYB11class import PYB11TemplateClass
    result = [x for x in dir(modobj) if isinstance(eval("modobj.%s" % x), PYB11TemplateClass)]
    result = [(x, eval("modobj.%s" % x)) for x in result]
    return result

#-------------------------------------------------------------------------------
# PYB11ClassMethods
#
# Get the methods to bind from a class
#-------------------------------------------------------------------------------
def PYB11ClassMethods(obj):
    result = inspect.getmembers(obj, predicate=inspect.ismethod)
    result.sort(key = PYB11sort_by_line)
    return result

#-------------------------------------------------------------------------------
# PYB11functions
#
# Get the functions to bind from a module
#-------------------------------------------------------------------------------
def PYB11functions(modobj):
    result = [(name, meth) for (name, meth) in inspect.getmembers(modobj, predicate=inspect.isfunction)
              if name[:5] != "PYB11"]
    result.sort(key = PYB11sort_by_line)
    return result

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
# PYB11virtualClass
#
# Test if the given class has virtual methods.
#-------------------------------------------------------------------------------
def PYB11virtualClass(klass):
    klassattrs = PYB11attrs(klass)
    allmethods = [(mname, meth) for (mname, meth) in PYB11ClassMethods(klass)
                  if not PYB11attrs(meth)["ignore"]]
    virtual = False
    for mname, meth in allmethods:
        methattrs = PYB11attrs(meth)
        if methattrs["virtual"] or methattrs["pure_virtual"]:
            virtual = True
    return virtual

#-------------------------------------------------------------------------------
# PYB11mangle
#
# Return the C++ template <...> description, and a mangled string thereof.
#-------------------------------------------------------------------------------
def PYB11mangle(templateargs):
    tt = "<"
    for i, arg in enumerate(templateargs):
        if i < len(templateargs) - 1:
            tt += "%s, " % arg
        else:
            tt += "%s>" % arg
    mt = tt.replace("<", "__").replace(">", "__").replace("::", "_").replace(", ", "_")
    return tt, mt

#-------------------------------------------------------------------------------
# PYB11indentedIO
#
# Add extra indentation to an output stream.
#-------------------------------------------------------------------------------
class PYB11indentedIO:
    def __init__(self, prefix):
        self.prefix = prefix
        self.fs = StringIO.StringIO()
        return
    def __call__(self, stuff):
        newstuff = stuff.replace("\n", "\n" + self.prefix)
        self.fs.write(newstuff)
        return
    def getvalue(self):
        return self.fs.getvalue()
    def close(self):
        self.fs.close()

#-------------------------------------------------------------------------------
# PYB11attrs
#
# Read the possible PYB11 generation attributes from the obj
#-------------------------------------------------------------------------------
def PYB11attrs(obj):
    d = {"pyname"         : obj.__name__,
         "cppname"        : obj.__name__,
         "ignore"         : False,
         "namespace"      : "",
         "singleton"      : False,
         "virtual"        : False,
         "pure_virtual"   : False,
         "const"          : False,
         "static"         : False,
         "readwrite"      : False,          # Attribute
         "readonly"       : False,          # Attribute
         "implementation" : None,
         "template"       : (),
         "template_dict"  : {}}
    for key in d:
        if hasattr(obj, "PYB11" + key):
            d[key] = eval("obj.PYB11%s" % key)
    return d
