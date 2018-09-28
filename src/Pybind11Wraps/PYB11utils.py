from PYB11Decorators import *
import inspect, StringIO

#-------------------------------------------------------------------------------
# PYB11inject
#
# Add methods defined in fromclass to tocls
#-------------------------------------------------------------------------------
def PYB11inject(fromcls, tocls,
                virtual = None,
                pure_virtual = None):
    assert not (virtual and pure_virtual), "PYB11inject: cannot specify both virtual and pure_virtual as True!"
    names = [x for x in dir(fromcls) if (x[:2] != "__" and inspect.ismethod(eval('fromcls.%s' % x)))]
    for name in names:
        exec('''tocls.%(name)s = eval("fromcls.__dict__['%(name)s']")''' % {"name": name})
        if not virtual is None:
            exec('tocls.%s.__dict__["PYB11virtual"] = %s' % (name, virtual))
        if not pure_virtual is None:
            exec('tocls.%s.__dict__["PYB11pure_virtual"] = %s' % (name, pure_virtual))
    return

#-------------------------------------------------------------------------------
# Key function to sort lists by source code order.
#-------------------------------------------------------------------------------
def PYB11sort_by_line(stuff):
    from PYB11class import PYB11TemplateClass
    name, obj = stuff
    if isinstance(obj, PYB11TemplateClass):
        return obj.order + 0
    else:
        try:
            source, lineno = inspect.findsource(obj)
        except:
            raise RuntimeError, "Cannot find source for %s?" % name
        return lineno

#-------------------------------------------------------------------------------
# PYB11classes
#
# Get the classes to bind from a module
#-------------------------------------------------------------------------------
def PYB11classes(modobj):
    result = [(name, cls) for (name, cls) in inspect.getmembers(modobj, predicate=inspect.isclass)
              if name[:5] != "PYB11"]
    return sorted(result, key = PYB11sort_by_line)

#-------------------------------------------------------------------------------
# PYB11classTemplateInsts
#
# Get the template class instantiations to bind from a module
#-------------------------------------------------------------------------------
def PYB11classTemplateInsts(modobj):
    from PYB11class import PYB11TemplateClass
    result = [x for x in dir(modobj) if isinstance(eval("modobj.%s" % x), PYB11TemplateClass)]
    result = [(x, eval("modobj.%s" % x)) for x in result]
    return sorted(result, key = PYB11sort_by_line)

#-------------------------------------------------------------------------------
# PYB11ClassMethods
#
# Get the methods to bind from a class
#-------------------------------------------------------------------------------
def PYB11ClassMethods(obj):
    result = inspect.getmembers(obj, predicate=inspect.ismethod)
    # It's nice to sort in the same order the user created, but not necessary
    try:
        result.sort(key = PYB11sort_by_line)
    except:
        pass
    return result

#-------------------------------------------------------------------------------
# PYB11ThisClassMethods
#
# Cull the methods found in PYB11ClassMethods to just those defined locally
# in obj.
#-------------------------------------------------------------------------------
def PYB11ThisClassMethods(obj):
    result = PYB11ClassMethods(obj)
    return [(name, meth) for (name, meth) in result if name in obj.__dict__]

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
# Mangle a string to a safe C++ variable name.
#-------------------------------------------------------------------------------
def PYB11mangle(name):
    result = name.replace("<", "__").replace(">", "__").replace("::", "_").replace(", ", "_").replace(",", "_")
    return result

#-------------------------------------------------------------------------------
# Union of dictionarys.
#-------------------------------------------------------------------------------
def PYB11union_dict(*args):
    result = {}
    for d in args:
        for key in d:
            result[key] = d[key]
    return result

#-------------------------------------------------------------------------------
# PYB11CPPsafe
#
# Mangle a string to make commas safe for CPP directives.
#-------------------------------------------------------------------------------
def PYB11CPPsafe(string):
    return string.replace(",", " PYB11COMMA ")

#-------------------------------------------------------------------------------
# PYB11cppname_exts
#
# Return the C++ template <...> description, and a mangled string thereof.
#-------------------------------------------------------------------------------
def PYB11cppname_exts(templateargs):
    tt, mt = "", ""
    if templateargs:
        tt = "<"
        for i, arg in enumerate(templateargs):
            if i < len(templateargs) - 1:
                tt += "%s," % arg
            else:
                tt += "%s>" % arg
        mt = PYB11mangle(tt)
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
         "implementation" : None,
         "returnpolicy"   : None,
         "template"       : (),
         "template_dict"  : {}}
    for key in d:
        if hasattr(obj, "PYB11" + key):
            d[key] = eval("obj.PYB11%s" % key)
    safeexts= PYB11cppname_exts(d["template"])
    d["full_cppname"] = d["cppname"] + safeexts[0]
    d["mangle_cppname"] = d["cppname"] + safeexts[1]
    return d
