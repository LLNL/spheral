from PYB11Decorators import *
import inspect, StringIO, types, itertools

#-------------------------------------------------------------------------------
# PYB11inject
#
# Add methods defined in fromclass to tocls
#-------------------------------------------------------------------------------
# We need this bit of trickery from stackoverflow to make a new copy of the method
# in question.  Necessary to avoid corrupting attributes of the original classes
# methods.
def PYB11copy_func(f, name=None):
    '''
    return a function with same code, globals, defaults, closure, and 
    name (or provide a new name)
    '''
    fn = types.FunctionType(f.__code__, f.__globals__, name or f.__name__,
                            f.__defaults__, f.__closure__)
    # in case f was given attrs (note this dict is a shallow copy):
    fn.__dict__.update(f.__dict__) 
    return fn

def PYB11inject(fromcls, tocls,
                virtual = None,
                pure_virtual = None):
    assert not (virtual and pure_virtual), "PYB11inject: cannot specify both virtual and pure_virtual as True!"

    # Methods
    names = [x for x in dir(fromcls) if (inspect.ismethod(eval('fromcls.%s' % x)))]
    for name in names:
        exec('''tocls.%(name)s = PYB11copy_func(fromcls.%(name)s)''' % {"name": name})
        #exec('''tocls.%(name)s = copy_func(eval("fromcls.__dict__['%(name)s']"))''' % {"name": name})
        if not virtual is None:
            exec('tocls.%s.__dict__["PYB11virtual"] = %s' % (name, virtual))
        if not pure_virtual is None:
            exec('tocls.%s.__dict__["PYB11pure_virtual"] = %s' % (name, pure_virtual))

    # Properties
    from PYB11class import PYB11TemplateMethod
    names = [x for x in dir(fromcls) if isinstance(eval('fromcls.%s' % x), PYB11TemplateMethod)]
    for name in names:
        exec('''tocls.%(name)s = PYB11TemplateMethod(func_template = fromcls.%(name)s.func_template,
                                                     template_parameters = [x[1] for x in fromcls.%(name)s.template_parameters],
                                                     cppname = fromcls.%(name)s.cppname,
                                                     pyname = fromcls.%(name)s.pyname,
                                                     docext = fromcls.%(name)s.docext)''' % {"name": name})

    # Properties
    from PYB11property import PYB11property
    names = [x for x in dir(fromcls) if isinstance(eval('fromcls.%s' % x), PYB11property)]
    for name in names:
        exec('''tocls.%(name)s = PYB11property(returnType = fromcls.%(name)s.returnType,
                                               getter = fromcls.%(name)s.getter,
                                               setter = fromcls.%(name)s.setter,
                                               doc = fromcls.%(name)s.doc,
                                               getterraw = fromcls.%(name)s.getterraw,
                                               setterraw = fromcls.%(name)s.setterraw,
                                               getterconst = fromcls.%(name)s.getterconst,
                                               setterconst = fromcls.%(name)s.setterconst,
                                               static = fromcls.%(name)s.static,
                                               returnpolicy = fromcls.%(name)s.returnpolicy)''' % {"name": name})

    # Attributes
    from PYB11ClassAttribute import PYB11ClassAttribute
    names = [x for x in dir(fromcls) if isinstance(eval('fromcls.%s' % x), PYB11ClassAttribute)]
    for name in names:
        exec('''tocls.%(name)s = PYB11ClassAttribute(static = fromcls.%(name)s.static,
                                                     pyname = fromcls.%(name)s.pyname,
                                                     cppname = fromcls.%(name)s.cppname,
                                                     doc = fromcls.%(name)s.doc,
                                                     deftype = fromcls.%(name)s.deftype)''' % {"name": name})

    return

#-------------------------------------------------------------------------------
# Return the base classes of a class
#
# Computes a dictionary giving the direct bases for all classes in the
# inheritance hierarchy of a class.
#-------------------------------------------------------------------------------
def PYB11getBaseClasses(klass):
    stuff = inspect.getclasstree(inspect.getmro(klass), unique=True)
    def flatten(s, result):
        if type(s) is list:
            for val in s:
                s = flatten(val, result)
        else:
            result.append(s)
    flatstuff = []
    flatten(stuff, flatstuff)
    result = { k[0] : k[1] for k in flatstuff }
    return result

#-------------------------------------------------------------------------------
# Key function to sort lists by source code order.
#-------------------------------------------------------------------------------
def PYB11sort_by_line(stuff):
    from PYB11class import PYB11TemplateClass
    name, obj = stuff
    if isinstance(obj, PYB11TemplateClass):
        #return obj.order + 0
        try:
            source, lineno = inspect.findsource(obj.klass_template)
        except:
            raise RuntimeError, "Cannot find source for %s?" % name
        #print " **> ", name, lineno
        return lineno
    else:
        try:
            source, lineno = inspect.findsource(obj)
        except:
            raise RuntimeError, "Cannot find source for %s?" % name
        #print " ==> ", name, lineno
        return lineno

#-------------------------------------------------------------------------------
# PYB11sort_by_inheritance
#
# Key sorting function to put base classes first.
#-------------------------------------------------------------------------------
class PYB11sort_by_inheritance:
    def __init__(self, klasses):
        from PYB11class import PYB11TemplateClass

        # First pass, order by line number
        self.keys = {}
        for (name, obj) in klasses:
            if isinstance(obj, PYB11TemplateClass):
                klass = obj.klass_template
            else:
                klass = obj
            self.keys[klass] = PYB11sort_by_line((name, klass))

        # Now make sure classes come after any of their bases
        changed = True
        while changed:
            changed = False
            for (name, obj) in klasses:
                if isinstance(obj, PYB11TemplateClass):
                    klass = obj.klass_template
                else:
                    klass = obj
                for bklass in inspect.getmro(klass)[1:]:
                    if self.keys[klass] <= self.keys[bklass]:
                        self.keys[klass] = self.keys[bklass] + 1
                        changed = True

    def __call__(self, stuff):
        from PYB11class import PYB11TemplateClass
        obj = stuff[1]
        if isinstance(obj, PYB11TemplateClass):
            klass = obj.klass_template
        else:
            klass = obj
        return self.keys[klass]

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
# PYB11othermods
#
# Get the modules we should import if any
#-------------------------------------------------------------------------------
def PYB11othermods(modobj):
    if hasattr(modobj, "import_modules"):
        return modobj.import_modules
    else:
        return []

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
    # It's nice to sort in the same order the user created, but not necessary
    try:
        result.sort(key = PYB11sort_by_line)
    except:
        pass
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
    allmethods = PYB11ClassMethods(klass)
    virtual = False
    for mname, meth in allmethods:
        methattrs = PYB11attrs(meth)
        if methattrs["virtual"] or methattrs["pure_virtual"]:
            virtual = True
    return virtual

#-------------------------------------------------------------------------------
# PYB11protectedClass
#
# Test if the given class has protected methods.
#-------------------------------------------------------------------------------
def PYB11protectedClass(klass):
    klassattrs = PYB11attrs(klass)
    allmethods = PYB11ThisClassMethods(klass)
    protected = False
    for mname, meth in allmethods:
        methattrs = PYB11attrs(meth)
        if methattrs["protected"]:
            protected = True
    return protected

#-------------------------------------------------------------------------------
# PYB11badchars
#
# Check if any of the forbidden characters for PYBIND11_OVERLOAD* are present.
#-------------------------------------------------------------------------------
def PYB11badchars(name):
    return any((c in ("<", ">", ",")) for c in name)

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
# PYB11docstring
#
# Generate a reasonably formatted doc string
#-------------------------------------------------------------------------------
def PYB11docstring(doc, ss):
    if doc:
        stuff = doc.split("\n")
        if len(stuff) == 1:
            ss('"%s"' % doc.replace('"', '\\"'))
        else:
            ss("\n")
            for i, line in enumerate(doc.split('\n')):
                ss('            "%s\\n"' % line.replace('"', '\\"'));
                if i < len(stuff) - 1:
                    ss("\n")
    return

#-------------------------------------------------------------------------------
# PYB11attrs
#
# Read the possible PYB11 generation attributes from the obj
#-------------------------------------------------------------------------------
def PYB11attrs(obj):
    d = {"pyname"                : obj.__name__,
         "cppname"               : obj.__name__,
         "ignore"                : False,
         "namespace"             : "",
         "singleton"             : False,
         "holder"                : None,
         "exposeBaseOverloads"   : True,
         "dynamic_attr"          : None,
         "virtual"               : False,
         "pure_virtual"          : False,
         "protected"             : False,
         "const"                 : False,
         "static"                : False,
         "implementation"        : None,
         "returnpolicy"          : None,
         "keepalive"             : None,
         "template"              : (),
         "template_dict"         : {},
         "module"                : {}}
    for key in d:
        if hasattr(obj, "PYB11" + key):
            d[key] = eval("obj.PYB11%s" % key)
    safeexts= PYB11cppname_exts(d["template"])
    d["full_cppname"] = d["cppname"] + safeexts[0]
    d["mangle_cppname"] = d["cppname"] + safeexts[1]
    return d
