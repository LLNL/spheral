#--------------------------------------------------------------------------------
# PYB11class
#
# Stuff for handling classes in pybind11
#--------------------------------------------------------------------------------
from PYB11utils import *
from PYB11property import *
from PYB11ClassAttribute import *
from PYB11Trampoline import *
from PYB11enum import PYB11enum
import copy, StringIO, inspect

#-------------------------------------------------------------------------------
# PYB11generateModuleClasses
#
# Bind the classes in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleClasses(modobj, ss):
    klasses = PYB11classes(modobj) + PYB11classTemplateInsts(modobj)
    klasses = sorted(klasses, key=PYB11sort_by_inheritance(klasses))
    for kname, klass in klasses:
        if isinstance(klass, PYB11TemplateClass):
            klass(kname, ss)
        else:
            klassattrs = PYB11attrs(klass)
            mods = klassattrs["module"]
            if ((not klassattrs["ignore"]) and                                   # ignore this class?
                ((klass not in mods) or mods[klass] == modobj.PYB11modulename)): # is this class imported from another mod?
                PYB11generateClass(klass, klassattrs, ss)

#--------------------------------------------------------------------------------
# Make a class template instantiation
#--------------------------------------------------------------------------------
class PYB11TemplateClass:
    __order = 0

    def __init__(self,
                 klass_template,
                 template_parameters,
                 cppname = None,
                 pyname = None,
                 docext = ""):
        self.klass_template = klass_template
        self.cppname = cppname
        self.pyname = pyname
        self.docext = docext

        # Create the template parameter dictionary
        self.template_parameters = {}
        klassattrs = PYB11attrs(self.klass_template)
        if isinstance(template_parameters, str):
            assert len(klassattrs["template"]) == 1
            self.template_parameters[klassattrs["template"][0]] = template_parameters
        elif isinstance(template_parameters, tuple):
            assert len(klassattrs["template"]) == len(template_parameters)
            for name, val in zip(klassattrs["template"], template_parameters):
                self.template_parameters[name] = val
        else:
            assert isinstance(template_parameters, dict)
            for key in klassattrs["template"]:
                if not key in template_parameters:
                    raise RuntimeError, "Template parameter dictionary spec error: %s is missing from %s" % (key, template_parameters)
            self.template_parameters = template_parameters
            
        # Record the order of instantiations
        self.order = PYB11TemplateClass.__order + 1
        PYB11TemplateClass.__order += 1
        return

    def __call__(self, pyname, ss):
        klassattrs = self.mangleNames(pyname)
        if self.klass_template.__doc__:
            doc0 = copy.deepcopy(self.klass_template.__doc__)
            self.klass_template.__doc__ += self.docext
        PYB11generateClass(self.klass_template, klassattrs, ss)
        if self.klass_template.__doc__:
            self.klass_template.__doc__ = doc0
        return

    def makeTrampoline(self, pyname, ss):
        klassattrs = self.mangleNames(pyname)
        PYB11generateTrampoline(self.klass_template, klassattrs, ss)
        return

    # Do some template mangling (and magically put the template parameters in scope).
    def mangleNames(self, pyname):
        klassattrs = PYB11attrs(self.klass_template)
        template_ext = "<"
        doc_ext = ""
        for name in klassattrs["template"]:
            val = self.template_parameters[name]
            exec("%s = '%s'" % (name, val))
            template_ext += "%s, " % val
            doc_ext += "_%s_" % val.replace("::", "_").replace("<", "_").replace(">", "_")
        template_ext = template_ext[:-2] + ">"

        if self.cppname:
            klassattrs["cppname"] = self.cppname
        else:
            klassattrs["cppname"] += template_ext
        if self.pyname:
            klassattrs["pyname"] = self.pyname
        else:
            klassattrs["pyname"] = pyname

        klassattrs["template_dict"] = self.template_parameters
        return klassattrs

#-------------------------------------------------------------------------------
# Make a class method template instantiation
#-------------------------------------------------------------------------------
class PYB11TemplateMethod:

    def __init__(self,
                 func_template,
                 template_parameters,
                 cppname = None,
                 pyname = None,
                 docext = ""):
        if isinstance(template_parameters, str):
            template_parameters = (template_parameters,)
        self.func_template = func_template
        funcattrs = PYB11attrs(self.func_template)
        assert len(funcattrs["template"]) == len(template_parameters)
        self.template_parameters = [(name, val) for (name, val) in zip(funcattrs["template"], template_parameters)]
        self.cppname = cppname
        self.pyname = pyname
        self.docext = docext
        return

    def __call__(self, pyname, klass, klassattrs, ss):
        # Do some template mangling (and magically put the template parameters in scope).
        template_ext = "<"
        doc_ext = ""

        for name, val in self.template_parameters:
            exec("%s = '%s'" % (name, val))
            template_ext += "%s, " % val
            doc_ext += "_%s_" % val.replace("::", "_").replace("<", "_").replace(">", "_")
        template_ext = template_ext[:-2] + ">"

        funcattrs = PYB11attrs(self.func_template)
        if self.cppname:
            funcattrs["cppname"] = self.cppname
        else:
            funcattrs["cppname"] += template_ext
        if self.pyname:
            funcattrs["pyname"] = self.pyname
        else:
            funcattrs["pyname"] = pyname

        funcattrs["template_dict"] = {}
        for name, val in self.template_parameters:
            funcattrs["template_dict"][name] = val

        if self.func_template.__doc__:
            doc0 = copy.deepcopy(self.func_template.__doc__)
            self.func_template.__doc__ += self.docext
        fs = StringIO.StringIO()
        PYB11generic_class_method(klass, klassattrs, self.func_template, funcattrs, fs.write)
        ss(fs.getvalue() % PYB11union_dict(klassattrs["template_dict"], funcattrs["template_dict"]))
        fs.close()
        if self.func_template.__doc__:
            self.func_template.__doc__ = doc0
        return
        

#-------------------------------------------------------------------------------
# Generic class method generation
#-------------------------------------------------------------------------------
def PYB11generic_class_method(klass, klassattrs, meth, methattrs, ss):
    klassinst = klass()
    args = PYB11parseArgs(meth)
    methattrs["returnType"] = meth(klassinst)

    methattrs["namespace"] = "%(namespace)s" % klassattrs
    methattrs["classcppname"] = "%(cppname)s" % klassattrs
    if methattrs["protected"]:
        methattrs["classcppname"] = "PYB11Publicist" + methattrs["classcppname"]

    if methattrs["static"]:
        ss('    obj.def_static("%(pyname)s", ' % methattrs)
    else:
        ss('    obj.def("%(pyname)s", ' % methattrs)

    # Check for argument specs
    argString = ""
    for i, (argType, argName, default) in enumerate(args):
        argString += ', "%s"_a' % argName
        if methattrs["noconvert"]:
            argString += '.noconvert()'
        if default:
            argString += "=%s" % default

    # If there is an implementation, short-circuit the rest.
    if methattrs["implementation"]:
        ss(methattrs["implementation"] + argString)

    elif methattrs["returnType"] is None:
        if methattrs["static"]:
            ss("&%(namespace)s%(cppname)s" % methattrs)
        else:
            ss("&%(namespace)s%(classcppname)s::%(cppname)s" % methattrs)
    else:
        ss("(%(returnType)s " % methattrs)
        if methattrs["static"]:
            ss("(%(namespace)s*)(" % methattrs)
        else:
            ss("(%(namespace)s%(cppname)s::*)(" % klassattrs)
        for i, (argType, argName, default) in enumerate(args):
            ss(argType)
            if i < len(args) - 1:
                ss(", ")
        if methattrs["const"]:
            ss(") const) &%(namespace)s%(classcppname)s::%(cppname)s" % methattrs + argString)
        else:
            ss(")) &%(namespace)s%(classcppname)s::%(cppname)s" % methattrs + argString)

    # Is there a return value policy?
    if methattrs["returnpolicy"]:
        ss(", py::return_value_policy::%s" % methattrs["returnpolicy"])

    # Is there a call guard?
    if methattrs["call_guard"]:
        ss(", py::call_guard<%s>()" % methattrs["call_guard"])

    # Is there a keep_alive policy?
    if methattrs["keepalive"]:
        assert isinstance(methattrs["keepalive"], tuple)
        assert len(methattrs["keepalive"]) == 2
        ss(", py::keep_alive<%i, %i>()" % methattrs["keepalive"])

    # Write the doc string
    doc = inspect.getdoc(meth)
    if doc:
        ss(", ")
        PYB11docstring(doc, ss)
    ss(");\n")

#-------------------------------------------------------------------------------
# PYB11generateClass
#
# Bind the methods for the given class
#-------------------------------------------------------------------------------
def PYB11generateClass(klass, klassattrs, ssout):
    klassinst = klass()

    fs = StringIO.StringIO()
    ss = fs.write

    #...........................................................................
    # Ignore a method
    def ignore(mesh, methattrs, args):
        pass

    #...........................................................................
    # pyinit<>
    def pyinit(meth, methattrs, args):
        if methattrs["implementation"]:
            ss("    obj.def(py::init(%(implementation)s)" % methattrs)
        else:
            ss("    obj.def(py::init<")
            argString = ""
            for i, (argType, argName, default) in enumerate(args):
                if i < len(args) - 1:
                    ss("%s, " % argType)
                else:
                    ss("%s" % argType)
                argString += ', "%s"_a' % argName
                if methattrs["noconvert"]:
                    argString += '.noconvert()'
                if default:
                    argString += "=%s" % default
            ss(">()%s" % argString)
        doc = inspect.getdoc(meth)
        if doc:
            ss(", ")
            PYB11docstring(doc, ss)
        ss(");\n")


        return

    #...........................................................................
    # Binary operators
    def binary_operator(meth, methattrs, args, op):
        assert len(args) in (0, 1)
        if len(args) == 0:
            argType = "py::self"
        else:
            argType = args[0][0]
        ss('    obj.def(py::self %s %s);\n' % (op, argType))

    #...........................................................................
    # Reverse binary operators
    def reverse_binary_operator(meth, methattrs, args, op):
        assert len(args) in (0, 1)
        if len(args) == 0:
            argType = "py::self"
        else:
            argType = args[0][0]
        ss('    obj.def(%s %s py::self);\n' % (argType, op))

    #...........................................................................
    # Unary operators
    def unary_operator(meth, methattrs, args, op):
        assert len(args) == 0
        ss('    obj.def(%spy::self);\n' % op)

    #...........................................................................
    # operator()
    def call_operator(meth, methattrs, args, op):
        ss('    obj.def("__call__", ')
        if methattrs["implementation"]:
            ss(methattrs["implementation"])
        elif methattrs["returnType"] is None:
            ss("&%(namespace)s%(cppname)s::operator()" % klassattrs)
        else:
            argString = ""
            ss(("(%(returnType)s " % methattrs) + ("(%(namespace)s%(cppname)s::*)(" % klassattrs))
            for i, (argType, argName, default) in enumerate(args):
                ss(argType)
                if i < len(args) - 1:
                    ss(", ")
                argString += ', "%s"_a' % argName
                if methattrs["noconvert"]:
                    argString += '.noconvert()'
                if default:
                    argString += "=%s" % default
            if methattrs["const"]:
                ss((") const) &%(namespace)s%(cppname)s::operator()" % klassattrs) + argString)
            else:
                ss((")) &%(namespace)s%(cppname)s::operator()" % klassattrs) + argString)
        doc = inspect.getdoc(meth)
        if doc:
            ss(", ")
            PYB11docstring(doc, ss)
        ss(", py::is_operator());\n")

    #...........................................................................
    # Tabulate the dispatch for special operations.
    special_operators =  {"__init__": (ignore, ""),

                          "__call__": (call_operator, ""),

                          "__add__" : (binary_operator, "+"),
                          "__sub__" : (binary_operator, "-"),
                          "__mul__" : (binary_operator, "*"),
                          "__div__" : (binary_operator, "/"),
                          "__mod__" : (binary_operator, "%"),
                          "__and__" : (binary_operator, "&"),
                          "__xor__" : (binary_operator, "^"),
                          "__or__"  : (binary_operator, "|"),
                          
                          "__radd__" : (reverse_binary_operator, "+"),
                          "__rsub__" : (reverse_binary_operator, "-"),
                          "__rmul__" : (reverse_binary_operator, "*"),
                          "__rdiv__" : (reverse_binary_operator, "/"),
                          "__rmod__" : (reverse_binary_operator, "%"),
                          "__rand__" : (reverse_binary_operator, "&"),
                          "__rxor__" : (reverse_binary_operator, "^"),
                          "__ror__"  : (reverse_binary_operator, "|"),
                          
                          "__iadd__" : (binary_operator, "+="),
                          "__isub__" : (binary_operator, "-="),
                          "__imul__" : (binary_operator, "*="),
                          "__idiv__" : (binary_operator, "/="),
                          "__imod__" : (binary_operator, "%="),
                          "__iand__" : (binary_operator, "&="),
                          "__ixor__" : (binary_operator, "^="),
                          "__ior__"  : (binary_operator, "|="),

                          "__neg__"    : (unary_operator, "-"),
                          "__invert__" : (unary_operator, "~"),

                          "__lt__" : (binary_operator, "<"),
                          "__le__" : (binary_operator, "<="),
                          "__eq__" : (binary_operator, "=="),
                          "__ne__" : (binary_operator, "!="),
                          "__gt__" : (binary_operator, ">"),
                          "__ge__" : (binary_operator, ">=")}

    # Start generating.
    ss("""
  //............................................................................
  // Class %(pyname)s
  {
""" % klassattrs)
    # If the class has specified any typedefs, do them.
    if hasattr(klass, "PYB11typedefs"):
        ss(klass.PYB11typedefs + "\n")

    ss("    py::class_<%(namespace)s%(cppname)s" % klassattrs)

    # Check for base classes.
    cppname = "%(namespace)s%(cppname)s" % klassattrs
    bklasses = PYB11getBaseClasses(klass)
    for bklass in bklasses[klass]:
        bklassattrs = PYB11attrs(bklass)
        bcppname = "%(namespace)s%(cppname)s" % bklassattrs
        if bklassattrs["template"]:
            bcppname += "<"
            for i, t in enumerate(bklassattrs["template"]):
                if i < len(bklassattrs["template"]) - 1:
                    bcppname += ("%(" + t + ")s, ")
                else:
                    bcppname += ("%(" + t + ")s>")
            bcppname = bcppname % klassattrs["template_dict"]
        if bcppname != cppname:
            ss(", " + bcppname)

    # Any trampoline?
    if PYB11virtualClass(klass):
        ss(", %(namespace)sPYB11Trampoline%(cppname)s" % klassattrs)

    # Is this a singleton?
    if klassattrs["singleton"]:
        ss(", std::unique_ptr<%(namespace)s%(cppname)s, py::nodelete>" % klassattrs)

    # Did we specify a holder type?
    elif klassattrs["holder"]:
        ss(", %(holder)s<%(namespace)s%(cppname)s>" % klassattrs)

    # Close the template declaration
    ss('> obj(m, "%(pyname)s"' % klassattrs)

    # Are we allowing dynamic attributes for the class?
    if klassattrs["dynamic_attr"]:
        ss(", py::dynamic_attr()")

    # Close the class declaration
    ss(");\n")

    # Is there a doc string?
    doc = inspect.getdoc(klass)
    if doc:
        ss("    obj.doc() = ")
        PYB11docstring(doc, ss)
        ss(";\n")

    # Grab all the methods
    allmethods = [(mname, meth) for (mname, meth) in PYB11ThisClassMethods(klass)
                  if not PYB11attrs(meth)["ignore"]]

    # Bind constructors of the class.
    ss("\n    // Constructors\n")
    kills = []
    for i, (mname, meth) in enumerate(allmethods):
        if mname[:6] == "pyinit":
            methattrs = PYB11attrs(meth)
            args = PYB11parseArgs(meth)
            pyinit(meth, methattrs, args)
            kills.append(i)
    for i in reversed(kills):
        del allmethods[i]

    # Bind special operators.
    ss("\n    // Operators\n")
    kills = []
    for i, (mname, meth) in enumerate(allmethods):
        methattrs = PYB11attrs(meth)
        methattrs["returnType"] = eval("klassinst." + mname + "()")
        args = PYB11parseArgs(meth)
        if methattrs["pyname"] in special_operators:
            func, op = special_operators[methattrs["pyname"]]
            func(meth, methattrs, args, op)
            kills.append(i)
    for i in reversed(kills):
        del allmethods[i]

    # Bind the remaining methods of the class.
    if allmethods:
        ss("\n    // Methods\n")
        for i, (mname, meth) in enumerate(allmethods):
            methattrs = PYB11attrs(meth)
            PYB11generic_class_method(klass, klassattrs, meth, methattrs, ss)

    # Bind attributes
    PYB11GenerateClassAttributes(klass, klassinst, klassattrs, ss)

    # Bind properties
    PYB11GenerateClassProperties(klass, klassinst, klassattrs, ss)

    # Bind any templated methods
    templates = [x for x in dir(klassinst) if isinstance(eval("klassinst.%s" % x), PYB11TemplateMethod) and x in klass.__dict__]
    if templates:
        ss("\n    // %(cppname)s template methods\n" % klassattrs)
        for tname in templates:
            inst = eval("klassinst.%s" % tname)
            inst(tname, klass, klassattrs, ss)

    # Helper method to check if the given method spec is already in allmethods
    def newOverloadedMethod(meth, allmethods, klassattrs):
        def extractArgs(mmeth):
            result = [x[0] for x in PYB11parseArgs(mmeth)]
            for i in xrange(len(result)):
                try:
                    result[i] = result[i] % klassattrs["template_dict"]
                except:
                    pass
            return result
        methattrs = PYB11attrs(meth)
        args = extractArgs(meth)
        overload = False
        for (othername, othermeth) in allmethods:
            othermethattrs = PYB11attrs(othermeth)
            if methattrs["cppname"] == othermethattrs["cppname"]:
                overload = True
                otherargs = extractArgs(othermeth)
                if otherargs == args:
                    return False
        return overload

    # If we're choosing to expose base hidden methods, check for those too.
    if klassattrs["exposeBaseOverloads"]:
        ss("\n    // Overloaded base methods\n")
        for bklass in inspect.getmro(klass)[1:]:
            bklassattrs = PYB11attrs(bklass)
            bcppname = "%(cppname)s" % bklassattrs
            if bklassattrs["template"]:
                bcppname += "<"
                for i, t in enumerate(bklassattrs["template"]):
                    if i < len(bklassattrs["template"]) - 1:
                        bcppname += ("%(" + t + ")s, ")
                    else:
                        bcppname += ("%(" + t + ")s>")
                bcppname = bcppname % klassattrs["template_dict"]
                bklassattrs["cppname"] = bcppname
            for mname, meth in PYB11ThisClassMethods(bklass):
                if ((not PYB11attrs(meth)["ignore"]) and                # Ignore the method?
                    (mname[:6] != "pyinit") and                         # Ignore constructors
                    newOverloadedMethod(meth, allmethods, klassattrs)): # New overload?
                    methattrs = PYB11attrs(meth)
                    PYB11generic_class_method(bklass, bklassattrs, meth, methattrs, ss)

            # Same thing with any base templated methods
            templates = [x for x in dir(bklass) if isinstance(eval("bklass.%s" % x), PYB11TemplateMethod) and x in bklass.__dict__]
            if templates:
                for tname in templates:
                    inst = eval("bklass.%s" % tname)
                    meth = inst.func_template
                    methattrs = PYB11attrs(meth)
                    if newOverloadedMethod(meth, allmethods, klassattrs):
                        try:
                            inst(tname, bklass, bklassattrs, ss)
                        except Exception as excpt:
                            #import sys
                            #sys.stderr.write("Encountered ERROR: %s\n" % excpt)
                            raise RuntimeError, "ERROR encountered processing (%s, %s) template instantiation\n    ERROR was: %s %s %s " % (tname, meth, excpt, type(excpt), excpt.args)

    # Look for any class scope enums and bind them
    enums = [x for x in dir(klassinst) if isinstance(eval("klassinst.%s" % x), PYB11enum) and x in klass.__dict__]
    if enums:
        ss("\n    // %(cppname)s enums\n  " % klassattrs)
        ssenum = PYB11indentedIO("  ")
        for ename in enums:
            inst = eval("klassinst.%s" % ename)
            inst(klass, ssenum, klassattrs)
        ss(ssenum.getvalue())
        ssenum.close()

    ss("  }\n\n")

    # Look for any class scope classes and bind them
    klasses = [(x, eval("klass.%s" % x)) for x in dir(klassinst) if (inspect.isclass(eval("klass.%s" % x)) and x in klass.__dict__)]
    klasses = sorted(klasses, key=PYB11sort_by_inheritance(klasses))
    for (kname, nklass) in klasses:
        #nklass = eval("klassinst.%s" % kname)
        nklassattrs = PYB11attrs(nklass)
        nklassattrs["pyname"] = klassattrs["pyname"] + "_" + nklassattrs["pyname"]
        nklassattrs["cppname"] = klassattrs["cppname"] + "::" + nklassattrs["cppname"]
        nklassattrs["template_dict"].update(klassattrs["template_dict"])
        PYB11generateClass(nklass, nklassattrs, ss)

    ssout(fs.getvalue() % klassattrs["template_dict"])
    fs.close()

    return
