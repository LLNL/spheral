#-------------------------------------------------------------------------------
# PYB11function
#
# Handle functions in pybind11
#-------------------------------------------------------------------------------
from PYB11utils import *
from PYB11property import *
import copy, StringIO

#-------------------------------------------------------------------------------
# PYB11generateModuleFunctions
#
# Bind the methods in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleFunctions(modobj, ss):
    methods = PYB11functions(modobj)
    if methods:
        ss("  //...........................................................................\n")
        ss("  // Methods\n")
        for name, meth in methods:
            methattrs = PYB11attrs(meth)
            if not methattrs["ignore"]:
                PYB11generateFunction(meth, methattrs, ss)

    # Now look for any template function instantiations.
    func_templates = [x for x in dir(modobj) if isinstance(eval("modobj.%s" % x), PYB11TemplateFunction)]
    for ftname in func_templates:
        func_template = eval("modobj.%s" % ftname)
        func_template(ftname, ss)
    return

#--------------------------------------------------------------------------------
# Make a function template instantiation
#--------------------------------------------------------------------------------
class PYB11TemplateFunction:

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

    def __call__(self, pyname, ss):

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
        PYB11generateFunction(self.func_template, funcattrs, ss)
        if self.func_template.__doc__:
            self.func_template.__doc__ = doc0
        return

#-------------------------------------------------------------------------------
# Generate a function definition
#-------------------------------------------------------------------------------
def PYB11generateFunction(meth, methattrs, ssout):
    fs = StringIO.StringIO()
    ss = fs.write

    # Arguments
    stuff = inspect.getargspec(meth)
    nargs = len(stuff.args)
    argNames = stuff.args
    argTypes, argDefaults = [], []
    if nargs > 0:
        for thing in stuff.defaults:
            if isinstance(thing, tuple):
                assert len(thing) == 2
                argTypes.append(thing[0])
                argDefaults.append(thing[1])
            else:
                argTypes.append(thing)
                argDefaults.append(None)
    assert len(argNames) == nargs
    assert len(argTypes) == nargs
    assert len(argDefaults) == nargs

    # Return type
    returnType = meth(*tuple(stuff.args))
    methattrs["returnType"] = returnType

    # Build the argument string if any
    argString = ""
    for argType, argName, default in zip(argTypes, argNames, argDefaults):
        argString += (', "%s"_a' % argName)
        if methattrs["noconvert"]:
            argString += '.noconvert()'
        if not default is None:
            argString += ("=" + default)

    # Write the binding
    ss('  m.def("%(pyname)s", ' % methattrs)

    # If there is an implementation, short-circuit the rest of our work.
    if methattrs["implementation"]:
        ss(methattrs["implementation"] + argString)

    elif returnType:
        ss("(%s (*)(" % returnType)
        for i, argType in enumerate(argTypes):
            ss(argType)
            if i < nargs - 1:
                ss(", ")
        ss(")) &%(namespace)s%(cppname)s" % methattrs + argString)
    else:
        ss("&%(namespace)s%(cppname)s" % methattrs)

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

    ssout(fs.getvalue() % methattrs["template_dict"])
    fs.close()
    return

