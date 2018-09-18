#-------------------------------------------------------------------------------
# PYB11TrampolineGenerator
#-------------------------------------------------------------------------------
import inspect
import sys

from PYB11utils import *

#-------------------------------------------------------------------------------
# PYB11generateTrampoline
#
# Generate the trampoline class, including pure virtual hooks.
#-------------------------------------------------------------------------------
def PYB11generateTrampoline(klass, klassattrs, ss):
    klassinst = klass()

    # Common preliminary code
    PYB11generateTrampolineClassStart(klass, klassattrs, ss)

    # Bind the virtual methods
    methods = [(mname, meth) for (mname, meth) in PYB11ClassMethods(klass)
               if not PYB11attrs(meth)["ignore"] and
               (PYB11attrs(meth)["virtual"] or PYB11attrs(meth)["pure_virtual"])]
    for mname, meth in methods:
        methattrs = PYB11attrs(meth)
        methattrs["returnType"] = eval("klassinst." + mname + "()")
        assert methattrs["returnType"]    # We require the full spec for virtual methods
        ss("  virtual %(returnType)s %(cppname)s(" % methattrs)

        # Fill out the argument list for this method
        args = PYB11parseArgs(meth)
        for i, (argType, argName, default) in enumerate(args):
            ss(argType)
            if i < len(args) - 1:
                ss(", ")
        if methattrs["const"]:
            ss(") const override { ")
        else:
            ss(") override { ")

        if methattrs["pure_virtual"]:
            ss("    PYBIND11_OVERLOAD_PURE(%s, " % methattrs["returnType"])
        else:
            ss("    PYBIND11_OVERLOAD(%s, " % methattrs["returnType"])
        ss(" %(namespace)s%(cppname)s," % klassattrs)
        if len(args) > 0:
            ss(" %(cppname)s, " % methattrs)
        else:
            ss(" %(cppname)s);" % methattrs)

        for i, (argType, argName, default) in enumerate(args):
            if i < nargs - 1:
                ss(argName + ", ")
            else:
                ss(argName + ");")
        ss(" }\n")

    # Closing
    PYB11generateTrampolineClassEnd(klass, klassattrs, ss)
    return

#-------------------------------------------------------------------------------
# generateBindingFunction
#
# Generate a function that provides pybind11 bindings for the virtual methods.
#-------------------------------------------------------------------------------
def PYB11generateBindingFunction(obj):
    ss = sys.stdout.write
    name = obj.__class__.__name__

    # Compiler guard.
    ss("""//------------------------------------------------------------------------------
// Pybind11 binding for virtual methods in %(name)s
//------------------------------------------------------------------------------
#ifndef __pybind11Bindings_%(name)s__
#define __pybing11Bindings_%(name)s__

""" % {"name" : name})

    # Includes
    for inc in obj.includes:
        ss('#include "%s"\n' % inc)
    ss("\n")

    # Preamble
    if obj.preamble:
        ss(obj.preamble + "\n")

    # Namespaces
    for ns in obj.namespaces:
        ss("namespace " + ns + " {\n")
    ss("\n")

    # Template parameters
    ss("template<")
    for tp in obj.templates:
        ss("typename %s, " % tp)
    ss("typename Obj, typename PB11Obj>\n")

    # Function spec
    ss("void virtual%sBindings(PB11Obj& obj) {\n\n" % name)

    # typedefs
    if "Dimension" in obj.templates:
        ss("""
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

""");

    # Bind methods.
    ss("  // Methods\n")
    methods = [(name, meth) for (name, meth) in inspect.getmembers(obj, predicate=inspect.ismethod)
               if name[:2] != "__"]
    for name, method in methods:

        # Get the return type and arguments.
        returnType = method()
        stuff = inspect.getargspec(method)
        if "args" in stuff.args:
            args = stuff.defaults[stuff.args.index("args") - 1]
        else:
            args = []
        nargs = len(args)

        # Is this method const?
        if "const" in stuff.args:
            const = stuff.defaults[stuff.args.index("const") - 1]
        else:
            const = False

        # Is there a doc string?
        if "doc" in stuff.args:
            doc = stuff.defaults[stuff.args.index("doc") - 1]
        else:
            doc = None

        # Because python does not have function overloading, we provide the ability
        # to rename the c++ method.
        if "name" in stuff.args:
            name = stuff.defaults[stuff.args.index("name") - 1]

        # Write the binding
        dvals = {"name" : name, "returnType" : returnType}
        ss('  obj.def("%(name)s", (%(returnType)s (Obj::*)(' % dvals)
        for i, (argType, argName, default) in enumerate(__parseArgs(args)):
            ss(argType)
            if i < nargs - 1:
                ss(", ")
        if const:
            ss(") const)")
        else:
            ss("))")
        ss(" &Obj::" + name)
        for argType, argName, default in __parseArgs(args):
            ss(', "%s"_a' % argName)
            if default:
                ss("=" + default)
        if doc:
            ss(',\n          "%s"' % doc)
        ss(");\n")

    # Closing
    ss("}\n\n")
    for ns in obj.namespaces:
        ss("}\n")
    ss("#endif\n")
    
    return

#-------------------------------------------------------------------------------
# PYB11generateTrampolineClassStart
#
# All the stuff up to the methods.
#-------------------------------------------------------------------------------
def PYB11generateTrampolineClassStart(klass, klassattrs, ss):

    # Compiler guard.
    ss("""//------------------------------------------------------------------------------
// Trampoline class for %(cppname)s
//------------------------------------------------------------------------------
#ifndef __trampoline_%(pyname)s__
#define __trampoline_%(pyname)s__

""" % klassattrs)

    # Namespaces
    for ns in klassattrs["namespace"].split("::")[:-1]:
        ss("namespace " + ns + " {\n")

    # Class name
    ss("""
class PYB11Trampoline%(cppname)s: public %(cppname)s {
public:
  using %(cppname)s::%(cppname)s;   // inherit constructors

 """ % klassattrs)

    return

#-------------------------------------------------------------------------------
# PYB11generateTrampolineClassEnd
#
# Finish up the code.
#-------------------------------------------------------------------------------
def PYB11generateTrampolineClassEnd(klass, klassattrs, ss):
    ss("};\n\n")
    for ns in klassattrs["namespace"].split("::")[:-1]:
        ss("}\n")
    ss("\n#endif\n")
    return
