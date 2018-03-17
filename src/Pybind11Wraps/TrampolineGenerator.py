#-------------------------------------------------------------------------------
# TrampolineGenerator
#-------------------------------------------------------------------------------
import inspect
import sys

class TrampolineGenerator:
    """Users should inherit from this base class and provide

Optional attributes to override:
  self.includes : a list of string include paths
"""

    def __init__(self):
        self.includes = []
        self.namespaces = []
        self.templates = []
        self.preamble = None
        return

#-------------------------------------------------------------------------------
# generateAbstractTrampoline
#
# Generate the abstract trampoline class
#-------------------------------------------------------------------------------
def generateAbstractTrampoline(obj):
    ss = sys.stdout.write
    name = obj.__class__.__name__ + "Abstract"
    __generateClassStart(obj, ss, name)

    # Bind the methods
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

        # Is this method abstract?
        if "pure" in stuff.args:
            pure = stuff.defaults[stuff.args.index("pure") - 1]
        else:
            pure = False

        # Generate the spec for this method
        dvals = {"name" : name, "returnType" : returnType}
        firstline = "  virtual %(returnType)s %(name)s(" % dvals
        offset = " "*len(firstline)
        ss(firstline)
        for i, (argType, argName, default) in enumerate(__parseArgs(args)):
            if i > 0:
                ss(offset)
            ss(argType + " " + argName)
            if i < nargs - 1:
                ss(",\n")
            else:
                ss(")")
        if const:
            ss(" const override {\n")
        else:
            ss(" override {\n")

        if pure:
            ss("    PYBIND11_OVERLOAD_PURE(" + returnType + ",\t// Return type\n")
            offset = "                           "
        else:
            ss("    PYBIND11_OVERLOAD(" + returnType + ",\t// Return type\n")
            offset = "                      "
        ss(offset + "Base,\t\t// Parent class\n")
        ss(offset + name + ",\t// name of method\n")

        for i, (argType, argName, default) in enumerate(__parseArgs(args)):
            if i < nargs - 1:
                ss(offset + argName + ",\t// argument %i\n" % (i + 1))
            else:
                ss(offset + argName + ");\t// argument %i\n" % (i + 1))
        ss("  }\n")

    # Closing
    __generateClassEnd(obj, ss)
    return

#-------------------------------------------------------------------------------
# generateConcreteTrampoline
#
# Overload just the abstract methods with non-abstract overrides.
#-------------------------------------------------------------------------------
def generateConcreteTrampoline(obj):
    ss = sys.stdout.write
    name = obj.__class__.__name__
    __generateClassStart(obj, ss, name)

    # Bind the methods
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

        # Is this method abstract?
        if "pure" in stuff.args:
            pure = stuff.defaults[stuff.args.index("pure") - 1]
        else:
            pure = False

        if pure:
            # Generate the spec for this method
            dvals = {"name" : name, "returnType" : returnType}
            firstline = "  virtual %(returnType)s %(name)s(" % dvals
            offset = " "*len(firstline)
            ss(firstline)
            for i, (argType, argName, default) in enumerate(__parseArgs(args)):
                if i > 0:
                    ss(offset)
                ss(argType + " " + argName)
                if i < nargs - 1:
                    ss(",\n")
                else:
                    ss(")")
            if const:
                ss(" const override {\n")
            else:
                ss(" override {\n")

            ss("    PYBIND11_OVERLOAD(" + returnType + ",\t// Return type\n")
            offset = "                      "
            ss(offset + "Base,\t\t// Parent class\n")
            ss(offset + name + ",\t// name of method\n")
        
            for i, (argType, argName, default) in enumerate(__parseArgs(args)):
                if i < nargs - 1:
                    ss(offset + argName + ",\t// argument %i\n" % (i + 1))
                else:
                    ss(offset + argName + ");\t// argument %i\n" % (i + 1))
            ss("  }\n")
        
    # Closing
    __generateClassEnd(obj, ss)

    return

#-------------------------------------------------------------------------------
# generateBindingFunction
#
# Generate a function that provides pybind11 bindings for the virtual methods.
#-------------------------------------------------------------------------------
def generateBindingFunction(obj):
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
  typedef typename Dimension::Tensor Tensor
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
# __generateClassStart
#
# All the stuff up to the methods.
#-------------------------------------------------------------------------------
def __generateClassStart(obj, ss, name):

    # Compiler guard.
    ss("""//------------------------------------------------------------------------------
// Trampoline class for %(name)s
//------------------------------------------------------------------------------
#ifndef __trampoline_%(name)s__
#define __trampoline_%(name)s__

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
    ss("typename Base>\n")

    # Class name
    ss("""
class %(name)s: public Base {
public:
  using Base::Base;   // inherit constructors

""" % {"name"     : name,
      })

    # typedefs
    if "Dimension" in obj.templates:
        ss("""
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

""");

    return

#-------------------------------------------------------------------------------
# __generateClassEnd
#
# Finish up the code.
#-------------------------------------------------------------------------------
def __generateClassEnd(obj, ss):
    ss("};\n\n")
    for ns in obj.namespaces:
        ss("}\n")
    ss("\n#endif\n")
    return

#-------------------------------------------------------------------------------
# __parseArgs
#
# Return (argType, argName, default_value (optional)
#-------------------------------------------------------------------------------
def __parseArgs(args):
    result = []
    for tup in args:
        if len(tup) == 2:
            result.append((tup[0], tup[1], None))
        else:
            result.append(tup)
    return result
