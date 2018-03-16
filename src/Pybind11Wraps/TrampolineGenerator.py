#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
import inspect
import sys

class TrampolineGenerator:

    def __init__(self):
        self.includes = []
        self.namespaces = []
        self.templates = []
        self.preamble = None
        return

    def __call__(self):
        name = self.__class__.__name__
        ss = sys.stdout.write

        # Compiler guard.
        ss("""//------------------------------------------------------------------------------
// Trampoline class for %(name)s
//------------------------------------------------------------------------------
#ifndef __trampoline_%(name)s__
#define __trampoline_%(name)s__

""" % {"name" : name})
        # Includes
        for inc in self.includes:
            ss('#include "%s"\n' % inc)
        ss("\n")

        # Preamble
        if self.preamble:
            ss(self.preamble + "\n")

        # Namespaces
        for ns in self.namespaces:
            ss("namespace " + ns + " {\n")
        ss("\n")

        # Template parameters
        ss("template<")
        for tp in self.templates:
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
        if "Dimension" in self.templates:
            ss("""  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

""");

        # Bind the methods
        methods = [(name, meth) for (name, meth) in inspect.getmembers(self, predicate=inspect.ismethod)
                   if name[:2] != "__"]
        for name, method in methods:

            # Get the return type and arguments.
            stuff = inspect.getargspec(method)
            assert "args" in stuff.args
            returnType = method()
            args = stuff.defaults[stuff.args.index("args") - 1]
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
            for i, (argType, argName) in enumerate(args):
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

            for i, (argType, argName) in enumerate(args):
                if i < nargs - 1:
                    ss(offset + argName + ",\t// argument %i\n" % i)
                else:
                    ss(offset + argName + ");\t// argument %i\n" % (i + 1))
            ss("  }\n")

        # Closing
        ss("};\n\n")
        for ns in self.namespaces:
            ss("}\n")
        ss("\n#endif\n")

        return
