#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
import inspect
import sys

class TrampolineGenerator:

    def __init__(self):
        self.namespaces = []
        self.templates = []
        return

    def __call__(self):
        ss = sys.stdout.write

        # Namespaces
        for ns in self.namespaces:
            ss("namespace " + ns + " {")
        ss("\n\n")

        # Template parameters
        ss("template<typename Base")
        for tp in self.templates:
            ss(", typename " + tp)
        ss(">\n")

        # Class name
        ss("""
class %(name)s: public Base {
public:
  using Base::Base;   // inherit constructors

""" % {"name"     : self.__class__.__name__,
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
            stuff = inspect.getargspec(method)
            assert "returnType" in stuff.args
            assert "args" in stuff.args
            returnType = stuff.defaults[stuff.args.index("returnType") - 1]
            args = stuff.defaults[stuff.args.index("args") - 1]

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

            dvals = {"name" : name, "returnType" : returnType}
            firstline = "  virtual %(returnType)s %(name)s(" % dvals
            offset = " "*len(firstline)
            ss(firstline)
            for i, (argType, argName) in enumerate(args):
                if i > 0:
                    ss(offset)
                ss(argType + " " + argName)
                if i < len(args) - 1:
                    ss(",\n")
                else:
                    ss(")")
            if const:
                ss(" const override {\n")
            else:
                ss(" override {\n")

            if pure:
                ss("    PYBIND11_OVERLOAD_PURE(" + returnType + ",       // Return type\n")
            else:
                ss("    PYBIND11_OVERLOAD(" + returnType + ",       // Return type\n")

        # Closing
        ss("};\n\n")
        for ns in self.namespaces:
            ss("}\n")

        return
