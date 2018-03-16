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

        # pure methods


        # Closing
        ss("};\n\n")
        for ns in self.namespaces:
            ss("}\n")

        return
