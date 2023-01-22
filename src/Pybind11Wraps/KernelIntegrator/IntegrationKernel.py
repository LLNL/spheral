from PYB11Generator import *

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class IntegrationKernel:
    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""
    def pyinit(self):
        "Evaluates functions for integration"

    @PYB11pure_virtual
    @PYB11const
    def extent(self,
               Hmult = "const Scalar"):
        "Get extent of kernel"
        return "double"
        
    @PYB11pure_virtual
    @PYB11const
    def evaluate(self,
                 xp = "const Vector&",
                 indices = "const std::vector<std::pair<int, int>>&",
                 position = "const FieldList<%(Dimension)s, Vector>&",
                 H = "const FieldList<%(Dimension)s, SymTensor>&",
                 volume = "const FieldList<%(Dimension)s, Scalar>&",
                 Hmult = "const Scalar",
                 values = "std::vector<Scalar>&",
                 dvalues = "std::vector<Vector>&"):
        "Evaluate some functions at a point"
        return "void"
