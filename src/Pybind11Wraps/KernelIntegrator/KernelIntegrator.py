from PYB11Generator import *

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class KernelIntegrator:
    def pyinit(self,
               integrationOrder = "const int",
               kernel = "const std::shared_ptr<IntegrationKernel<%(Dimension)s>>",
               dataBase = "const DataBase<%(Dimension)s>&",
               flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Performs integrals"

    @PYB11virtual
    def setState(time = "const double",
                 state = "const State<%(Dimension)s>&"):
        "Set the state, which should include Voronoi information"
        return "void"
    
    @PYB11virtual
    def addIntegral(self,
                    integral = "std::shared_ptr<KernelIntegralBase<%(Dimension)s>>"):
        "Add an integral to the integrator"
        return "void"

    @PYB11virtual
    def performIntegration(self):
        "Perform the integrals"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("getFlatConnectivity")
    def flatConnectivity(self):
        "Return flat connectivity"
        return "const FlatConnectivity<%(Dimension)s>&"
    
    @PYB11template("DataType")
    @PYB11const
    def coefficientsToValue(self,
                            coeffs = "const FieldList<%(Dimension)s, %(DataType)s>&",
                            value = "FieldList<%(Dimension)s, %(DataType)s>&"):
        "Given coefficients, return values at centers"
        return "void"

    @PYB11const
    def totalNumSubcells(self):
        "How many subcells did we create to integrate over?"
        return "int"

    @PYB11const
    def totalNumSubfacets(self):
        "How many subfacets did we create to integrate over?"
        return "int"
    
    coefficientsToValueScalar = PYB11TemplateMethod(coefficientsToValue, "typename %(Dimension)s::Scalar", pyname="coefficientsToValue")
    coefficientsToValueVector = PYB11TemplateMethod(coefficientsToValue, "typename %(Dimension)s::Vector", pyname="coefficientsToValue")
    coefficientsToValueTensor = PYB11TemplateMethod(coefficientsToValue, "typename %(Dimension)s::Tensor", pyname="coefficientsToValue")
    coefficientsToValueSymTensor = PYB11TemplateMethod(coefficientsToValue, "typename %(Dimension)s::SymTensor", pyname="coefficientsToValue")
