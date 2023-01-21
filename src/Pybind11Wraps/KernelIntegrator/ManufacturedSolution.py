from PYB11Generator import *
from IntegrationCoefficient import *

@PYB11template("Dimension")
@PYB11template_dict({"CoefficientType" : "double"})
@PYB11holder("std::shared_ptr")
class ManufacturedFunction(IntegrationCoefficient):
    PYB11typedefs = '''
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
'''

    def pyinit(self):
        "Manufactured function"

    @PYB11const
    @PYB11pure_virtual
    def evaluate(self,
                 t = "const double",
                 x = "const Vector&"):
        return "double"

    @PYB11const
    @PYB11pure_virtual
    def evaluateSpatialGradient(self,
                                t = "const double",
                                x = "const Vector&"):
        return "Vector"

    @PYB11const
    @PYB11pure_virtual
    def evaluateSpatialHessian(self,
                               t = "const double",
                               x = "const Vector&"):
        return "SymTensor"

    @PYB11const
    @PYB11pure_virtual
    def evaluateTimeDerivative(self,
                               t = "const double",
                               x = "const Vector&"):
        return "double"

    @PYB11const
    @PYB11virtual
    def evaluateCoefficient(self,
                            kid = "const KernelIntegrationData<%(Dimension)s>&"):
        return "double"

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class ManufacturedSteadyStateFunction(ManufacturedFunction):
    PYB11typedefs = '''
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
'''

    def pyinit(self,
               t = "const double",
               func = "std::shared_ptr<ManufacturedFunction<%(Dimension)s>>"):
        "Manufactured steady state function"
        
    @PYB11const
    @PYB11virtual
    def evaluate(self,
                 t = "const double",
                 x = "const Vector&"):
        return "double"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialGradient(self,
                                t = "const double",
                                x = "const Vector&"):
        return "Vector"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialHessian(self,
                               t = "const double",
                               x = "const Vector&"):
        return "SymTensor"

    @PYB11const
    @PYB11virtual
    def evaluateTimeDerivative(self,
                               t = "const double",
                               x = "const Vector&"):
        return "double"

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class ManufacturedConstantFunction(ManufacturedFunction):
    PYB11typedefs = '''
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
'''

    def pyinit(self,
               coefficient = "const double"):
        "Manufactured constant function"

    @PYB11const
    @PYB11virtual
    def evaluate(self,
                 t = "const double",
                 x = "const Vector&"):
        return "double"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialGradient(self,
                                t = "const double",
                                x = "const Vector&"):
        return "Vector"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialHessian(self,
                               t = "const double",
                               x = "const Vector&"):
        return "SymTensor"

    @PYB11const
    @PYB11virtual
    def evaluateTimeDerivative(self,
                               t = "const double",
                               x = "const Vector&"):
        return "double"

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class ManufacturedSinusoidalFunction(ManufacturedFunction):
    PYB11typedefs = '''
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
'''

    def pyinit(self,
               coefficients = "const std::vector<double>&"):
        "Manufactured sinusoidal function"

    @PYB11const
    @PYB11virtual
    def evaluate(self,
                 t = "const double",
                 x = "const Vector&"):
        return "double"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialGradient(self,
                                t = "const double",
                                x = "const Vector&"):
        return "Vector"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialHessian(self,
                               t = "const double",
                               x = "const Vector&"):
        return "SymTensor"

    @PYB11const
    @PYB11virtual
    def evaluateTimeDerivative(self,
                               t = "const double",
                               x = "const Vector&"):
        return "double"

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class ManufacturedWaveFunction(ManufacturedFunction):
    PYB11typedefs = '''
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::SymTensor SymTensor;
'''

    def pyinit(self,
               coefficients = "const std::vector<double>&"):
        "Manufactured wave function"

    @PYB11const
    @PYB11virtual
    def evaluate(self,
                 t = "const double",
                 x = "const Vector&"):
        return "double"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialGradient(self,
                                t = "const double",
                                x = "const Vector&"):
        return "Vector"

    @PYB11const
    @PYB11virtual
    def evaluateSpatialHessian(self,
                               t = "const double",
                               x = "const Vector&"):
        return "SymTensor"

    @PYB11const
    @PYB11virtual
    def evaluateTimeDerivative(self,
                               t = "const double",
                               x = "const Vector&"):
        return "double"
    
@PYB11template("Dimension")
@PYB11template_dict({"CoefficientType" : "std::vector<double>"})
@PYB11holder("std::shared_ptr")
class ManufacturedTransportSolution(IntegrationCoefficient):
    PYB11typedefs = '''
    typedef typename %(Dimension)s::Vector Vector;
'''

    def pyinit(self,
               c = "const double",
               numOrdinates = "const int",
               angularNorm = "const double",
               ordinates = "const std::vector<Vector>&",
               phiFunc = "std::shared_ptr<ManufacturedFunction<%(Dimension)s>>",
               sigmaAFunc = "std::shared_ptr<ManufacturedFunction<%(Dimension)s>>"):
        "Manufactured transport solution"

    @PYB11const
    @PYB11virtual
    def evaluatePhi(self,
                    t = "const double",
                    x = "const Vector&"):
        return "double"

    @PYB11const
    @PYB11virtual
    def evaluatePsi(self,
                    t = "const double",
                    x = "const Vector&"):
        return "std::vector<double>"


    @PYB11const
    @PYB11virtual
    def evaluateSource(self,
                       t = "const double",
                       x = "const Vector&"):
        return "std::vector<double>"


    @PYB11const
    @PYB11virtual
    def evaluateCoefficient(self,
                            t = "const double",
                            x = "const Vector&"):
        return "std::vector<double>"
    
    @PYB11const
    @PYB11virtual
    def evaluateCoefficient(self,
                            kid = "const KernelIntegrationData<%(Dimension)s>&"):
        return "std::vector<double>"
