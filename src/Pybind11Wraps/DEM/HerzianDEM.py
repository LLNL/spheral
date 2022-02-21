#-------------------------------------------------------------------------------
# SPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DEMBase import *

@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class HerzianDEM(DEMBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename DEMBase<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(dataBase = "const DataBase<%(Dimension)s>&",
               YoungsModulus = "const Scalar",
               restitutionCoefficient = "const Scalar",
               stepsPerCollision = "const Scalar",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "DEMBase constructor"

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "calculate the derivatives for Linear Spring DEM."
        return "void"

    YoungsModulus = PYB11property("Scalar", "YoungsModulus", "YoungsModulus", doc="elastic modulus")
    restitutionCoefficient = PYB11property("Scalar", "restitutionCoefficient", "restitutionCoefficient", doc="linear restitution coefficient")
    beta = PYB11property("Scalar", "beta", "beta", doc="a damping parameter")
    
