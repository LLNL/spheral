#-------------------------------------------------------------------------------
# SPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DEMBase import *

@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class LinearSpringDEM(DEMBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename DEMBase<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(dataBase = "const DataBase<%(Dimension)s>&",
               normalSpringConstant = "const Scalar",
               normalRestitutionCoefficient = "const Scalar",
               tangentialSpringConstant = "const Scalar",
               tangentialRestitutionCoefficient = "const Scalar",
               dynamicFrictionCoefficient = "const Scalar",
               staticFrictionCoefficient = "const Scalar",
               rollingFrictionCoefficient = "const Scalar",
               torsionalFrictionCoefficient = "const Scalar",
               cohesiveTensileStrength = "const Scalar",
               shapeFactor = "const Scalar",
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

    normalSpringConstant = PYB11property("Scalar", "normalSpringConstant", "normalSpringConstant", doc="normal spring constant")
    normalRestitutionCoefficient = PYB11property("Scalar", "normalRestitutionCoefficient", "normalRestitutionCoefficient", doc="normal restitution coefficient")
    tangentialSpringConstant = PYB11property("Scalar", "tangentialSpringConstant", "tangentialSpringConstant", doc="tangential spring constant")
    tangentialRestitutionCoefficient = PYB11property("Scalar", "tangentialRestitutionCoefficient", "tangentialRestitutionCoefficient", doc="tangential restitution coefficient")
    cohesiveTensileStrength = PYB11property("Scalar", "cohesiveTensileStrength", "cohesiveTensileStrength", doc="constant for normal cohesion constant")
    
    dynamicFrictionCoefficient = PYB11property("Scalar", "dynamicFrictionCoefficient", "dynamicFrictionCoefficient", doc="sliding friction coefficient - dynamic")
    staticFrictionCoefficient = PYB11property("Scalar", "staticFrictionCoefficient", "staticFrictionCoefficient", doc="sliding friction coefficient - static")
    rollingFrictionCoefficient = PYB11property("Scalar", "rollingFrictionCoefficient", "rollingFrictionCoefficient", doc="rolling friction coefficient")
    torsionalFrictionCoefficient = PYB11property("Scalar", "torsionalFrictionCoefficient", "torsionalFrictionCoefficient", doc="torsional friction coefficient")
    
    shapeFactor = PYB11property("Scalar", "shapeFactor", "shapeFactor", doc="shape factor - simple approach to non-spherical particles")
    normalBeta = PYB11property("Scalar", "normalBeta", "normalBeta", doc="a damping parameter")
    tangentialBeta = PYB11property("Scalar", "tangentialBeta", "tangentialBeta", doc="a damping parameter")
    timeStep = PYB11property("Scalar", "timeStep", "timeStep", doc="constant time-step for this model")