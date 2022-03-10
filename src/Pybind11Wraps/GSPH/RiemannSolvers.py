from PYB11Generator import *
#from RiemannSolverBaseAbstractMethods import *
#-------------------------------------------------------------------------------
# Base class for riemann solver wave speeds
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
class RiemannSolverBase:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    """

    def pyinit(slopeLimiter = "LimiterBase<%(Dimension)s>&",
               waveSpeed = "WaveSpeedBase<%(Dimension)s>&",
               linearReconstruction = "const bool"):
        "slope limiter constructor"


    linearReconstruction = PYB11property("bool", "linearReconstruction", "linearReconstruction", doc="linear reconstruction of L and R state for riemann problem") 

    waveSpeed = PYB11property("WaveSpeedBase<%(Dimension)s>&", "waveSpeed",returnpolicy="reference_internal", doc="wave speed object")
    limiter = PYB11property("LimiterBase<%(Dimension)s>&", "limiter",returnpolicy="reference_internal", doc="slope limiter object")

    #DpDx = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DpDx",returnpolicy="reference_internal")
    #DvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "DvDx",returnpolicy="reference_internal")
    
#PYB11inject(RiemannSolverBaseAbstractMethods, RiemannSolverBase, pure_virtual=True)

#-------------------------------------------------------------------------------
# HLLC Approximate Riemann Solver
#-------------------------------------------------------------------------------
class HLLC(RiemannSolverBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    """

    def pyinit(slopeLimiter = "LimiterBase<%(Dimension)s>&",
               waveSpeed = "WaveSpeedBase<%(Dimension)s>&",
               linearReconstruction = "const bool"):
        "slope limiter constructor"

#-------------------------------------------------------------------------------
# HLLC Approximate Riemann Solver with constant grav acceleration
#-------------------------------------------------------------------------------
class GHLLC(HLLC):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    """

    def pyinit(slopeLimiter = "LimiterBase<%(Dimension)s>&",
               waveSpeed = "WaveSpeedBase<%(Dimension)s>&",
               linearReconstruction = "const bool",
               gravitationalAcceleration = "const Vector"):
        "slope limiter constructor"

    
    gravitationalAcceleration = PYB11property("Vector", "gravitationalAcceleration", "gravitationalAcceleration", doc="constant gravitational acceleration vector") 
