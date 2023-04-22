from PYB11Generator import *
from SolidBoundaryAbstractMethods import *

#-------------------------------------------------------------------------------
# Base class for riemann solver wave speeds
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class SolidBoundary:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  """
    def pyinit():
        "constructor for base class DEM solid boundary conditions"

PYB11inject(SolidBoundaryAbstractMethods, SolidBoundary, pure_virtual=True)



#-------------------------------------------------------------------------------
# VanLeer limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class PlanarWall(SolidBoundary):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
  """

    def pyinit(point  = "const Vector&", 
               normal = "const Vector&"):
        "solid planar boundary"

    @PYB11virtual
    @PYB11const
    def velocity(self,
                 position = "const Vector&"):
        "velocity of bc."
        return "Vector"

    @PYB11virtual
    @PYB11const
    def distance(self,
                 position = "const Vector&"):
        "distance vector to bc."
        return "Vector"

    point  = PYB11property("const Vector&", "point",  "point", returnpolicy="reference_internal", doc="point  in plane definition")
    normal = PYB11property("const Vector&", "normal", "normal",returnpolicy="reference_internal", doc="normal in plane definition")
    