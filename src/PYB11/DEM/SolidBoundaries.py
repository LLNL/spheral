from PYB11Generator import *
from SolidBoundaryBaseAbstractMethods import *

#-------------------------------------------------------------------------------
# Base class for riemann solver wave speeds
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class SolidBoundaryBase:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  """
    def pyinit():
        "constructor for base class DEM solid boundary conditions"

PYB11inject(SolidBoundaryBaseAbstractMethods, SolidBoundaryBase, pure_virtual=True)

#-------------------------------------------------------------------------------
# Infinite planar solid boundary
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class InfinitePlaneSolidBoundary(SolidBoundaryBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
  """

    def pyinit(point  = "const Vector&", 
               normal = "const Vector&"):
        "solid planar boundary"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state solid bc expects to use and evolve."
        return "void"
    
    @PYB11virtual
    def update(self,
               multiplier = "const double",
               t = "const double",
               dt = "const double",):
        "distance vector to bc."
        return "void"

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

    velocity = PYB11property("const Vector&", "velocity",  "velocity", returnpolicy="reference_internal", doc="velocity of plane")
    point  = PYB11property("const Vector&", "point",  "point", returnpolicy="reference_internal", doc="point  in plane definition")
    normal = PYB11property("const Vector&", "normal", "normal",returnpolicy="reference_internal", doc="normal in plane definition")
    
#-------------------------------------------------------------------------------
# Finite rectangular solid boundary
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class RectangularPlaneSolidBoundary(SolidBoundaryBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
  """

    def pyinit(point  = "const Vector&",
               extent = "const Vector&",
               basis  = "const Tensor&"):
        "solid planar boundary"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state solid bc expects to use and evolve."
        return "void"
    
    @PYB11virtual
    def update(self,
               multiplier = "const double",
               t = "const double",
               dt = "const double",):
        "distance vector to bc."
        return "void"

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

    velocity = PYB11property("const Vector&", "velocity",  "velocity", returnpolicy="reference_internal", doc="velocity of plane")
    point  = PYB11property("const Vector&", "point",  "point", returnpolicy="reference_internal", doc="point in plane definition")
    extent  = PYB11property("const Vector&", "extent",  "extent", returnpolicy="reference_internal", doc="extent of rectangle")
    basis = PYB11property("const Tensor&", "basis", "basis",returnpolicy="reference_internal", doc="basis vectors for rectangle")
    
#-------------------------------------------------------------------------------
# Disk shaped solid boundary
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class CircularPlaneSolidBoundary(SolidBoundaryBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
  """

    def pyinit(point  = "const Vector&",
               normal  = "const Vector&",
               extent = "const Scalar"):
        "solid planar boundary"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state solid bc expects to use and evolve."
        return "void"
    
    @PYB11virtual
    def update(self,
               multiplier = "const double",
               t = "const double",
               dt = "const double",):
        "distance vector to bc."
        return "void"

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

    velocity = PYB11property("const Vector&", "velocity",  "velocity", returnpolicy="reference_internal", doc="velocity of plane")
    point  = PYB11property("const Vector&", "point",  "point", returnpolicy="reference_internal", doc="point in plane definition")
    extent  = PYB11property("Scalar", "extent",  "extent", returnpolicy="reference_internal", doc="extent of rectangle")
    normal = PYB11property("const Vector&", "normal", "normal",returnpolicy="reference_internal", doc="normal in plane definition")

#-------------------------------------------------------------------------------
# Cylinder solid boundary. In 2d this would be two planes.
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class CylinderSolidBoundary(SolidBoundaryBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
  """

    def pyinit(point  = "const Vector&",
               axis  = "const Vector&",
               radius = "const Scalar",
               length = "const Scalar"):
        "solid planar boundary"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state solid bc expects to use and evolve."
        return "void"
    
    @PYB11virtual
    def update(self,
               multiplier = "const double",
               t = "const double",
               dt = "const double",):
        "distance vector to bc."
        return "void"

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

    velocity = PYB11property("const Vector&", "velocity",  "velocity", returnpolicy="reference_internal", doc="velocity of plane")
    point  = PYB11property("const Vector&", "point",  "point", returnpolicy="reference_internal", doc="point in plane definition")
    axis  = PYB11property("const Vector&", "axis",  "axis", returnpolicy="reference_internal", doc="extent of rectangle")
    radius = PYB11property("Scalar", "radius", "radius", doc="normal in plane definition")
    length = PYB11property("Scalar", "length", "length", doc="normal in plane definition")


#-------------------------------------------------------------------------------
# Cylinder solid boundary. In 2d this would be two planes.
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class SphereSolidBoundary(SolidBoundaryBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
  """

    def pyinit(center  = "const Vector&",
               radius = "const Scalar",
               clipPoint  = "const Vector&",
               clipAxis = "const Vector&"):
        "solid planar boundary"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state solid bc expects to use and evolve."
        return "void"
    
    @PYB11virtual
    def update(self,
               multiplier = "const double",
               t = "const double",
               dt = "const double",):
        "distance vector to bc."
        return "void"

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

    velocity = PYB11property("const Vector&", "velocity",  "velocity", returnpolicy="reference_internal", doc="velocity of plane")
    center  = PYB11property("const Vector&", "center",  "center", returnpolicy="reference_internal", doc="center of sphere")
    radius = PYB11property("Scalar", "radius", "radius", doc="radius of sphere")
    clipPoint  = PYB11property("const Vector&", "clipPoint",  "clipPoint", returnpolicy="reference_internal", doc="point on clip plane")
    clipAxis = PYB11property("const Vector&", "clipAxis", "clipAxis", returnpolicy="reference_internal", doc="normal in clip plane")