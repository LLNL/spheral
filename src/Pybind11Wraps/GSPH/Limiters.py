from PYB11Generator import *
from SlopeLimiterBaseAbstractMethods import *

#-------------------------------------------------------------------------------
# Base class for riemann solver wave speeds
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
class SlopeLimiterBase:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "slope limiter constructor"

    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

PYB11inject(SlopeLimiterBaseAbstractMethods, SlopeLimiterBase, pure_virtual=True)

#-------------------------------------------------------------------------------
# minmod limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
class MinModLimiter:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "minmod slope limiter constructor"

    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# VanLeer limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
class VanLeerLimiter:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "VanLeer slope limiter constructor"

    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# VanAlba limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
class VanAlbaLimiter:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "VanAlba slope limiter constructor"

    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# Superbee limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
class SuperbeeLimiter:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "Superbee slope limiter constructor"

    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# Ospre limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralGSPH")
class OspreLimiter:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "Ospre slope limiter constructor"

    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"



