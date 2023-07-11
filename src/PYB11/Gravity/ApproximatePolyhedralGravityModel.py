#-------------------------------------------------------------------------------
# Approximate Polyhedral Gravity Model
#-------------------------------------------------------------------------------
from PYB11Generator import *

class ApproximatePolyhedralGravityModel:

    def pyinit(self,
               poly = "const GeomPolyhedron &",
               M = "const double",
               G = "const double"):
        """approximate polyhedral grav model constructor."""

    @PYB11const
    def potential(self,
                     position = "const Dim<3>::Vector&"):
        """Return the gravitation acceleration at a specified point"""
        return "Dim<3>::Scalar"

    @PYB11const
    def acceleration(self,
                     position = "const Dim<3>::Vector&"):
        """Return the gravitation acceleration at a specified point"""
        return "Dim<3>::Vector"

    numQuadraturePoints = PYB11property("unsigned int", "numQuadraturePoints", doc="number of quadrature points making up the model")
    quadraturePoints =    PYB11property("const std::vector<Dim<3>::Vector>&", returnpolicy="reference_internal")
    values =              PYB11property("const std::vector<Dim<3>::Vector>&", returnpolicy="reference_internal")
    resolutions =         PYB11property("const std::vector<Dim<3>::Scalar>&", returnpolicy="reference_internal")
