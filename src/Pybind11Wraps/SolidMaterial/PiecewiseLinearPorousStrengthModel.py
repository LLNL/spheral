#-------------------------------------------------------------------------------
# PorousStrengthModel
#-------------------------------------------------------------------------------
from PYB11Generator import *
from PorousStrengthModel import *

@PYB11template("Dimension")
@PYB11module("SpheralSolidMaterial")
class PiecewiseLinearPorousStrengthModel(PorousStrengthModel):
    """PorousStrengthModel

    The shear modulus and yield strength ratios are treated as piecewise
    linear functions of the porosity. Here ratio means the porous value
    divided by the value from the solidStrength model."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               solidStrength = "const StrengthModel<%(Dimension)s>&",
               porosityAbscissa = "const vector<Scalar>",
               shearModulusRatios = "const vector<Scalar>",
               yieldStrengthRatios = "const vector<Scalar>",):
        "Construct with the strength model we're modifying"


    def getPorousShearModulusRatio(self,
                                   alphai = "const Scalar"):
        "returns the ratio porous shear modulus to solid shear modulus"
        return "Scalar"

    def getPorousYieldStrengthRatio(self,
                                    alphai = "const Scalar"):
        "returns the ratio porous shear modulus to solid shear modulus"
        return "Scalar"
    #...........................................................................
    # Properties
    porosityAbscissa = PYB11property("vector<Scalar>&","porosityAbscissa", returnpolicy="reference_internal")
    shearModulusRatios = PYB11property("vector<Scalar>&","shearModulusRatios", returnpolicy="reference_internal")
    yieldStrengthRatios = PYB11property("vector<Scalar>&","yieldStrengthRatios", returnpolicy="reference_internal")



