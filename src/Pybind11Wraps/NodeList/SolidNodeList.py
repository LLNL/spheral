from PYB11Generator import *
from FluidNodeList import FluidNodeList

#-------------------------------------------------------------------------------
# SolidNodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class SolidNodeList(FluidNodeList):
    "Spheral SolidNodeList base class in %(Dimension)s, i.e.,  the NodeList for solid dynamics."

    typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef Field<%(Dimension)s, int> IntField;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef Field<%(Dimension)s, Vector> VectorField;
    typedef Field<%(Dimension)s, Tensor> TensorField;
    typedef Field<%(Dimension)s, SymTensor> SymTensorField;
"""

    def pyinit(self,
               name = "std::string",
               eos = "EquationOfState<%(Dimension)s>&",
               strength = "StrengthModel<%(Dimension)s>&",
               numInternal = ("int", "0"),
               numGhost = ("int", "0"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("int", "500"),
               rhoMin = ("double", "1e-10"),
               rhoMax = ("double", "1e100")):
        "Constructor for a SolidNodeList class."
        return

    @PYB11virtual
    @PYB11const
    def soundSpeed(self, result="ScalarField&"):
        "Compute the current sound speed, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def bulkModulus(self, result="ScalarField&"):
        "Compute the current bulk modulus, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def yieldStrength(self, result="ScalarField&"):
        "Compute the current yield strength, storing the result in the argument ScalarField."
        return "void"

    @PYB11const
    def deviatoricStress(self):
        "The mass density field"
        return "const SymTensorField&"

    @PYB11const
    def plasticStrain(self):
        "The plastic strain field"
        return "const ScalarField&"

    @PYB11const
    def plasticStrainRate(self):
        "The plastic strain rate field"
        return "const ScalarField&"

    @PYB11const
    def damage(self):
        "The damage field"
        return "const SymTensorField&"

    @PYB11const
    def effectiveDamage(self):
        "The effective damage field"
        return "const SymTensorField&"

    @PYB11const
    def damageGradient(self):
        "The damage gradient field"
        return "const VectorField&"

    @PYB11const
    def fragmentIDs(self):
        "The fragment IDs field"
        return "const IntField&"

    @PYB11const
    def particleTypes(self):
        "The particle type field"
        return "const IntField&"

    @PYB11const
    def strengthModel(self):
        "Return the strength model object this SolidNodeList is associated with."
        return "const StrengthModel<%(Dimension)s>&"

    @PYB11virtual
    @PYB11const
    def label(self):
        "Label for restart files"
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self, file="FileIO&", pathName="const std::string&"):
        "Serialize under the given path in a FileIO object"
        return "void"

    @PYB11virtual
    def restoreState(self, file="const FileIO&", pathName="const std::string&"):
        "Restore state from the given path in a FileIO object"
        return "void"

    # Comparison
    def __eq__(self):
        "Equivalence test with another SolidNodeList"

    def __ne__(self):
        "Inequivalence test with another SolidNodeList"

