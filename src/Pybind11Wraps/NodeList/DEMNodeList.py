from PYB11Generator import *
from NodeList import NodeList
from RestartMethods import *

#-------------------------------------------------------------------------------
# FluidNodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralNodeList")
@PYB11dynamic_attr
class DEMNodeList(NodeList):
    "Spheral DEMNodeList base class in %(Dimension)s, i.e.,  the NodeList for Discrete Element Modelling."

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef Field<%(Dimension)s, int> intField;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
"""

    def pyinit(self,
               name = "std::string",
               numInternal = ("int", "0"),
               numGhost = ("int", "0"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               neighborSearchBuffer = ("Scalar","0.1"),
               maxNumNeighbors = ("int", "500")):
        "Constructor for a DEMNodeList class."
        return

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def particleRadius(self):
        "The particle radius field"
        return "const ScalarField&"


    @PYB11pycppname("particleRadius")
    def setparticleRadius(self, val="const ScalarField&"):
        "Set the particle radii"
        return "void"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def compositeParticleIndex(self):
        "the composite particle index field"
        return "const intField&"

    @PYB11pycppname("compositeParticleIndex")
    def setCompositeParticleIndex(self, val="const intField&"):
        "set the composite particle indices"
        return "void"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def uniqueIndex(self):
        "the unique particle index field"
        return "const intField&"

    @PYB11pycppname("uniqueIndex")
    def setUniqueIndex(self, val="const intField&"):
        "set the unique particle indices"
        return "void"

    def setHfieldFromParticleRadius(uniqueIndex = "const int"):
        "set a good H value for the neighbor search based on the particle radius"
        return "void"

    #...........................................................................
    # Comparison
    def __eq__(self):
        "Equivalence test with another DEMNodeList"

    def __ne__(self):
        "Inequivalence test with another DEMNodeList"

    #...........................................................................
    # Properties


#-------------------------------------------------------------------------------
# Inject the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, DEMNodeList)
