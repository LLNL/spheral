from PYB11Generator import *

#-------------------------------------------------------------------------------
# NodeListRegistrar
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11singleton
class NodeListRegistrar:

    @PYB11const
    def valid(self):
        return "bool"

    # The instance attribute.  We expose this as a property of the class.
    @PYB11static
    @PYB11cppname("instance")
    @PYB11ignore
    def getinstance(self):
        return "NodeListRegistrar<%(Dimension)s>&"
    instance = property(getinstance, doc="The static NodeListRegistrar<%(Dimension)s> instance.")

