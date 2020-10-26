#-------------------------------------------------------------------------------
# Group
#-------------------------------------------------------------------------------
from PYB11Generator import *
# from axom import sidre
# from sidre import Group

@PYB11namespace("axom::sidre")
@PYB11singleton
class Group:

    #...........................................................................
    # Constructors
    # def pyinit(self):
    #     "Default constructor"

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11pycppname("print")
    def printGroup(self):
        "JSON description of data Group to stdout."
        return "void"

    @PYB11namespace("axom::sidre::Group")
    def getView(path = "const std::string&"):
        "Return pointer to non-const View with given name or path."
        return "axom::sidre::View*"

    # @PYB11namespace("axom::sidre::View")
    # def getData():
    #     "Return data held by view and cast it to any compatible type allowed by Conduit (return type depends on type caller assigns it to)."
    #     return "axom::sidre::Node::Value"

    #...........................................................................
    # Attributes
