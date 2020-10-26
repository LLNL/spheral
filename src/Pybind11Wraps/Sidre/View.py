#-------------------------------------------------------------------------------
# Group
#-------------------------------------------------------------------------------
from PYB11Generator import *
# from axom import sidre
# from sidre import Group

@PYB11namespace("axom::sidre")
@PYB11singleton
class View:

    #...........................................................................
    # Constructors
    # def pyinit(self):
    #     "Default constructor"

    #...........................................................................
    # Methods

    # @PYB11namespace("axom::sidre::View")
    # def getData(self):
    #     "Return data held by view and cast it to any compatible type allowed by Conduit (return type depends on type caller assigns it to)."
    #     return "axom::sidre::Node::Value"

    #@PYB11namespace("axom::sidre::View")
    @PYB11implementation("""[](axom::sidre::View &self, int n) {
                                                                int* viewData = self.getData();
                                                                py::list result;
                                                                for (int i; i < n; i++)
                                                                     result.append(viewData[i]);
                                                                return result;
                                                               }""")
    def getDataA(self, n = "int"):
        "Return data held by view and cast it to any compatible type allowed by Conduit (return type depends on type caller assigns it to)."
        return "py::list"
