#-------------------------------------------------------------------------------
# Inject the virtual method bindings for Neighbor descendants.
#-------------------------------------------------------------------------------
from PYB11Generator import *
import inspect, types

@PYB11ignore
def injectNeighborVirtualMethods(cls):
    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList1(self,
                       position = "const Vector&",
                       H = "const Scalar&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList2(self,
                       position = "const Vector&",
                       H = "const SymTensor&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList1(self,
                               position = "const Vector&",
                               H = "const Scalar&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList2(self,
                               position = "const Vector&",
                               H = "const SymTensor&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given (position, H)"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList3(self,
                       position = "const Vector&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given position"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighobrList3(self,
                               position = "const Vector&",
                               coarseNeighbors = "const std::vector<int>&",
                               refineNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (coarse, refine) neighbor info for the given position"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList4(self,
                       enterPlane = "const Plane&",
                       exitPlane = "const Plane&",
                       masterList = "std::vector<int>&",
                       coarseNeighbors = "std::vector<int>&"):
        "Fill the given arrays with (master, coarse) neighbor info for the given (enter, exit) plane proximity"
        return "void"

    @PYB11virtual
    def updateNodes(self):
        "Update the internal connectivity information based on the state of associated NodeList"
        return "void"

    @PYB11virtual
    @PYB11pycppname("updateNodes")
    def updateNodes1(self,
                     nodeIDs = "const std::vector<int>&"):
        "Update the internal connectivity information for the given nodes based on the state of associated NodeList"
        return "void"

    @PYB11virtual
    @PYB11const
    def valid(self):
        "Test if the Neighbor is valid, i.e., ready to be queried for connectivity information."
        return "bool"

    # Inject em...
    names = [x for x in dir() if inspect.isfunction(eval(x))]
    for name in names:
        exec('cls.%(name)s = eval("%(name)s")' % {"name": name})
    return
