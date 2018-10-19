#-------------------------------------------------------------------------------
# DataBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DataBase import *

@PYB11template("Dimension")
class DataBase:

    typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef typename DataBase<%(Dimension)s>::ConnectivityMapPtr ConnectivityMapPtr;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def reinitializeNeighbors(self):
        "Optimize all Neighbor objects for the current state."
        return "void"

    @PYB11const
    def updateConnectivityMap(self, computeGhostConnectivity="const bool"):
        "Update the internal connectivity map."
        return "void"

    @PYB11const
    def patchConnectivityMap(self,
                             flags = "const FieldList<%(Dimension)s, int>&",
                             old2new = "const FieldList<%(Dimension)s, int>&"):
        "Update the internal connectivity map."
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def connectivityMap(self, computeGhostConnectivity=("const bool", "false")):
        "Get the connectivity map, optionally including ghost connectivity"
        return "const ConnectivityMap<%(Dimension)s>&"

    @PYB11const
    def connectivityMapPtr(self, computeGhostConnectivity=("const bool", "false")):
        "Get the connectivity map as a std::shared_ptr, optionally including ghost connectivity"
        return "ConnectivityMapPtr"

    def appendNodeList(self, nodeList="SolidNodeList<%(Dimension)s>&"):
        "Add a SolidNodeList"
        return "void"

    @PYB11pycppname("appendNodeList")
    def appendNodeList1(self, nodeList="FluidNodeList<%(Dimension)s>&"):
        "Add a FluidNodeList"
        return "void"

    @PYB11pycppname("appendNodeList")
    def appendNodeList2(self, nodeList="NodeList<%(Dimension)s>&"):
        "Add a NodeList"
        return "void"

    def deleteNodeList(self, nodeList="SolidNodeList<%(Dimension)s>&"):
        "Remove a SolidNodeList"
        return "void"

    @PYB11pycppname("deleteNodeList")
    def deleteNodeList1(self, nodeList="FluidNodeList<%(Dimension)s>&"):
        "Remove a FluidNodeList"
        return "void"

    @PYB11pycppname("deleteNodeList")
    def deleteNodeList2(self, nodeList="NodeList<%(Dimension)s>&"):
        "Remove a NodeList"
        return "void"

    @PYB11const
    def haveNodeList(self, nodeList="NodeList<%(Dimension)s>&"):
        "Check if a NodeList is in the DataBase"
        return "bool"

    @PYB11const
    def setMasterNodeLists(self,
                           position = "const Vector&",
                           H = "const SymTensor&",
                           masterLists = "std::vector<std::vector<int>>&",
                           coarseNeighbors = "std::vector<std::vector<int>>&"):
        "Set the master/coarse neighbor lists for all NodeLists"
        return "void"

    @PYB11const
    def setMasterFluidNodeLists(self,
                                position = "const Vector&",
                                H = "const SymTensor&",
                                masterLists = "std::vector<std::vector<int>>&",
                                coarseNeighbors = "std::vector<std::vector<int>>&"):
        "Set the master/coarse neighbor lists for all NodeLists"
        return "void"

    @PYB11const
    def setRefineNodeLists(self,
                           position = "const Vector&",
                           H = "const SymTensor&",
                           coarseNeighbors = "const std::vector<std::vector<int>>&",
                           refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor lists for all NodeLists"
        return "void"

    @PYB11const
    def setRefineFluidNodeLists(self,
                                position = "const Vector&",
                                H = "const SymTensor&",
                                coarseNeighbors = "const std::vector<std::vector<int>>&",
                                refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor lists for all FluidNodeLists"
        return "void"

    #...........................................................................
    # Template methods
    @PYB11template("DataType")
    @PYB11const
    def newGlobalFieldList(self,
                           value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                           name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"')):
        "Construct a new FieldList<%(DataType)s> for all NodeLists in DataBase"
        return "FieldList<%(Dimension)s, %(DataType)s>"

    @PYB11template("DataType")
    @PYB11const
    def newFluidFieldList(self,
                          value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                          name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"')):
        "Construct a new FieldList<%(DataType)s> for all FluidNodeLists in DataBase"
        return "FieldList<%(Dimension)s, %(DataType)s>"

    @PYB11template("DataType")
    @PYB11const
    def newSolidFieldList(self,
                          value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                          name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"')):
        "Construct a new FieldList<%(DataType)s> for all SolidNodeLists in DataBase"
        return "FieldList<%(Dimension)s, %(DataType)s>"

    newGlobalIntFieldList              = PYB11TemplateMethod(newGlobalFieldList, template_parameters="int")
    newGlobalScalarFieldList           = PYB11TemplateMethod(newGlobalFieldList, template_parameters="double")
    newGlobalVectorFieldList           = PYB11TemplateMethod(newGlobalFieldList, template_parameters="Vector")
    newGlobalTensorFieldList           = PYB11TemplateMethod(newGlobalFieldList, template_parameters="Tensor")
    newGlobalSymTensorFieldList        = PYB11TemplateMethod(newGlobalFieldList, template_parameters="SymTensor")
    newGlobalThirdRankTensorFieldList  = PYB11TemplateMethod(newGlobalFieldList, template_parameters="ThirdRankTensor")
    newGlobalFourthRankTensorFieldList = PYB11TemplateMethod(newGlobalFieldList, template_parameters="FourthRankTensor")
    newGlobalFifthRankTensorFieldList  = PYB11TemplateMethod(newGlobalFieldList, template_parameters="FifthRankTensor")
    newGlobalvector_of_doubleFieldList = PYB11TemplateMethod(newGlobalFieldList, template_parameters="std::vector<double>")
    newGlobalvector_of_VectorFieldList = PYB11TemplateMethod(newGlobalFieldList, template_parameters="std::vector<Vector>")

    newFluidIntFieldList              = PYB11TemplateMethod(newFluidFieldList, template_parameters="int")
    newFluidScalarFieldList           = PYB11TemplateMethod(newFluidFieldList, template_parameters="double")
    newFluidVectorFieldList           = PYB11TemplateMethod(newFluidFieldList, template_parameters="Vector")
    newFluidTensorFieldList           = PYB11TemplateMethod(newFluidFieldList, template_parameters="Tensor")
    newFluidSymTensorFieldList        = PYB11TemplateMethod(newFluidFieldList, template_parameters="SymTensor")
    newFluidThirdRankTensorFieldList  = PYB11TemplateMethod(newFluidFieldList, template_parameters="ThirdRankTensor")
    newFluidFourthRankTensorFieldList = PYB11TemplateMethod(newFluidFieldList, template_parameters="FourthRankTensor")
    newFluidFifthRankTensorFieldList  = PYB11TemplateMethod(newFluidFieldList, template_parameters="FifthRankTensor")
    newFluidvector_of_doubleFieldList = PYB11TemplateMethod(newFluidFieldList, template_parameters="std::vector<double>")
    newFluidvector_of_VectorFieldList = PYB11TemplateMethod(newFluidFieldList, template_parameters="std::vector<Vector>")

    newSolidIntFieldList              = PYB11TemplateMethod(newSolidFieldList, template_parameters="int")
    newSolidScalarFieldList           = PYB11TemplateMethod(newSolidFieldList, template_parameters="double")
    newSolidVectorFieldList           = PYB11TemplateMethod(newSolidFieldList, template_parameters="Vector")
    newSolidTensorFieldList           = PYB11TemplateMethod(newSolidFieldList, template_parameters="Tensor")
    newSolidSymTensorFieldList        = PYB11TemplateMethod(newSolidFieldList, template_parameters="SymTensor")
    newSolidThirdRankTensorFieldList  = PYB11TemplateMethod(newSolidFieldList, template_parameters="ThirdRankTensor")
    newSolidFourthRankTensorFieldList = PYB11TemplateMethod(newSolidFieldList, template_parameters="FourthRankTensor")
    newSolidFifthRankTensorFieldList  = PYB11TemplateMethod(newSolidFieldList, template_parameters="FifthRankTensor")
    newSolidvector_of_doubleFieldList = PYB11TemplateMethod(newSolidFieldList, template_parameters="std::vector<double>")
    newSolidvector_of_VectorFieldList = PYB11TemplateMethod(newSolidFieldList, template_parameters="std::vector<Vector>")

    #...........................................................................
    # Properties
    numNodeLists = PYB11property("int", "numNodeLists", doc="Number of NodeLists in DataBase")
    numFluidNodeLists = PYB11property("int", "numFluidNodeLists", doc="Number of FluidNodeLists in DataBase")
    numSolidNodeLists = PYB11property("int", "numSolidNodeLists", doc="Number of SolidNodeLists in DataBase")

    numNodes = PYB11property("int", "numNodes", doc="Number of nodes in all NodeLists in DataBase")
    numInternalNodes = PYB11property("int", "numInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase")
    numGhostNodes = PYB11property("int", "numGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase")

    globalNumNodes = PYB11property("int", "globalNumNodes", doc="Number of nodes in all NodeLists in DataBase across all processors")
    globalNumInternalNodes = PYB11property("int", "globalNumInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase across all processors")
    globalNumGhostNodes = PYB11property("int", "globalNumGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase across all processors")

    nodeListPtrs = PYB11property("const std::vector<NodeList<%(Dimension)s>*>&", "nodeListPtrs", doc="The set of NodeLists in the DataBase")
    fluidNodeListPtrs = PYB11property("const std::vector<FluidNodeList<%(Dimension)s>*>&", "fluidNodeListPtrs", doc="The set of FluidNodeLists in the DataBase")
    solidNodeListPtrs = PYB11property("const std::vector<SolidNodeList<%(Dimension)s>*>&", "solidNodeListPtrs", doc="The set of SolidNodeLists in the DataBase")
