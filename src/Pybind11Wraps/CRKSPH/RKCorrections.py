#-------------------------------------------------------------------------------
# RKCorrections
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class RKCorrections(Physics):
    "Computes RK correction terms"
    
    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    def pyinit(self,
               dataBase = "const DataBase<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               correctionOrder = "const CRKOrder",
               volumeType = "const CRKVolumeType"):
        "Constructor"
        
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Evaluate derivatives"
        return "void"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state"
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register derivatives"
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields"
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields"
        return "void"
    
    @PYB11virtual
    def initializeProblemStartup(self, dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup"
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step"
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "void"
                  
    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&",
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the hydro at the completion of an integration step."
        return "void"
                  
    #...........................................................................
    # Properties
    correctionOrder = PYB11property("CRKOrder", "correctionOrder", "correctionOrder",
                                    doc="Flag to choose CRK Correction Order")
    volumeType = PYB11property("CRKVolumeType", "volumeType", "volumeType",
                               doc="Flag for the CRK volume weighting definition")
    volume = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "volume", returnpolicy="reference_internal")

    A = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "A", returnpolicy="reference_internal")
    B = PYB11property("const FieldList<%(Dimension)s, Vector>&", "B", returnpolicy="reference_internal")
    C = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "C", returnpolicy="reference_internal")
    D = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "D", returnpolicy="reference_internal")
    gradA = PYB11property("const FieldList<%(Dimension)s, Vector>&", "gradA", returnpolicy="reference_internal")
    gradB = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "gradB", returnpolicy="reference_internal")
    gradC = PYB11property("const FieldList<%(Dimension)s, ThirdRankTensor>&", "gradC", returnpolicy="reference_internal")
    gradD = PYB11property("const FieldList<%(Dimension)s, FourthRankTensor>&", "gradD", returnpolicy="reference_internal")
    hessA = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "hessA", returnpolicy="reference_internal")
    hessB = PYB11property("const FieldList<%(Dimension)s, ThirdRankTensor>&", "hessB", returnpolicy="reference_internal")
    hessC = PYB11property("const FieldList<%(Dimension)s, FourthRankTensor>&", "hessC", returnpolicy="reference_internal")
    hessD = PYB11property("const FieldList<%(Dimension)s, FifthRankTensor>&", "hessD", returnpolicy="reference_internal")
    
    surfacePoint = PYB11property("const FieldList<%(Dimension)s, int>&", "surfacePoint", returnpolicy="reference_internal")
    etaVoidPoints = PYB11property("const FieldList<%(Dimension)s, std::vector<Vector>>&", "etaVoidPoints", returnpolicy="reference_internal")
    cells = PYB11property("const FieldList<%(Dimension)s, FacetedVolume>&", "cells", returnpolicy="reference_internal")
    cellFaceFlags = PYB11property("const FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&", "cellFaceFlags", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, RKCorrections, virtual=True)
PYB11inject(RestartMethods, RKCorrections)

