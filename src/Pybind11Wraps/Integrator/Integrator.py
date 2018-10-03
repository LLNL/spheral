#-------------------------------------------------------------------------------
# Integrator base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class Integrator:
    "Base class for all Spheral time integration algorithms"

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Construct an itegrator"

    def pyinit1(self, dataBase = "DataBase<DIM>&"):
        "Construct an integrator with a DataBase"

    def pyinit1(self,
                dataBase = "DataBase<DIM>&",
                physicsPackages = "const std::vector<Physics<DIM>*>&"):
        "Construct an integrator with a DataBase and physics packages"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def step(self, maxTime="Scalar"):
        "Take a step"
        return "void"

    @PYB11virtual
    @PYB11const
    def selectDt(self,
                 dtMin = "const Scalar",
                 dtMax = "const Scalar",
                 state = "const State<DIM>&",
                 derivs = "const StateDerivatives<DIM>&"):
        "Provide a method of looping over the physics packages and picking a time step."
        return "Scalar"

    @PYB11virtual
    def preStepInitialize(self,
                          state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        """Perform generic initializations at the beginning of a timestep.
To be called once per advance cycle."""
        return "void"

    @PYB11virtual
    def initializeDerivatives(self,
                              t = "const double",
                              dt = "const double",
                              state = "State<DIM>&",
                              derivs = "StateDerivatives<DIM>&"):
        """Prepare all physics packages for calls to evaluateDerivatives.
To be called before any call to Physics::evaluateDerivatives, therefore potentially
several times during a time step."""
        return "void"

    @PYB11virtual
    def postStepFinalize(self,
                         t = "const double",
                         dt = "const double",
                         state = "State<DIM>&",
                         derivs = "StateDerivatives<DIM>&"):
        "Finalize at the end a timestep, therefore called once at the end of a timestep."
        return "void"

    @PYB11virtual
    def advance(self, goalTime="Scalar"):
        "Advance the set of Physics packages to the given time."
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def evaluateDerivatives(self,
                            t = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Iterate over all physics packages and call evaluateDerivatives."
        return "void"

    @PYB11const
    def finalizeDerivatives(self,
                            t = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Iterate over all physics packages and call finalizeDerivatives."
        return "void"

    @PYB11const
    def postStateUpdate(self,
                        t = "const Scalar",
                        dt = "const Scalar",
                        dataBase = "const DataBase<DIM>&",
                        state = "State<DIM>&",
                        derivs = "StateDerivatives<DIM>&"):
        "Iterate over all physics packages and call postStateUpdate"
        return "void"

    def appendPhysicsPackage(self, package="Physics<DIM>&"):
        "Add a Physics package."
        return "void"

    @PYB11const
    def havePhysicsPackage(self, package="const Physics<DIM>&"):
        "Test if the given Physics package is listed in the integrator."
        return "bool"

    @PYB11const
    def uniqueBoundaryConditions(self):
        "Get the unique set of boundary conditions across all physics packages."
        return "std::vector<Boundary<DIM>*>"

    def setGhostNodes(self):
        "Set the ghost nodes for all node lists according to the boundary conditions."
        return "void"

    def applyGhostBoundaries(self,
                             state = "State<DIM>&",
                             derivs = "StateDerivatives<DIM>&"):
        "Set the ghost node values on the Fields of the nodes lists in the data base."
        return "void"

    def finalizeGhostBoundaries(self):
        "Finalize the ghost node boundary conditions."
        return "void"

    def setViolationNodes(self):
        "Find the nodes in violation of the boundary conditions."
        return "void"

    def enforceBoundaries(self,
                          state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        """Reset any internal nodes in violation of boundary conditions to be brought 
into compliance."""
        return "void"

    @PYB11const
    def copyGhostState(self,
                       state0 = "const State<DIM>&",
                       state1 = "State<DIM>&"):
        "Copy the ghost positions and H's from one state to another."
        return "void"

    #...........................................................................
    # Properties
    currentTime = PYB11property("Scalar", "currentTime", "currentTime", doc="Simulation time")
    currentCycle = PYB11property("int", "currentCycle", "currentCycle", doc="Simulation cycle")
    dtMin = PYB11property("Scalar", "dtMin", "dtMin", doc="Minimum allowed time step")
    dtMax = PYB11property("Scalar", "dtMax", "dtMax", doc="Maximum allowed time step")
    lastDt = PYB11property("Scalar", "lastDt", "lastDt", doc="Last timestep used")
    dtGrowth = PYB11property("Scalar", "dtGrowth", "dtGrowth", doc="Maximum allowed fractional time step growth")
    dataBase = PYB11property("DataBase<DIM>&", "dataBase", returnpolicy="reference_internal", doc="The DataBase of NodeLists")
    physicsPackages = PYB11property("const std::vector<Physics<DIM>*>&", returnpolicy="reference_internal", doc="The set of physics packages")
    rigorousBoundaries = PYB11property("bool", doc="Toggle if ghost nodes should be recomputed every derivative estimate")
    updateBoundaryFrequency = PYB11property("int", "updateBoundaryFrequency", "updateBoundaryFrequency", doc="Optionally update the boundary ghost nodes only on this frequency of cycles")
    verbose = PYB11property("bool", "verbose", "verbose", doc="Verbose time step information every step")
    domainDecompositionIndependent = PYB11property("bool", "domainDecompositionIndependent", "domainDecompositionIndependent", doc="Order operations to be bit perfect reproducible regardless of domain decomposition")
    cullGhostNodes = PYB11property("bool", "cullGhostNodes", "cullGhostNodes", doc="Cull ghost nodes to just active set")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, Integrator, pure_virtual=True)
PYB11inject(RestartMethods, Integrator)
