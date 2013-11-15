#-------------------------------------------------------------------------------
# KidderIsentropicCapsuleBoundary
#
# Enforces the boundary conditions appropriate for the Kidder isentropic ICF
# capsule.
#-------------------------------------------------------------------------------
from Spheral import *
from math import *

class KidderIsentropicCapsuleBoundary1d(Boundary1d):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self, 
                 innerBoundary,  # True => inner boundary, False => outer boundary
                 integrator,     # Needed to get the time
                 answer,         # An instance of the KidderAnalytic... class
                 nodeList,       # The NodeList we're going to be applying to
                 nrGhostNodes,   # The desired radial number of ghost nodes
                 dr0,            # Initial desired spacing in r
                 ):
        Boundary1d.__init__(self)
        self.innerBoundary = innerBoundary
        self.integrator = integrator
        self.answer = answer
        self.nodes = nodeList
        self.nr = nrGhostNodes
        self.dr0 = dr0

        # Create a mapping for each type of field we're going to do boundary
        # conditions for to the method that does the job.
        self.fieldToMethodMap = {HydroFieldNames.mass : self.massBoundary,
                                 HydroFieldNames.velocity : self.velocityBoundary,
                                 HydroFieldNames.massDensity : self.rhoBoundary,
                                 HydroFieldNames.specificThermalEnergy : self.noopBoundary,
                                 HydroFieldNames.weight : self.weightBoundary,
                                 HydroFieldNames.pressure : self.pressureBoundary,
                                 HydroFieldNames.soundSpeed : self.noopBoundary,
                                 HydroFieldNames.positionWeight : self.positionWeightBoundary,
                                 HydroFieldNames.omegaGradh : self.omegaBoundary,
                                 HydroFieldNames.specificThermalEnergy + "0" : self.noopBoundary,
                                 HydroFieldNames.A_CSPH : self.noopBoundary,
                                 HydroFieldNames.B_CSPH : self.noopBoundary,
                                 HydroFieldNames.C_CSPH : self.noopBoundary,
                                 HydroFieldNames.D_CSPH : self.noopBoundary,
                                 HydroFieldNames.gradA_CSPH : self.noopBoundary,
                                 HydroFieldNames.gradB_CSPH : self.noopBoundary,
                                 }

        return

    #---------------------------------------------------------------------------
    # Create our ghost nodes.
    #---------------------------------------------------------------------------
    def setGhostNodes(self, nodeList):
        self.addNodeList(self.nodes)
        self.nodeIDs = range(self.nodes.numNodes,
                             self.nodes.numNodes + self.nr)
        nodeList.numGhostNodes = self.nodes.numGhostNodes + self.nr
        boundNodes = self.accessBoundaryNodes(nodeList)
        for i in self.nodeIDs:
            boundNodes.ghostNodes().append(i)

        self.updateGhostNodes(nodeList)
        return

    #---------------------------------------------------------------------------
    # Update the minimal ghost node info (positions and H's).
    #---------------------------------------------------------------------------
    def updateGhostNodes(self, nodeList):

        # Fields we're going to update.
        pos = nodeList.positions()
        H = nodeList.Hfield()
        t = self.integrator.currentTime

        # Initial values for r and H.
        if self.innerBoundary:
            dr = -self.dr0 * self.answer.hfrac(t)
            ri = self.answer.rInner(t) - 0.5*dr
        else:
            dr = self.dr0 * self.answer.hfrac(t)
            ri = self.answer.rOuter(t) - 0.5*dr
        self.dr = abs(dr)
        Hi = SymTensor1d(1.0/(self.dr*nodeList.nodesPerSmoothingScale))

        # Step out from the Boundary incrementally, updating each ghost node.
        for i in self.nodeIDs:
            ri += dr
            pos[i].x = ri
            H[i] = Hi

        mass = nodeList.mass()
        self.massBoundary(mass)
        print "mass: ", [mass[i] for i in self.nodeIDs]

        nodeList.neighbor().updateNodes()
        return

    #---------------------------------------------------------------------------
    # Update ghost field values.
    # This is specialized for each of the fields that gets passed in for a 
    # hydro run.
    #---------------------------------------------------------------------------
    def applyGhostBoundary(self, field):

        if field.name not in self.fieldToMethodMap:
            print "KidderBoundary WARNING: unable to apply boundary condition to %s" % field.name
            return

        self.fieldToMethodMap[field.name](field)
        return

    #---------------------------------------------------------------------------
    # Violation nodes -- not defined in this case.
    #---------------------------------------------------------------------------
    def setViolationNodes(self, nodeList):
        return

    def updateViolationNodes(self, nodeList):
        return

    def enforceBoundary(self, field):
        return

    #---------------------------------------------------------------------------
    # No op boundary -- do nothing.
    #---------------------------------------------------------------------------
    def noopBoundary(self, field):
        return

    #---------------------------------------------------------------------------
    # Mass.
    #---------------------------------------------------------------------------
    def massBoundary(self, field):
        assert field.name == HydroFieldNames.mass
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        for i in self.nodeIDs:
            rhoi = self.answer.rho(t, pos[i].x)
            field[i] = rhoi * self.dr
            assert field[i] > 0.0
        return

    #---------------------------------------------------------------------------
    # Velocity.
    #---------------------------------------------------------------------------
    def velocityBoundary(self, field):
        assert field.name == HydroFieldNames.velocity
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        for i in self.nodeIDs:
            field[i].x = self.answer.vr(t, pos[i].x)
        return

    #---------------------------------------------------------------------------
    # Mass density.
    #---------------------------------------------------------------------------
    def rhoBoundary(self, field):
        assert field.name == HydroFieldNames.massDensity
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        for i in self.nodeIDs:
            field[i] = self.answer.rho(t, pos[i].x)
        return

    #---------------------------------------------------------------------------
    # Weight.
    #---------------------------------------------------------------------------
    def weightBoundary(self, field):
        assert field.name == HydroFieldNames.weight
        for i in self.nodeIDs:
            field[i] = self.dr
        return

    #---------------------------------------------------------------------------
    # Pressure.
    #---------------------------------------------------------------------------
    def pressureBoundary(self, field):
        assert field.name == HydroFieldNames.pressure
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        for i in self.nodeIDs:
            field[i] = self.answer.P(t, pos[i].x)
        return

    #---------------------------------------------------------------------------
    # Position weight.
    #---------------------------------------------------------------------------
    def positionWeightBoundary(self, field):
        assert field.name == HydroFieldNames.positionWeight
        for i in self.nodeIDs:
            field[i] = 1.0
        return

    #---------------------------------------------------------------------------
    # Omega.
    #---------------------------------------------------------------------------
    def omegaBoundary(self, field):
        assert field.name == HydroFieldNames.omegaGradh
        for i in self.nodeIDs:
            field[i] = 1.0
        return

#-------------------------------------------------------------------------------
# This class directly sets the state of a set of internal nodes, rather than
# create ghost nodes.  This is necessary for the inner boundary, where ghost
# nodes would often need a negative mass density.
#-------------------------------------------------------------------------------
class KidderIsentropicCapsuleEnforcementBoundary1d(Physics1d):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self, 
                 integrator,     # Needed to get the time
                 answer,         # An instance of the KidderAnalytic... class
                 nodeList,       # The NodeList we're going to be applying to
                 nodeIDs,        # The set of internal nodes we're going to control
                 interiorNodeIDs,# The interior nodes not controlled by a boundary condition
                 hinitial,       # Initial smoothing scale that we're going to impose
                 hmin,           # minimum smoothing scale
                 hmax,           # maximum smoothing scale
                 dr0,            # Intial node spacing
                 ):
        Physics1d.__init__(self)
        self.integrator = integrator
        self.answer = answer
        self.nodes = nodeList
        self.nodeIDs = nodeIDs
        self.interiorNodeIDs = interiorNodeIDs
        self.hinitial = hinitial
        self.hmin = hmin
        self.hmax = hmax
        self.dr0 = dr0
        self.gamma = self.answer.gamma
        self.gamma1 = self.answer.gamma1

        # Store some of the initial values for nodes we're controlling.
        pos = self.nodes.positions()
        eps = self.nodes.specificThermalEnergy()
        self.initialRadii = [pos[i].x for i in self.nodeIDs]
        self.initialEps = [eps[i] for i in self.nodeIDs]

        # Store the total mass in interior nodes.
        mass = self.nodes.mass()
        self.interiorMass = mpi.allreduce(sum([mass[i] for i in interiorNodeIDs]), mpi.SUM)

        return

    #---------------------------------------------------------------------------
    # Override the state values for the physics variables, imposing the analytic
    # solution.
    #---------------------------------------------------------------------------
    def overrideState(self, state):

        t = self.integrator.currentTime

        # Extract the state we're going to set.
        Key = pair_NodeList1d_string
        rhoKey = Key(self.nodes, HydroFieldNames.massDensity)
        posKey = Key(self.nodes, HydroFieldNames.position)
        epsKey = Key(self.nodes, HydroFieldNames.specificThermalEnergy)
        velKey = Key(self.nodes, HydroFieldNames.velocity)
        Hkey = Key(self.nodes, HydroFieldNames.H)
        Pkey = Key(self.nodes, HydroFieldNames.pressure)
        csKey = Key(self.nodes, HydroFieldNames.soundSpeed)
        wKey = Key(self.nodes, HydroFieldNames.weight)
        mKey = Key(self.nodes, HydroFieldNames.mass)
        rho = state.scalarField(rhoKey)
        pos = state.vectorField(posKey)
        eps = state.scalarField(epsKey)
        vel = state.vectorField(velKey)
        H = state.symTensorField(Hkey)
        P = state.scalarField(Pkey)
        cs = state.scalarField(csKey)
        weight = state.scalarField(wKey)
        mass = state.scalarField(mKey)
        
        hfrac = self.answer.hfrac(t)
        Hi = SymTensor1d(1.0/max(self.hmin, min(self.hmax, self.hinitial*hfrac)))

        # Now set these variables for the nodes we're controlling.
        for k in xrange(len(self.nodeIDs)):
            i = self.nodeIDs[k]
            ri = self.initialRadii[k]*hfrac
            rhoi = mass[i]/(self.dr0*hfrac)   # self.answer.rho(t, ri)
            epsi = self.initialEps[k]/(hfrac*hfrac)
            Pi = self.gamma1*rhoi*epsi
            #Pi = self.answer.P(t, ri)
            #epsi = Pi/(self.gamma1*rhoi)
            rho[i] = rhoi
            pos[i].x = ri
            eps[i] = epsi
            vel[i].x = self.answer.vr(t, ri)
            H[i] = Hi
            P[i] = Pi
            cs[i] = sqrt(self.gamma * self.gamma1 * epsi)
            weight[i] = mass[i]/rhoi
        
        n = self.nodes.neighbor()
        n.updateNodes()        

        return

    #---------------------------------------------------------------------------
    # Override the derivative values for the physics variables, imposing the
    # analytic solution.
    #---------------------------------------------------------------------------
    def overrideDerivatives(self, derivs):

        t = self.integrator.currentTime

        # Extract the state we're going to set.
        Key = pair_NodeList1d_string
        DxDtKey = Key(self.nodes, "delta " + HydroFieldNames.position)
        DvDtKey = Key(self.nodes, "delta " + HydroFieldNames.velocity)
        DvDxKey = Key(self.nodes, HydroFieldNames.velocityGradient)
        internalDvDxKey = Key(self.nodes, HydroFieldNames.internalVelocityGradient)
        DHDtKey = Key(self.nodes, "delta " + HydroFieldNames.H)
        DrhoDtKey = Key(self.nodes, "delta " + HydroFieldNames.massDensity)
        DepsDtKey = Key(self.nodes, "delta " + HydroFieldNames.specificThermalEnergy)
        rhoSumKey = Key(self.nodes, "new " + HydroFieldNames.massDensity)
        HidealKey = Key(self.nodes, "new " + HydroFieldNames.H)
        XSPHDeltaVkey = Key(self.nodes, HydroFieldNames.XSPHDeltaV)
        DxDt = derivs.vectorField(DxDtKey)
        DvDt = derivs.vectorField(DvDtKey)
        DvDx = derivs.tensorField(DvDxKey)
        internalDvDx = derivs.tensorField(internalDvDxKey)
        DHDt = derivs.symTensorField(DHDtKey)
        DrhoDt = derivs.scalarField(DrhoDtKey)
        DepsDt = derivs.scalarField(DepsDtKey)
        rhoSum = derivs.scalarField(rhoSumKey)
        Hideal = derivs.symTensorField(HidealKey)
        XSPHDeltaV = derivs.vectorField(XSPHDeltaVkey)
        mass = self.nodes.mass()

        hfrac = self.answer.hfrac(t)
        Hi = SymTensor1d(1.0/max(self.hmin, min(self.hmax, self.hinitial*hfrac)))

        # Now set these variables for the nodes we're controlling.
        for k in xrange(len(self.nodeIDs)):
            i = self.nodeIDs[k]
            ri = self.initialRadii[k]*hfrac
            rhoi = mass[i]/(self.dr0*hfrac)   # self.answer.rho(t, ri)
            epsi = self.initialEps[k]/(hfrac*hfrac)
            Pi = self.gamma1*rhoi*epsi
            #Pi = self.answer.P(t, ri)
            DvDxi = Tensor1d(self.answer.DvrDr(t, ri))
            rhoDoti = self.answer.rhoDot(t, ri)
            Pdoti = self.answer.Pdot(t, ri)

            DxDt[i].x = self.answer.vr(t, ri)
            DvDt[i].x = self.answer.vrDot(t, ri)
            DvDx[i] = DvDxi
            internalDvDx[i] = DvDxi
            DHDt[i] = Hi*rhoDoti/rhoi
            DrhoDt[i] = rhoDoti
            DepsDt[i] = (Pdoti - Pi*rhoDoti/rhoi)/(self.gamma1*rhoi)
            rhoSum[i] = rhoi
            Hideal[i] = Hi
            XSPHDeltaV.Zero()
        
        return

    #---------------------------------------------------------------------------
    # Physics::evaluateDerivatives
    #---------------------------------------------------------------------------
    def evaluateDerivatives(self, currentTime, dt, dataBase, state, derivs):
        self.overrideDerivatives(derivs)
        return

    #---------------------------------------------------------------------------
    # Physics::dt
    #---------------------------------------------------------------------------
    def dt(self, dataBase, state, derivs, currentTime):
        return pair_double_string(1.0e50, "No vote.")
        return

    #---------------------------------------------------------------------------
    # Physics::registerState
    #---------------------------------------------------------------------------
    def registerState(self, dataBase, state):
        return

    #---------------------------------------------------------------------------
    # Physics::registerDerivatives
    #---------------------------------------------------------------------------
    def registerDerivatives(self, dataBase, derivs):
        return

    #---------------------------------------------------------------------------
    # Physics::initialize
    #---------------------------------------------------------------------------
    def initialize(self, currentTime, dt, dataBase, state, derivs):
        self.overrideState(state)
        self.overrideDerivatives(derivs)
        return

    #---------------------------------------------------------------------------
    # Physics::finalize
    #---------------------------------------------------------------------------
    def finalize(self, currentTime, dt, dataBase, state, derivs):
        self.overrideState(state)
        self.overrideDerivatives(derivs)
        n = self.nodes.neighbor()
        n.updateNodes()        
        return

    #---------------------------------------------------------------------------
    # Physics::finalizeDerivatives
    #---------------------------------------------------------------------------
    def finalizeDerivatives(self, currentTime, dt, dataBase, state, derivs):
        self.overrideDerivatives(derivs)
        return

    #---------------------------------------------------------------------------
    # Physics::postStateUpdate
    #---------------------------------------------------------------------------
    def postStateUpdate(self, dataBase, state, derivs):

##         # Take a snapshot of the internal energies.
##         Key = pair_NodeList1d_string
##         mKey = Key(self.nodes, HydroFieldNames.mass)
##         epsKey = Key(self.nodes, HydroFieldNames.specificThermalEnergy)
##         mass = state.scalarField(mKey)
##         eps = state.scalarField(epsKey)
##         self.preEps = [eps[i] for i in self.nodeIDs]

        self.overrideState(state)
        n = self.nodes.neighbor()
        n.updateNodes()        

##         # Redistribute the energy change we imposed to the interior nodes.
##         deltaE = sum([mass[i] * (eps[i] - self.preEps[k]) for i, k in zip(self.nodeIDs, range(len(self.nodeIDs)))])
##         print "Energy changed by ", deltaE
##         deltaEps = -deltaE/self.interiorMass
##         for i in self.interiorNodeIDs:
##             eps[i] += deltaEps

        return
