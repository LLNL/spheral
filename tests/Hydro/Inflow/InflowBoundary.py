#-------------------------------------------------------------------------------
# Inflow Boundary
#
# Creates particles at an inflow boundary
#-------------------------------------------------------------------------------
from Spheral import *
from math import *

class InflowBoundary2d(Boundary2d):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self, 
                 integrator,     # Needed to get the time
                 nodeList,       # The NodeList we're going to be applying to                
                 ):
        Boundary1d.__init__(self)
        self.integrator = integrator
        self.nodes = nodeList

        # Create a mapping for each type of field we're going to do boundary
        # conditions for to the method that does the job.
        self.fieldToMethodMap = {HydroFieldNames.mass : self.massBoundary,
                                 HydroFieldNames.velocity : self.velocityBoundary,
                                 HydroFieldNames.massDensity : self.rhoBoundary,
                                 HydroFieldNames.specificThermalEnergy : self.noopBoundary,
                                 HydroFieldNames.pressure : self.pressureBoundary,
                                 HydroFieldNames.soundSpeed : self.noopBoundary,
                                 HydroFieldNames.omegaGradh : self.omegaBoundary,
                                 HydroFieldNames.specificThermalEnergy + "0" : self.noopBoundary,
                                 HydroFieldNames.A_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.B_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.C_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.gradA_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.gradB_CRKSPH : self.noopBoundary,
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
            boundNodes.controlNodes().append(0)

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

        print "Apply GHOST to ", field.name
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
        print "Applying mass boundary"
        for i in self.nodeIDs:
            rhoi = self.answer.rho(t, pos[i].x)
            field[i] = rhoi * self.dr
            assert field[i] > 0.0
            print " --> ", i, rhoi, field[i]
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
    # Omega.
    #---------------------------------------------------------------------------
    def omegaBoundary(self, field):
        assert field.name == HydroFieldNames.omegaGradh
        for i in self.nodeIDs:
            field[i] = 1.0
        return