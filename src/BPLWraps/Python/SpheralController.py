#------------------------------------------------------------------------------
# A simple class to control simulation runs for Spheral.
#------------------------------------------------------------------------------
import sys
from Spheral import AccessType
from Spheral import State1d, State2d, State3d
from Spheral import StateDerivatives1d, StateDerivatives2d, StateDerivatives3d
from Spheral import iterateIdealH
from Spheral import TableKernel1d, TableKernel2d, TableKernel3d
from Spheral import BSplineKernel1d, BSplineKernel2d, BSplineKernel3d
from Spheral import RestartableObject, RestartRegistrar
from SpheralTimer import SpheralTimer
from SpheralConservation import SpheralConservation
from ExtendFlatFileIO import FlatFileIO
from GzipFileIO import GzipFileIO

import loadmpi, fakempi
mpi, rank, nprocs = loadmpi.loadmpi()

class SpheralController(RestartableObject):

    #--------------------------------------------------------------------------
    # Constuctor.
    #--------------------------------------------------------------------------
    def __init__(self, integrator, kernel,
                 statsStep = 1,
                 printStep = 1,
                 redistributeStep = None,
                 restartStep = None,
                 restartBaseName = "restart",
                 restartObjects = [],
                 restartFileConstructor = GzipFileIO,
                 initializeMassDensity = False):
        RestartableObject.__init__(self)
        self.integrator = integrator
        self.kernel = kernel
        self.restartObjects = restartObjects
        self.restartFileConstructor = restartFileConstructor
        self.initializeMassDensity = initializeMassDensity

        # Determine the dimensionality of this run, based on the integrator.
        self.dim = "%id" % self.integrator.dataBase().nDim

        # If this is a parallel run, automatically construct and insert
        # a DistributedBoundaryCondition into each physics package.
        self.insertDistributedBoundary(integrator.physicsPackages())

        # Generic initialization work.
        self.reinitializeProblem(restartBaseName,
                                 statsStep = statsStep,
                                 printStep = printStep,
                                 redistributeStep = redistributeStep,
                                 restartStep = restartStep)

        # Add the dynamic redistribution object to the controller.
        self.addRedistributeNodes(self.kernel)

        return

    #--------------------------------------------------------------------------
    # (Re)initialize the current problem (and controller state).
    # This method is intended to be called before the controller begins a new
    # problem from time 0.
    #--------------------------------------------------------------------------
    def reinitializeProblem(self, restartBaseName,
                            initialTime = 0.0,
                            statsStep = 1,
                            printStep = 1,
                            redistributeStep = None,
                            restartStep = None):

        # Intialize the cycle count.
        self.totalSteps = 0

        # Construct a timer to track the cycle step time.
        self.stepTimer = SpheralTimer("Time per integration cycle.")

        # Construct a fresh conservation check object.
        self.conserve = SpheralConservation(self.integrator.dataBase(),
                                            self.integrator.physicsPackages())

        # Prepare an empty set of periodic work.
        self._periodicWork = []
        
        # Set the restart file base name.
        self.setRestartBaseName(restartBaseName)
        
        # Sum the mass density if requested
        self.computeMassDensity()

        # Set the simulation time.
        self.integrator.setCurrentTime(initialTime)
##         state = eval("State%s(self.integrator.dataBase(), self.integrator.physicsPackages())" % (self.dim))
##         derivs = eval("StateDerivatives%s(self.integrator.dataBase(), self.integrator.physicsPackages())" % (self.dim))
##         self.integrator.initialize(state, derivs)
##         self.integrator.initialize(state, derivs)

        # Set up the default periodic work.
        self.appendPeriodicWork(self.printCycleStatus, printStep)
        self.appendPeriodicWork(self.updateConservation, statsStep)
        self.appendPeriodicWork(self.updateDomainDistribution, redistributeStep)
        self.appendPeriodicWork(self.updateRestart, restartStep)

        return

    #--------------------------------------------------------------------------
    # Initialize the mass density, weight, and correction terms.
    #--------------------------------------------------------------------------
    def computeMassDensity(self):
        db = self.integrator.dataBase()
        state = eval("State%s(db, self.integrator.physicsPackages())" % (self.dim))
        derivs = eval("StateDerivatives%s(db, self.integrator.physicsPackages())" % (self.dim))
        self.integrator.initialize(state, derivs)
        db.updateFluidWeight()
        self.integrator.initialize(state, derivs)
        db.updateFluidCorrections(True)
        if self.initializeMassDensity:
            db.updateFluidMassDensity()
        self.integrator.applyGhostBoundaries(state, derivs)
        return

    #--------------------------------------------------------------------------
    # Set the restart base name.
    #--------------------------------------------------------------------------
    def setRestartBaseName(self, name):
        self.restartBaseName = name

        # If we're running parallel then add the domain info to the restart
        # base name.
        try:
            import mpi
            self.restartBaseName += '_rank%i_of_%idomains' % (mpi.rank, mpi.procs)
        except ImportError:
            pass

        return

    #--------------------------------------------------------------------------
    # Return the current time.
    #--------------------------------------------------------------------------
    def time(self):
        return self.integrator.currentTime

    #--------------------------------------------------------------------------
    # Return the last timestep.
    #--------------------------------------------------------------------------
    def lastDt(self):
        return self.integrator.lastDt

    #--------------------------------------------------------------------------
    # Smooth the physical variables.
    #--------------------------------------------------------------------------
    def smoothState(self, smoothIters=1):
        import FieldOperations
        db = self.integrator.dataBase()
        scalarSmooth = eval("FieldOperations.smoothScalarFields%id" % db.nDim)
        vectorSmooth = eval("FieldOperations.smoothVectorFields%id" % db.nDim)
        tensorSmooth = eval("FieldOperations.smoothSymTensorFields%id" % db.nDim)
        for iter in xrange(smoothIters):
            state = eval("State%id(db, self.integrator.physicsPackages())" % db.nDim)
            derivs = eval("StateDerivatives%id(db, self.integrator.physicsPackages())" % db.nDim)
            self.integrator.setGhostNodes()
            self.integrator.applyGhostBoundaries(state, derivs)
            smoothedVelocity = vectorSmooth(db.fluidVelocity,
                                            db.fluidPosition,
                                            db.fluidWeight,
                                            db.fluidMass,
                                            db.fluidMassDensity,
                                            db.fluidHfield,
                                            self.kernel)
            smoothedSpecificThermalEnergy = scalarSmooth(db.fluidSpecificThermalEnergy,
                                                         db.fluidPosition,
                                                         db.fluidWeight,
                                                         db.fluidMass,
                                                         db.fluidMassDensity,
                                                         db.fluidHfield,
                                                         self.kernel)
            smoothedHfield = tensorSmooth(db.fluidHfield,
                                          db.fluidPosition,
                                          db.fluidWeight,
                                          db.fluidMass,
                                          db.fluidMassDensity,
                                          db.fluidHfield,
                                          self.kernel)
            db.fluidVelocity.assignFields(smoothedVelocity)
            db.fluidSpecificThermalEnergy.assignFields(smoothedSpecificThermalEnergy)
            db.fluidHfield.assignFields(smoothedHfield)
            for nodeList in db.fluidNodeLists():
                nodeList.neighbor().updateNodes()
            db.updateConnectivityMap()
            for nodeList in db.fluidNodeLists():
                nodeList.updateWeight(db.connectivityMap())

        return

    #--------------------------------------------------------------------------
    # Advance the system to the given simulation time.  The user can also
    # specify a max number of steps to take.
    #--------------------------------------------------------------------------
    def advance(self, goalTime, maxSteps=None):
        currentSteps = 0
        while (self.time() < goalTime and
               (maxSteps == None or currentSteps < maxSteps)):
            self.stepTimer.start()
            self.integrator.step(goalTime)
            self.stepTimer.stop()
            currentSteps = currentSteps + 1
            self.totalSteps = self.totalSteps + 1

            # Do any periodic work, as determined by the number of steps.
            for tup in self._periodicWork:
                method = tup[0]
                frequency = tup[1]
                if frequency is not None and self.totalSteps % frequency == 0:
                    method(self.totalSteps, self.time(), self.lastDt())

        # Print how much time was spent per integration cycle.
        self.stepTimer.printStatus()

        return

    #--------------------------------------------------------------------------
    # Step method, where the user can specify to take a given number of steps.
    #--------------------------------------------------------------------------
    def step(self, steps=1):
        self.advance(1e40, steps)

    #--------------------------------------------------------------------------
    # Add a (method, frequency) tuple to the set 'o stuff the controller should
    # call during advance.
    #--------------------------------------------------------------------------
    def appendPeriodicWork(self, method, frequency):
        self._periodicWork.append((method, frequency))

    #--------------------------------------------------------------------------
    # Remove all instances of given method from the set of periodic work.
    #--------------------------------------------------------------------------
    def removePeriodicWork(self, method):
        cache = self._periodicWork[:]
        self._periodicWork = []
        for tup in cache:
            if tup[0] != method:
                self._periodicWork.append(tup)

    #--------------------------------------------------------------------------
    # Change the frequency at which the given method is called in periodic
    # work.
    #--------------------------------------------------------------------------
    def setFrequency(self, method, frequency):
        i = 0
        while (i < len(self._periodicWork) and
               self._periodicWork[i][0] != method):
            i += 1

        if i == len(self._periodicWork):
            print "Error, could not find periodic work calling ", method
            return
        else:
            self._periodicWork[i] = (method, frequency)

    #--------------------------------------------------------------------------
    # Get the frequency at which the given method is called in periodic
    # work.
    #--------------------------------------------------------------------------
    def getFrequency(self, method):
        i = 0
        while (i < len(self._periodicWork) and
               self._periodicWork[i][0] != method):
            i += 1

        if i == len(self._periodicWork):
            raise "Error, could not find periodic work calling %s" % str(method)

        return self._periodicWork[i][1]

    #--------------------------------------------------------------------------
    # A method for printing the cycle information.
    #--------------------------------------------------------------------------
    def printCycleStatus(self, cycle, Time, dt):
        print "Cycle=%i, \tTime=%g, \tTimeStep=%g" % (cycle, Time, dt)
        return

    #--------------------------------------------------------------------------
    # Periodically update the conservation statistics.
    #--------------------------------------------------------------------------
    def updateConservation(self, cycle, Time, dt):
        self.conserve.updateHistory(cycle, Time)
        return
    
    #--------------------------------------------------------------------------
    # Periodically drop a restart file.
    #--------------------------------------------------------------------------
    def updateRestart(self, cycle, Time, dt):
        self.dropRestartFile()
        return

    #--------------------------------------------------------------------------
    # Periodically redistribute the nodes between domains.
    #--------------------------------------------------------------------------
    def updateDomainDistribution(self, cycle, Time, dt):
        if self.redistribute:
            self.redistributeTimer.start()
            self.redistribute.redistributeNodes(self.integrator.dataBase(),
                                                self.integrator.uniqueBoundaryConditions())
            self.redistributeTimer.stop()
            self.redistributeTimer.printStatus()
        return

    #--------------------------------------------------------------------------
    # Find the name associated with the given object.
    #--------------------------------------------------------------------------
    def findname(thing):
        for mod in sys.modules.values():
            for name, val in mod.__dict__.items():
                if val is thing:
                    return name

    #--------------------------------------------------------------------------
    # Iterate over all the restartable objects and drop their state to a file.
    #--------------------------------------------------------------------------
    def dropRestartFile(self):

        # First find out if the requested directory exists.
        import os
        dire = os.path.dirname(os.path.abspath(self.restartBaseName))
        if not os.path.exists(dire):
            raise RuntimeError("Directory %s does not exist or is inaccessible." %
                               dire)

        # Now we can invoke the restart!
        import time
        start = time.clock()
        fileName = self.restartBaseName + "_cycle%i" % self.totalSteps
        file = self.restartFileConstructor(fileName, AccessType.Create)
        RestartRegistrar.instance().dumpState(file)
        print "Wrote restart file in %0.2f seconds" % (time.clock() - start)

        file.close()
        return

    #--------------------------------------------------------------------------
    # Iterate over all the restartable objects and restore their state.
    #--------------------------------------------------------------------------
    def loadRestartFile(self, restoreCycle,
                        frameDict=None):

        # Find out if the requested file exists.
        import os
        fileName = self.restartBaseName + "_cycle%i" % restoreCycle
        if self.restartFileConstructor is GzipFileIO:
            fileName += ".gz"
        if not os.path.exists(fileName):
            raise RuntimeError("File %s does not exist or is inaccessible." %
                               fileName)

        # Read that sucker.
        print 'Reading from restart file', fileName
        import time
        start = time.clock()
        if self.restartFileConstructor is GzipFileIO:
            file = self.restartFileConstructor(fileName, AccessType.Read,
                                               readToMemory = True)
        else:
            file = self.restartFileConstructor(fileName, AccessType.Read)
        RestartRegistrar.instance().restoreState(file)
        print "Finished: required %0.2f seconds" % (time.clock() - start)

        file.close()
        return

    #--------------------------------------------------------------------------
    # Generate a label for our restart info.
    #--------------------------------------------------------------------------
    def label(self):
        return "Controller"

    #--------------------------------------------------------------------------
    # Store the controllers necessary data.
    #--------------------------------------------------------------------------
    def dumpState(self, file, path):
        controlPath = path + "/self"
        file.writeObject(self.totalSteps, controlPath + "/totalSteps")
        return

    #--------------------------------------------------------------------------
    # Restore the controllers state.
    #--------------------------------------------------------------------------
    def restoreState(self, file, path):
        controlPath = path + "/self"
        self.totalSteps = file.readObject(controlPath + "/totalSteps")
        return

    #--------------------------------------------------------------------------
    # If necessary create and add a distributed boundary condition to each
    # physics package
    #--------------------------------------------------------------------------
    def insertDistributedBoundary(self, physicsPackages):
        import Boundary

        # This is the list of boundary types that need to precede the distributed
        # boundary, since their ghost nodes need to be communicated.
        precedeDistributed = [Boundary.PeriodicBoundary1d,
                              Boundary.PeriodicBoundary2d,
                              Boundary.PeriodicBoundary3d,
                              Boundary.ConstantBoundary1d,
                              Boundary.ConstantBoundary2d,
                              Boundary.ConstantBoundary3d,
                              Boundary.CylindricalBoundary,
                              Boundary.SphericalBoundary]

        # Check if this is a pyMPI process or not.
        if isinstance(mpi, fakempi.fakempi):
            self.domainbc = None

        # If this is a parallel run, then we need to create a distributed
        # boundary condition and insert it into the list of boundaries for each physics
        # package.
        else:
            import Distributed
            self.domainbc = eval("Distributed.NestedGridDistributedBoundary%s.instance()" % self.dim)

            # Iterate over each of the physics packages.
            for package in physicsPackages:

                # Make a copy of the current set of boundary conditions for this package,
                # and clear out the set in the physics package.
                boundaryConditions = list(package.boundaryConditions())
                package.clearBoundaries()

                # Sort the boundary conditions into two lists: those that need
                # to precede the distributed boundary condition, and those that
                # should follow it.
                precedeBoundaries = []
                followBoundaries = []
                for boundary in boundaryConditions:
                    precede = False
                    for btype in precedeDistributed:
                        if isinstance(boundary, btype):
                            precede = True
                    if precede:
                        precedeBoundaries.append(boundary)
                    else:
                        followBoundaries.append(boundary)

                assert len(precedeBoundaries) + len(followBoundaries) == len(boundaryConditions)
                for boundary in boundaryConditions:
                    assert (boundary in precedeBoundaries) or (boundary in followBoundaries)

                # Add in the distributed boundary condition.
                #precedeBoundaries.append(self.domainbc)
                followBoundaries.append(self.domainbc)

                # Now put the boundaries back into the package.
                for bc in precedeBoundaries + followBoundaries:
                    package.appendBoundary(bc)

        # That's it.
        if not (self.domainbc is None):
            for package in physicsPackages:
                assert package.haveBoundary(self.domainbc)
        return

    #--------------------------------------------------------------------------
    # If this is a parallel run, then add a domain repartitioner to dynamically
    # balance the workload.
    #--------------------------------------------------------------------------
    def addRedistributeNodes(self, W):
        assert W.kernelExtent > 0.0
        self.redistribute = None
        self.redistributeTimer = SpheralTimer("Time for redistributing nodes.")
        if not isinstance(mpi, fakempi.fakempi):
            import Distributed
            if mpi.procs > 1:
                try:
                    #self.redistribute = eval("Distributed.ParmetisRedistributeNodes%s(W.kernelExtent())" % self.dim)
                    #self.redistribute = eval("Distributed.SortAndDivideRedistributeNodes%s(W.kernelExtent())" % self.dim)
                    self.redistribute = eval("Distributed.PeanoHilbertOrderRedistributeNodes%s(W.kernelExtent())" % self.dim)
                    self.redistribute.workBalance = True
                except:
                    print "Warning: this appears to be a parallel run, but Controller cannot construct"
                    print "         dynamic redistributer."
                    pass
        return

    #--------------------------------------------------------------------------
    # Iteratively set the H tensors, until the desired convergence criteria
    # are met.
    #--------------------------------------------------------------------------
    def iterateIdealH(self,
                      hydro = None,
                      hmin = None,
                      hmax = None,
                      hminratio = None,
                      maxIdealHIterations = 100,
                      idealHTolerance = 1.0e-4):
        print "SpheralController: Initializing H's..."

        db = self.integrator.dataBase()
        bcs = self.integrator.uniqueBoundaryConditions()
        if hydro:
            WT = hydro.kernel()
        else:
            WT = eval("TableKernel%s(BSplineKernel%s(), 1000)" % (self.dim, self.dim))
        iterateIdealH(db, bcs, WT, maxIdealHIterations, idealHTolerance, 0.0, False, False)
        self.computeMassDensity()

        return
