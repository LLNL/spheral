#------------------------------------------------------------------------------
# A simple class to control simulation runs for Spheral.
#------------------------------------------------------------------------------
import sys, os, gc, warnings, mpi

from SpheralModules.Spheral import FileIOSpace
from SpheralModules.Spheral.DataOutput import RestartableObject, RestartRegistrar
from SpheralModules.Spheral import BoundarySpace
from SpheralModules.Spheral.FieldSpace import *
from SpheralModules.Spheral.FileIOSpace import *
from SpheralModules import Timer
from SpheralTimer import SpheralTimer
from SpheralConservation import SpheralConservation
from GzipFileIO import GzipFileIO
from SpheralTestUtilities import globalFrame
from NodeGeneratorBase import ConstantRho
from findLastRestart import findLastRestart

from spheralDimensions import spheralDimensions
dims = spheralDimensions()
for dim in dims:
    exec("""
from SpheralModules.Spheral import State%(dim)sd
from SpheralModules.Spheral import StateDerivatives%(dim)sd
from SpheralModules.Spheral import iterateIdealH%(dim)sd
from SpheralModules.Spheral.NodeSpace import ASPHSmoothingScale%(dim)sd
from SpheralModules.Spheral.NodeSpace import SPHSmoothingScale%(dim)sd
from SpheralModules.Spheral.KernelSpace import TableKernel%(dim)sd
from SpheralModules.Spheral.KernelSpace import BSplineKernel%(dim)sd
from SpheralModules import vector_of_Physics%(dim)sd
""" % {"dim" : dim})

class SpheralController:

    #--------------------------------------------------------------------------
    # Constuctor.
    #--------------------------------------------------------------------------
    def __init__(self, integrator, kernel,
                 statsStep = 1,
                 printStep = 1,
                 garbageCollectionStep = 100,
                 redistributeStep = None,
                 restartStep = None,
                 restartBaseName = "restart",
                 restartObjects = [],
                 restartFileConstructor = SiloFileIO,
                 restoreCycle = None,
                 initializeDerivatives = False,
                 vizBaseName = None,
                 vizDir =  ".",
                 vizStep = None,
                 vizTime = None,
                 vizMethod = None,
                 vizFields = [],
                 vizFieldLists = [],
                 vizGhosts = False,
                 vizDerivs = False,
                 initialTime = 0.0,
                 SPH = False,
                 skipInitialPeriodicWork = False,
                 iterateInitialH = True,
                 numHIterationsBetweenCycles = 0):
        self.restart = RestartableObject(self)
        self.integrator = integrator
        self.kernel = kernel
        self.restartObjects = restartObjects
        self.restartFileConstructor = restartFileConstructor
        self.SPH = SPH
        self.numHIterationsBetweenCycles = numHIterationsBetweenCycles
        self._break = False

        # Determine the dimensionality of this run, based on the integrator.
        self.dim = "%id" % self.integrator.dataBase().nDim

        # Determine the visualization method.
        dumpPhysicsStatePoints, dumpPhysicsStateCells = None, None
        if self.dim == "1d":
            from Spheral1dVizDump import dumpPhysicsState as dumpPhysicsStatePoints
        else:
            from SpheralPointmeshSiloDump import dumpPhysicsState as dumpPhysicsStatePoints
            from SpheralVoronoiSiloDump import dumpPhysicsState as dumpPhysicsStateCells
        if vizMethod:
            dumpPhysicsStatePoints = vizMethod
        self.vizMethodPoints = dumpPhysicsStatePoints
        self.vizMethodCells = dumpPhysicsStateCells
        self.vizGhosts = vizGhosts
        self.vizDerivs = vizDerivs

        # If this is a parallel run, automatically construct and insert
        # a DistributedBoundaryCondition into each physics package.
        self.insertDistributedBoundary(integrator.physicsPackages())

        # Should we look for the last restart set?
        if restoreCycle == -1:
            restoreCycle = findLastRestart(restartBaseName)

        # Generic initialization work.
        self.reinitializeProblem(restartBaseName,
                                 vizBaseName,
                                 initialTime = initialTime,
                                 statsStep = statsStep,
                                 printStep = printStep,
                                 garbageCollectionStep = garbageCollectionStep,
                                 redistributeStep = redistributeStep,
                                 restartStep = restartStep,
                                 restoreCycle = restoreCycle,
                                 initializeDerivatives = initializeDerivatives,
                                 vizDir = vizDir,
                                 vizStep = vizStep,
                                 vizTime = vizTime,
                                 vizFields = vizFields,
                                 vizFieldLists = vizFieldLists,
                                 skipInitialPeriodicWork = skipInitialPeriodicWork,
                                 iterateInitialH = True)

        # Read the restart information if requested.
        if not restoreCycle is None:
            print "Reading from restart file ", restoreCycle
            self.loadRestartFile(restoreCycle)
        
        return

    #--------------------------------------------------------------------------
    # (Re)initialize the current problem (and controller state).
    # This method is intended to be called before the controller begins a new
    # problem from time 0.
    #--------------------------------------------------------------------------
    def reinitializeProblem(self, restartBaseName, vizBaseName,
                            initialTime = 0.0,
                            statsStep = 1,
                            printStep = 1,
                            garbageCollectionStep = 100,
                            redistributeStep = None,
                            restartStep = None,
                            restoreCycle = None,
                            initializeDerivatives = False,
                            vizDir = None,
                            vizStep = None,
                            vizTime = None,
                            vizFields = [],
                            vizFieldLists = [],
                            skipInitialPeriodicWork = False,
                            iterateInitialH = True):

        # Intialize the cycle count.
        self.totalSteps = 0

        # Construct a timer to track the cycle step time.
        self.stepTimer = SpheralTimer("Time per integration cycle.")

        # Prepare an empty set of periodic work.
        self._periodicWork = []
        self._periodicTimeWork = []
        
        # Set the restart file base name.
        self.setRestartBaseName(restartBaseName)
        
        # Set the simulation time.
        self.integrator.currentTime = initialTime

        # Create ghost nodes for the physics packages to initialize with.
        self.integrator.setGhostNodes()

        # Initialize the integrator and packages.
        db = self.integrator.dataBase()
        packages = self.integrator.physicsPackages()
        for package in packages:
            package.initializeProblemStartup(db)
            for bc in package.boundaryConditions():
                bc.initializeProblemStartup()
        state = eval("State%s(db, packages)" % (self.dim))
        derivs = eval("StateDerivatives%s(db, packages)" % (self.dim))
        self.integrator.preStepInitialize(state, derivs)
        self.integrator.initializeDerivatives(initialTime, 0.0, state, derivs)

        # If requested, initialize the derivatives.
        if initializeDerivatives:
            self.integrator.evaluateDerivatives(initialTime, 0.0, db, state, derivs)

        # If we're starting from scratch, initialize the H tensors.
        if restoreCycle is None and not skipInitialPeriodicWork and iterateInitialH:
            self.iterateIdealH()

        # Set up the default periodic work.
        self.appendPeriodicWork(self.printCycleStatus, printStep)
        self.appendPeriodicWork(self.garbageCollection, garbageCollectionStep)
        self.appendPeriodicWork(self.updateConservation, statsStep)
        self.appendPeriodicWork(self.updateRestart, restartStep)

        # Add the dynamic redistribution object to the controller.
        self.addRedistributeNodes(self.kernel)

        # Set up any visualization required.
        if not vizBaseName is None:
            assert not vizDir is None
            assert not (vizStep is None and vizTime is None)
            self.addVisualizationDumps(vizBaseName, vizDir, vizStep, vizTime,
                                       vizFields, vizFieldLists)

        # Construct a fresh conservation check object.
        # Hopefully by this time all packages have initialized their own extra energy bins.
        self.conserve = SpheralConservation(self.integrator.dataBase(),
                                            self.integrator.physicsPackages())

        # Force the periodic work to fire at problem initalization.
        if (not skipInitialPeriodicWork) and (restoreCycle is None):
            self.doPeriodicWork(force=True)

        # We add this one after forcing periodic work so it's not always fired right at the beginning of a calculation.
        self.appendPeriodicWork(self.updateDomainDistribution, redistributeStep)

        return

    #--------------------------------------------------------------------------
    # Set the restart base name.
    #--------------------------------------------------------------------------
    def setRestartBaseName(self,
                           name,
                           rank = mpi.rank,
                           procs = mpi.procs):
        self.restartBaseName = name

        # If we're running parallel then add the domain info to the restart
        # base name.
        if procs > 1:
            self.restartBaseName += '_rank%i_of_%idomains' % (rank, procs)

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
        db = self.integrator.dataBase()
        scalarSmooth = eval("smoothScalarFields%id" % db.nDim)
        vectorSmooth = eval("smoothVectorFields%id" % db.nDim)
        tensorSmooth = eval("smoothSymTensorFields%id" % db.nDim)
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
    # Allow Spheral to smoothly exit out of the controller loop.
    #--------------------------------------------------------------------------
    def stop(self):
        self._break = True
        return

    #--------------------------------------------------------------------------
    # Advance the system to the given simulation time.  The user can also
    # specify a max number of steps to take.
    #--------------------------------------------------------------------------
    def advance(self, goalTime, maxSteps=None):
        currentSteps = 0
        while (self.time() < goalTime and
               (maxSteps == None or currentSteps < maxSteps) and (self._break == False)):
            self.stepTimer.start()
            self.integrator.step(goalTime)
            if self.numHIterationsBetweenCycles > 0:
                self.iterateIdealH(maxIdealHIterations = self.numHIterationsBetweenCycles)
            self.stepTimer.stop()
            currentSteps = currentSteps + 1
            self.totalSteps = self.totalSteps + 1

            # Do the periodic work.
            self.doPeriodicWork()

        # Force the periodic work to fire at the end of an advance (except for any redistribution).
        if maxSteps != 0:
            thpt = self.redistribute
            self.redistribute = None
            self.doPeriodicWork(force=True)
            self.redistribute = thpt

        db = self.integrator.dataBase()
        bcs = self.integrator.uniqueBoundaryConditions()
        numActualGhostNodes = 0
        for bc in bcs:
            numActualGhostNodes += bc.numGhostNodes
        print "Total number of (internal, ghost, active ghost) nodes : (%i, %i, %i)" % (mpi.allreduce(db.numInternalNodes, mpi.SUM),
                                                                                        mpi.allreduce(db.numGhostNodes, mpi.SUM),
                                                                                        mpi.allreduce(numActualGhostNodes, mpi.SUM))

        # Print how much time was spent per integration cycle.
        self.stepTimer.printStatus()

        # Output any timer info
        Timer.TimerSummary()

        return

    #--------------------------------------------------------------------------
    # Step method, where the user can specify to take a given number of steps.
    #--------------------------------------------------------------------------
    def step(self, steps=1):
        self.advance(1e40, steps)

    #--------------------------------------------------------------------------
    # Do the periodic work.
    #--------------------------------------------------------------------------
    def doPeriodicWork(self, force=False):
        dt = self.lastDt()
        t1 = self.time()
        t0 = t1 - dt

        # Do any periodic cycle work, as determined by the number of steps.
        for method, frequency in self._periodicWork:
            if frequency is not None and (force or self.totalSteps % frequency == 0):
                method(self.totalSteps, t1, dt)

        # Do any periodic time work, as determined by the current time.
        for method, frequency in self._periodicTimeWork:
            if frequency is not None and (force or ((t0 // frequency) != (t1 // frequency))):
                method(self.totalSteps, t1, dt)

        return

    #--------------------------------------------------------------------------
    # Add a (method, frequency) tuple to the cyclic periodic work.
    # call during advance.
    #--------------------------------------------------------------------------
    def appendPeriodicWork(self, method, frequency):
        self._periodicWork.append((method, frequency))

    #--------------------------------------------------------------------------
    # Add a (method, frequency) tuple to the time based periodic work.
    #--------------------------------------------------------------------------
    def appendPeriodicTimeWork(self, method, frequency):
        self._periodicTimeWork.append((method, frequency))

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
    # Periodically force garbage collection.
    #--------------------------------------------------------------------------
    def garbageCollection(self, cycle, Time, dt):
##         everything = globalFrame().f_globals
##         stuff = [x for x in everything if hasattr(everything[x], "__wards__")]
##         print "SpheralController.garbageCollection: Stuff I found with wards: ", stuff
##         for x in stuff:
##             y = everything[x]
##             for z in y.__wards__:
##                 del z
##             del y.__wards__
        gc.collect()
        return

    #--------------------------------------------------------------------------
    # Periodically drop a restart file.
    # We do a garbage collection clean up pass here since sometimes temporary
    # objects are restartable and wind up dumping data here.  Doesn't really
    # hurt anything, but it can be wasteful.
    #--------------------------------------------------------------------------
    def updateRestart(self, cycle, Time, dt):
        import gc
        while gc.collect():
            pass
        self.dropRestartFile()
        return

    #--------------------------------------------------------------------------
    # Periodically redistribute the nodes between domains.
    #--------------------------------------------------------------------------
    def updateDomainDistribution(self, cycle, Time, dt):
        if self.redistribute:

            # It is *critical* that each NodeList have the same number of fields
            # registered against it on each processor, therefore we pause to
            # garbage collect here and make sure any temporaries are gone.
            import gc
            while gc.collect():
                pass

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
        file = self.restartFileConstructor(fileName, FileIOSpace.Create)
        RestartRegistrar.instance().dumpState(file)
        print "Wrote restart file in %0.2f seconds" % (time.clock() - start)

        file.close()
        del file
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
        if self.restartFileConstructor is SiloFileIO:
            fileName += ".silo"
        if not os.path.exists(fileName):
            raise RuntimeError("File %s does not exist or is inaccessible." %
                               fileName)

        # Read that sucker.
        print 'Reading from restart file', fileName
        import time
        start = time.clock()
        if self.restartFileConstructor is GzipFileIO:
            file = self.restartFileConstructor(fileName, FileIOSpace.Read)
                                               #readToMemory = True)
        else:
            file = self.restartFileConstructor(fileName, FileIOSpace.Read)
        RestartRegistrar.instance().restoreState(file)
        print "Finished: required %0.2f seconds" % (time.clock() - start)

        # Do we need to force a boundary update to create ghost nodes?
        if (self.integrator.updateBoundaryFrequency > 1 and
            self.integrator.currentCycle % self.integrator.updateBoundaryFrequency != 0):
            print "Creating ghost nodes."
            start = time.clock()
            self.integrator.setGhostNodes()
            print "Finished: required %0.2f seconds" % (time.clock() - start)

        file.close()
        del file
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

        # This is the list of boundary types that need to precede the distributed
        # boundary, since their ghost nodes need to be communicated.
        precedeDistributed = []
        if 3 in dims:
            precedeDistributed += [BoundarySpace.CylindricalBoundary,
                                   BoundarySpace.SphericalBoundary]
        for dim in dims:
            exec("""
precedeDistributed += [BoundarySpace.PeriodicBoundary%(dim)sd,
                       BoundarySpace.ConstantBoundary%(dim)sd]
""" % {"dim" : dim})

        # Check if this is a parallel process or not.
        if mpi.procs == 1:
            self.domainbc = None

        # If this is a parallel run, then we need to create a distributed
        # boundary condition and insert it into the list of boundaries for each physics
        # package.
        else:
            # exec("from SpheralModules.Spheral.BoundarySpace import NestedGridDistributedBoundary%s" % self.dim)
            # self.domainbc = eval("NestedGridDistributedBoundary%s.instance()" % self.dim)
            # from SpheralModules.Spheral.BoundarySpace import BoundingVolumeDistributedBoundary1d, \
            #                                                  BoundingVolumeDistributedBoundary2d, \
            #                                                  BoundingVolumeDistributedBoundary3d
            # self.domainbc = eval("BoundingVolumeDistributedBoundary%s.instance()" % self.dim)
            exec("from SpheralModules.Spheral.BoundarySpace import TreeDistributedBoundary%s" % self.dim)
            self.domainbc = eval("TreeDistributedBoundary%s.instance()" % self.dim)

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

                # Now put the boundaries back into the package.
                # NOTE!  We currently force the parallel boundary condition to the end of the list.
                # This is required because Boundary conditions are combinatorial -- i.e., ghost nodes
                # from prior boundary conditions can be used as controls in later.  If we put the parallel
                # boundary at the beginning of the list it's ghost state is not valid until finalizeBoundary
                # is called.
                for bc in precedeBoundaries + followBoundaries + [self.domainbc]:
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
        if mpi.procs > 1:
            from SpheralModules.Spheral import PartitionSpace
            try:
                #self.redistribute = eval("PartitionSpace.ParmetisRedistributeNodes%s(W.kernelExtent)" % self.dim)
                #self.redistribute = eval("PartitionSpace.SortAndDivideRedistributeNodes%s(W.kernelExtent)" % self.dim)
                #self.redistribute = eval("PartitionSpace.PeanoHilbertOrderRedistributeNodes%s(W.kernelExtent)" % self.dim)
                self.redistribute = eval("PartitionSpace.VoronoiRedistributeNodes%s(W.kernelExtent)" % self.dim)
            except:
                print "Warning: this appears to be a parallel run, but Controller cannot construct"
                print "         dynamic redistributer."
                pass
        return

    #---------------------------------------------------------------------------
    # Add visualization.
    #---------------------------------------------------------------------------
    def addVisualizationDumps(self, vizBaseName, vizDir, vizStep, vizTime,
                              vizFields, vizFieldLists):
        self.vizBaseName = vizBaseName
        self.vizDir = vizDir
        self.vizStep = vizStep
        self.vizTime = vizTime
        self.vizFields = vizFields
        self.vizFieldLists = vizFieldLists
        if vizStep is not None:
            self.appendPeriodicWork(self.dropViz, vizStep)
        if vizTime is not None:
            self.appendPeriodicTimeWork(self.updateViz, vizTime)
        return

    #---------------------------------------------------------------------------
    # Periodically drop viz files.
    #---------------------------------------------------------------------------
    def updateViz(self, cycle, Time, dt):
        self.dropViz(cycle, Time, dt)
        return

    #---------------------------------------------------------------------------
    # Actually does the viz dump.
    #---------------------------------------------------------------------------
    def dropViz(self,
                cycle = None,
                Time = None,
                dt = None):
        mpi.barrier()
        if self.vizMethodPoints:
            self.vizMethodPoints(self.integrator,
                                 baseFileName = self.vizBaseName,
                                 baseDirectory = os.path.join(self.vizDir, "points"),
                                 fields = self.vizFields,
                                 fieldLists = self.vizFieldLists,
                                 currentTime = self.time(),
                                 currentCycle = self.totalSteps,
                                 dumpGhosts = self.vizGhosts,
                                 dumpDerivatives = self.vizDerivs,
                                 boundaries = self.integrator.uniqueBoundaryConditions())
        if self.vizMethodCells:
            self.vizMethodCells(self.integrator,
                                baseFileName = self.vizBaseName,
                                baseDirectory = os.path.join(self.vizDir, "cells"),
                                fields = self.vizFields,
                                fieldLists = self.vizFieldLists,
                                currentTime = self.time(),
                                currentCycle = self.totalSteps,
                                dumpGhosts = self.vizGhosts,
                                dumpDerivatives = self.vizDerivs,
                                boundaries = self.integrator.uniqueBoundaryConditions())
        return

    #--------------------------------------------------------------------------
    # Iteratively set the H tensors, until the desired convergence criteria
    # are met.
    #--------------------------------------------------------------------------
    def iterateIdealH(self, 
                      maxIdealHIterations = 50,
                      idealHTolerance = 1.0e-4):
        print "SpheralController: Initializing H's..."
        db = self.integrator.dataBase()
        bcs = self.integrator.uniqueBoundaryConditions()
        if self.SPH:
            method = eval("SPHSmoothingScale%s()" % self.dim)
        else:
            method = eval("ASPHSmoothingScale%s()" % self.dim)
        iterateIdealH = eval("iterateIdealH%s" % self.dim)
        iterateIdealH(db, bcs, self.kernel, method, maxIdealHIterations, idealHTolerance, 0.0, False, False)

        return

    #---------------------------------------------------------------------------
    # Reinitialize the mass of each node such that the Voronoi mass density
    # matches the expected values.
    #---------------------------------------------------------------------------
    def voronoiInitializeMass(self):
        from generateMesh import generateLineMesh, generatePolygonalMesh, generatePolyhedralMesh
        db = self.integrator.dataBase()
        nodeLists = db.fluidNodeLists()
        boundaries = self.integrator.uniqueBoundaryConditions()
        method = eval("generate%sMesh" % {1 : "Line", 2 : "Polygonal", 3 : "Polyhedral"}[db.nDim])
        mesh, void = method(nodeLists,
                            boundaries = boundaries,
                            generateParallelConnectivity = False)
        for nodes in nodeLists:
            mass = nodes.mass()
            rho = nodes.massDensity()
            for i in xrange(nodes.numInternalNodes):
                mass[i] = rho[i]*mesh.zone(nodes, i).volume()
        return

    #---------------------------------------------------------------------------
    # Use a given pre-relaxation object (like the Voronoi hourglass control) to
    # relax a node distribution to some converged state.
    #---------------------------------------------------------------------------
    def prerelaxNodeDistribution(self,
                                 hourglass,
                                 rho,
                                 maxIterations = 100,
                                 tol = 1.0e-5):

        # What did the user pass in for rho?
        if type(rho) == type(1.0):
            rho = ConstantRho(rho)

        db = self.integrator.dataBase()
        nodeLists = db.fluidNodeLists()
        boundaries = self.integrator.uniqueBoundaryConditions()
        allpackages = self.integrator.physicsPackages()
        packages = eval("vector_of_Physics%id()" % db.nDim)
        packages.append(hourglass)

        # A helpful method for setting the density.
        def setRho():
            for nodes in nodeLists:
                pos = nodes.positions()
                massDensity = nodes.massDensity()
                for i in xrange(nodes.numInternalNodes):
                    massDensity[i] = rho(pos[i])
        setRho()

        # Iterate until we're done.
        iter = 0
        maxDisp = 2.0*tol
        while (iter < maxIterations and
               maxDisp > tol):
            iter += 1
            oldpos = db.fluidPosition
            oldpos.copyFields()
            state = eval("State%id(db, allpackages)" % db.nDim)
            derivs = eval("StateDerivatives%id(db, allpackages)" % db.nDim)
            self.integrator.initialize(state, derivs)
            state.update(derivs, 1.0, 0.0, 1.0)
            self.integrator.enforceBoundaries()
            self.integrator.applyGhostBoundaries()
            self.integrator.postStateUpdate()
            self.integrator.finalizeGhostBoundaries()
            self.integrator.finalize(0.0, 0.0, db, state, derivs)

            # Check the displacements.
            maxDisp = 0.0
            newpos = db.fluidPosition
            for (oldf, newf) in zip(oldpos, newpos):
                maxDisp = max([(old - new).magnitude() for old, new in zip(oldf.internalValues(), newf.internalValues())] + [0.0])
            maxDisp = mpi.allreduce(maxDisp, mpi.MAX)
            print " --> Iteration %i : max change = %g" % (iter, maxDisp)

        # That's about it.
        setRho()
        return
